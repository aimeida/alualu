// utils function for alu_insert 
#include "insert_utils.h"

bool alu_mate_flag( BamFileHandler * bam_fh, MapFO &fileMap ){
  seqan::BamAlignmentRecord record;
  int rID_pre = -1;
  while ( bam_fh -> fetch_a_read(record) ) { 
    if ( !QC_insert_read(record) ) continue;
    if ( bam_fh->rID_chrn.find(record.rID) == bam_fh->rID_chrn.end() ) continue;
    if ( record.rID != rID_pre) {
      if (rID_pre > -1) {
	cerr << "done with " << rID_pre << endl; 
      }
      rID_pre = record.rID;	
    }       
    //if (ia++ % 1000000 == 0) cerr << ia << " " << endl;
    int len_seq = length(record.seq);
    if ( getAlignmentLengthInRef(record) < (size_t)(len_seq - 2 * CLIP_BP)) continue; // want very good quality reads 
    if (record.rID < record.rNextId and bam_fh->rID_chrn.find(record.rNextId) != bam_fh->rID_chrn.end() ) { 
      *(fileMap[record.rNextId]) << "dif_chr " << record.qName << " " << record.rNextId << " " << record.pNext << " " 
				 << record.rID << " " << record.beginPos << " " << len_seq << endl;      
      *(fileMap[record.rID]) << "dif_chr " << record.qName << " " << record.rID << " " << record.beginPos << " " 
				 << record.rNextId << " " << record.pNext << " " << len_seq << endl;            
      continue;
    } 

    if (record.rID == record.rNextId and record.beginPos < record.pNext) { 
      if (abs(record.tLen) > DISCORDANT_LEN) {
	*(fileMap[record.rID]) << "ilong " << record.qName << " " << record.rID << " " << record.pNext << " " 
			       << record.rID << " " << record.beginPos << " " << len_seq << endl;      
	*(fileMap[record.rID]) << "ilong " << record.qName << " " << record.rID << " " << record.beginPos << " " 
			       << record.rID << " " << record.pNext << " " << len_seq << endl;      
      }    
    }
  }
  return 0;
}

// only support 4 types for now 
string parse_alu_type(string alu_name){
  assert ( !alu_name.empty() );
  if ( alu_name.substr(0,4) == "AluY") return "AluY";
  if ( alu_name.substr(0,4) == "AluS") return "AluSx";
  if ( alu_name.substr(0,4) == "AluJ") {
    if ( alu_name.substr(0,5) == "AluJo") return "AluJo";
    if ( alu_name.substr(0,5) == "AluJb") return "AluJb";
    return "AluJo";
  }
  return "AluY";
}

/** move read *S*M to left until can't move */
bool clipRight_move_left(seqan::CharString & read_seq, seqan::CharString & ref_fa, list <int> & cigar_cnts, int refBegin, int & clipPos, int & align_len) {
  int i = 1;  // i = 0 is the last match by BWA
  while ( *cigar_cnts.begin() - i >= 0 and 
	  ref_fa[clipPos - refBegin - i] == read_seq[*cigar_cnts.begin() - i] )
    i++;
  if ( *cigar_cnts.begin() - (i-1) < CLIP_BP )
    return false;
  clipPos -= (i-1);  
  align_len = length(read_seq) - *cigar_cnts.begin() + (i-1);
  //if (movePos) cout << "move "  << movePos << endl;
  return true;  
}

/** move read *M*S to right until can't move */
bool clipLeft_move_right(seqan::CharString & read_seq, seqan::CharString & ref_fa, list <int> & cigar_cnts, int refBegin, int & clipPos, int & align_len){
  align_len = length(read_seq) - cigar_cnts.back();
  int i = 0;
  while (align_len + i < (int) length(read_seq) and ref_fa[clipPos - refBegin + i] == read_seq[align_len + i]) i++;
  if ( cigar_cnts.back() - i < CLIP_BP)
    return false;
  clipPos += i;
  align_len += i;
  return true;
}

bool global_align_insert(int hasRCFlag, seqan::CharString seq_read, seqan::CharString seq_ref, int &score, int cutEnd, float th_score, bool verbose){
  score = 0;
  if (verbose) {
    cout << "#1 " << seq_read << endl;
    cout << "#2 " << seq_ref << endl;
  }

  size_t align_len = length(seq_read);
  assert ( align_len == length(seq_ref) );
  if ( align_len < CLIP_BP) return false;
  seqan::Score<int> scoringScheme(1, -2, -2, -2); // match, mismatch, gap extend, gap open
  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align,0), seq_read); // 2,3, true means free gap
  assignSource(row(align,1), seq_ref);   // 1,4
  if (!hasRCFlag)
    score = globalAlignment(align, scoringScheme, seqan::AlignConfig<false, false, true, true>()); 
  else
    score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, false, false>());   
  if (score >= round(th_score * align_len)  )
    return true;
  
  if (verbose)  cout << align << endl;
  
  if ( (int)align_len <= CLIP_BP + cutEnd) return false;
  if (!hasRCFlag){
    assignSource(row(align,0), infix(seq_read, 0, align_len - cutEnd )); 
    assignSource(row(align,1), infix(seq_ref, 0, align_len - cutEnd ));  
  } else {
    assignSource(row(align,0), infix(seq_read, cutEnd, align_len)); 
    assignSource(row(align,1), infix(seq_ref, cutEnd, align_len));      
  }
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<false, false, false, false>()); 
  
  if (verbose)  cout << align << endl;

  return score >= round(th_score * (align_len - cutEnd)) ;
}
