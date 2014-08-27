#include "delete_utils.h"

string tread_toStr(T_READ td){
  if ( td == unknow_read ) return "unknow_read";
  if ( td == mid_read ) return "mid_read";
  if ( td == clip_read ) return "clip_read";
  return "uesless_read?";
}

bool get_align_pos(int aluBegin, int aluEnd, int beginPos, int endPos, int &ref_a, int &ref_b, int &read_a, int &read_b, seqan::BamAlignmentRecord &record){
  unsigned nl = length(record.cigar) - 1; // ignore complicated alignment between S and M
  int len_read = length(record.seq);
  int len_clip_toAlign = 0;
  if ( abs(beginPos - aluEnd) <= BOUNDARY_OFFSET and record.cigar[0].operation == 'S') {  //eg. S41M60, has_soft_first 
    len_clip_toAlign = record.cigar[0].count;
    ref_a = aluBegin - len_clip_toAlign;
    ref_b = aluBegin;
    read_a = 0;
    read_b = len_clip_toAlign;
  } else if (abs(endPos - aluBegin) <= BOUNDARY_OFFSET and record.cigar[nl].operation == 'S' ) { // eg. M27S93, has_soft_last
    len_clip_toAlign = record.cigar[nl].count;
    ref_a = aluEnd;
    ref_b = aluEnd + len_clip_toAlign;
    read_a = len_read - len_clip_toAlign;
    read_b = len_read;
  }
  if ( !len_clip_toAlign ) return false;
  /* beginning of reads is bad quality. we cut a few bp before aligning
     eg: 10bp cut 1bp, 110bp(90bp) cut 10bp;*/
  int cut_bp = len_read > 110 ? round( 0.11 * len_clip_toAlign - 0.13 ) : round( 0.1 * len_clip_toAlign ); 
  if ( hasFlagRC(record) and read_a == 0) read_a += cut_bp;    
  if ( (!hasFlagRC(record)) and read_b == len_read) read_b -= cut_bp;
  return read_a < read_b;
}

bool split_global_align(seqan::CharString &fa_seq, seqan::BamAlignmentRecord &record, int read_a, int read_b){
  // no need to use seed or chain, since don't know exact position to start mapping
  seqan::Score<int> scoringScheme(1, -3, -2, -5); // match, mismatch, gap extend, gap open
  TAlign align;
  int align_len, score;
  int min_align_len = max(CLIP_BP, (int) (0.85 * (read_b - read_a - BOUNDARY_OFFSET) ) ); 
  resize(rows(align), 2);
  assignSource(row(align,0), infix(record.seq, read_a, read_b) );  // 2,3
  assignSource(row(align,1), fa_seq);   // 1,4

  // free gap both ends
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>()); 
  align_len = min(toViewPosition(row(align, 0), read_b - read_a), toViewPosition(row(align, 1), length(fa_seq)))
    - max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  //cout << clippedBeginPosition(row(align, 0)) << " " << clippedBeginPosition(row(align, 1)) << " " <<  endl;
  //cout << "score1 " << score << " " << align_len << " " << min_align_len << " " << endl;
  if ( align_len >= min_align_len and score >= round( 0.75 * align_len - 2 ) ) return true; 
  
  // AlignConfig 2=true: free end gaps for record.seq
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<false, true, true, false>()); 
  align_len = min(toViewPosition(row(align, 0), read_b - read_a), toViewPosition(row(align, 1), length(fa_seq)))
    - max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  return ( align_len >= min_align_len and score >= round( 0.75 * align_len - 2 ));
}

T_READ classify_read(seqan::BamAlignmentRecord & record, int aluBegin, int aluEnd, FastaFileHandler *fasta_fh, bool debug){
  bool read_is_left = left_read(record);
  int beginPos = record.beginPos;
  int endPos = record.beginPos + getAlignmentLengthInRef(record);    
  int pair_begin = read_is_left ? beginPos : record.pNext;
  int pair_end = pair_begin + abs(record.tLen);
  if ( (has_soft_first(record, CLIP_BP) and abs(beginPos - aluEnd) <= BOUNDARY_OFFSET ) or 
       ( has_soft_last(record, CLIP_BP) and abs(endPos - aluBegin) <= BOUNDARY_OFFSET ) ) {    
    int ref_a, ref_b, read_a, read_b;
    if ( get_align_pos(aluBegin, aluEnd, beginPos, endPos, ref_a, ref_b, read_a, read_b, record)) {
      seqan::CharString fa_seq;
      fasta_fh->fetch_fasta_upper(ref_a - BOUNDARY_OFFSET, ref_b + BOUNDARY_OFFSET, fa_seq);          
      if (split_global_align(fa_seq, record, read_a, read_b)) 
	return clip_read;
    }
  }  
  if ( debug )
    cout << "debug# " <<  beginPos <<  " " << endPos << " " << pair_begin << " " << pair_end << endl;
  
  // only consider as mid_read if very certain, otherwise classify as unknown read
  if ( beginPos < aluEnd - BOUNDARY_OFFSET and endPos > aluBegin + BOUNDARY_OFFSET and count_non_match(record) <= 5)       
    return mid_read;    
  if ( pair_begin < aluBegin + BOUNDARY_OFFSET and pair_end > aluEnd - BOUNDARY_OFFSET )
    return unknow_read;  
  return useless_read;  // alu_flank is too large, we have a lot reads not useful 
}

bool combine_pns_vcf(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> & chrns, map <string, std::set<int> > & chrn_aluBegin, float llh_th, string ref_name) {
  vector <string>::iterator ci, pi;
  //seqan::VcfStream vcfout("-", seqan::VcfStream::WRITE);
  seqan::VcfStream vcfout(seqan::toCString(f_out), seqan::VcfStream::WRITE);
  for ( ci = chrns.begin(); ci != chrns.end(); ci++)
    appendValue(vcfout.header.sequenceNames, *ci);
  for ( pi = pns.begin(); pi != pns.end(); pi++) 
    appendValue(vcfout.header.sampleNames, *pi);
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("fileformat", "VCFv4.1"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("fileDate", "201402"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("reference", ref_name));  // eg: hg19
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=Integer, Description=\"Genotype\">"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FORMAT", "<ID=PL,Number=3,Type=Integer, Description=\"Phred-scaled likelihoods for genotypes\">"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FILTER", "<ID=highCov, Description=\"Region coverage too high\">"));
  seqan::VcfRecord record;    
  record.ref = ".";
  record.alt = "1";
  record.info = ".";
  record.format = "GT:PL";

  ifstream fin;
  stringstream ss;
  string line, chrn;
  int aluBegin, aluEnd, cii, flag;
  int midCnt, clipCnt, unknowCnt;
  float p0, p1, p2;      
  string phred_str_00 = "0,255,255";
  string phred_str_missing = ".,.,.";
  
  map < pair<int, int>, map<string, string> > pos_pnProb;
  map < pair<int, int>, map<string, string> > pos_pnGeno;
  for ( cii = 0, ci = chrns.begin(); ci != chrns.end(); ci++, cii++) {
    pos_pnProb.clear();
    pos_pnGeno.clear();
    for ( pi = pns.begin(); pi != pns.end(); pi++) {
      fin.open( (path0 + *pi + f_in_suffix).c_str() );
      assert(fin);
      getline(fin, line); 
      flag = 1;
      while ( getline(fin, line) ) {
	ss.clear(); ss.str( line );
	ss >> chrn;
	if (chrn != *ci) {
	  if (flag) continue;
	  else break;
	}
	flag = 0; 
	ss >>  aluBegin >> aluEnd >> midCnt;
	if ( midCnt < 0 ) {  // evidence of G00
	  pos_pnProb[ make_pair(aluBegin, aluEnd) ][*pi] = phred_str_00;
	} else {
	  ss >> clipCnt >> unknowCnt >> p0 >> p1 >> p2 ;
	  string phred_scaled_str = phred_scaled(p0, p1, p2);
	  pos_pnProb[ make_pair(aluBegin, aluEnd) ][*pi] = phred_scaled_str;
	  if (p1 > p0 or p2 > p0) {
	    pos_pnGeno[ make_pair(aluBegin, aluEnd) ][*pi] = ( p1 > p2) ? "1" : "2"; // only G01 and G11
	  }
	}
      }
      fin.close();
    }
    
    // Write out the records.    
    map < pair<int, int>, map<string, string> >::iterator pi3, pj3;
    map<string, string>::iterator pi2, pj2;
    for (pi3 = pos_pnProb.begin(); pi3 != pos_pnProb.end(); pi3++) {
      pj3 = pos_pnGeno.find ( pi3->first );
      if ( pj3 == pos_pnGeno.end() )
	continue;
      record.beginPos = (pi3->first).first;
      record.rID = cii;   // change ???
      record.id = int_to_string((pi3->first).second);      
      record.qual = 0;  // change ??
      for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
	string ginfo;
	pi2 = pi3->second.find(*pi);
	if ( pi2 == pi3->second.end() ) { 
	  ginfo = ".:" + phred_str_missing ;
	} else {
	  pj2 = pj3->second.find(*pi);
	  if ( pj2 == pj3->second.end() )
	    ginfo = "0:" + pi2->second;
	  else 
	    ginfo = pj2->second + ":" + pi2->second;
	}
	appendValue(record.genotypeInfos, ginfo);
      }
      if ( chrn_aluBegin[*ci].find( (pi3->first).first ) == chrn_aluBegin[*ci].end() )
	record.filter = "PASS";
      else
	record.filter = "highCov";
      writeRecord(vcfout, record);  
      clear(record.genotypeInfos);
    } 
    cout << *ci << " size " << pos_pnProb.size() << endl;
  }  // chrn finished
  clear(record);
  seqan::close(vcfout);
  return true;
}

