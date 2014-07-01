// utils function for alu_insert 
#include "insert_utils.h"

bool AlumateINFO::sort_pos2(const AlumateINFO* a, const AlumateINFO* b){
  return a-> pos2 < b-> pos2;
}

void AlumateINFO::delete_list(list <AlumateINFO *> & alumate_list){
  for ( list <AlumateINFO *> ::iterator ai = alumate_list.begin(); ai != alumate_list.end(); ai++ ) 
    delete *ai;
  alumate_list.clear();
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


bool read_match_clipLeft(string & line, int clipLeft, string & pn, string & qName) {
  stringstream ss;
  string tmp1, tmp2, sleft_right;
  int clipPos;
  ss.str( line );
  ss >> pn >> tmp1 >> tmp2 >> sleft_right >> clipPos >> qName;
  int match_offset =  ( sleft_right[0] == 'L') ? CLIP_BP_LEFT : CLIP_BP_RIGHT ;
  return  abs(clipPos - clipLeft ) <= match_offset;
}

bool read_match_clipLeft(string & line, int clipLeft, string & pn, string & qName, string & cigar, string & seq) {
  stringstream ss;
  string tmp1, tmp2, tmp3, sleft_right, tLen;
  int clipPos;
  ss.str( line );
  ss >> pn >> tmp1 >> tmp2 >> sleft_right >> clipPos >> qName >> tmp3 >> cigar >> tLen >> seq;
  int match_offset =  ( sleft_right[0] == 'L') ? CLIP_BP_LEFT : CLIP_BP_RIGHT ;
  return  abs(clipPos - clipLeft ) <= match_offset;
}

bool align_clip_to_ref(char left_right, int adj_clipPos,  int clipPos, int align_len, seqan::BamAlignmentRecord &record, FastaFileHandler *fasta_fh, ofstream &fout, string  header) {
  const int ref_plus_bp = 10; // allow small indels
  int score;
  seqan::CharString ref_fa;
  seqan::CharString read_clip_fa;
  int read_len = length(record.seq);
  if (left_right == 'R') {
    fasta_fh->fetch_fasta_upper(adj_clipPos, adj_clipPos + align_len + ref_plus_bp, ref_fa);    
    read_clip_fa = infix(record.seq, read_len - align_len, read_len);
  } else if (left_right == 'L') {
    fasta_fh->fetch_fasta_upper(adj_clipPos - align_len - ref_plus_bp, adj_clipPos, ref_fa);
    read_clip_fa = infix(record.seq, 0, align_len);
  }
//      debug_print_read(record, cout);
//      cout << "#R " << get_cigar(record) << " " << adj_clipPos << " " << adj_clipPos + align_len << endl;
//      cout << record.seq << endl;
//      fasta_fh->fetch_fasta_upper(adj_clipPos - 10, adj_clipPos + align_len + 10, ref_fa);    
//      cout << infix(ref_fa, 0, 10) << " " <<  infix(ref_fa, 10, 10 + align_len) << " " << infix(ref_fa, 10 + align_len, align_len + 20) << endl;
//      global_align_insert( hasFlagRC(record), infix(record.seq, read_len - align_len, read_len), ref_fa, score, ALIGN_END_CUT, 0.7, true);

  if (global_align_insert( hasFlagRC(record), read_clip_fa, ref_fa, score, ALIGN_END_CUT, 0.7) ) {
    fout << header << left_right << " " << adj_clipPos << " " << record.qName << " " << clipPos << " " <<  get_cigar(record) << endl;
    return true;
  }
  return false;
}

bool global_align_insert(const int hasRCFlag, seqan::CharString & seq_read, seqan::CharString & seq_ref, int &score, int cutEnd, float th_score, bool verbose){
  score = 0;
  if (verbose) {
    cout << "verbose #1 " << seq_read << endl;
    cout << "verbose #2 " << seq_ref << endl;
  }
  size_t align_len = length(seq_read);
  if ( align_len < CLIP_BP) return false;
  seqan::Score<int> scoringScheme(1, -2, -2, -2); // match, mismatch, gap extend, gap open
  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align,0), seq_read); // 2,3, true means free gap
  assignSource(row(align,1), seq_ref);   // 1,4
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>()); 
  if (score >= round(th_score * align_len)  )
    return true;
  // otherwise cut the end and realign
  if (verbose)  cout << align << endl;  
  if ( (int)align_len <= CLIP_BP + cutEnd) return false;
  if (!hasRCFlag){
    assignSource(row(align,0), infix(seq_read, 0, align_len - cutEnd )); 
    assignSource(row(align,1), infix(seq_ref, 0, align_len - cutEnd ));  
  } else {
    assignSource(row(align,0), infix(seq_read, cutEnd, align_len)); 
    assignSource(row(align,1), infix(seq_ref, cutEnd, align_len));      
  }
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>());   
  if (verbose)  cout << align << endl;
  return score >= round(th_score * (align_len - cutEnd)) ;
}

bool align_alu_cons(string &ref_fa, seqan::CharString alucons, float & sim_rate,float sim_th){
  TAlign align;
  seqan::Score<int> scoringScheme(0, -1, -1, -2); 
  resize(rows(align), 2);
  assignSource(row(align,0), ref_fa);  // 2,3
  assignSource(row(align,1), alucons);   // 1,4, free gap at end
  globalAlignment(align, scoringScheme, seqan::AlignConfig<true, false, false, true>());
  int align_start = max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  int align_end = min(toViewPosition(row(align, 0), ref_fa.size()), toViewPosition(row(align, 1), length(alucons)));
  sim_rate = 0;
  int align_len = align_end - align_start;
  if ( align_len <= CLIP_BP or align_len <= sim_th * ref_fa.size() )
    return false;
  TRow &row0 = row(align,0);
  TRowIterator it0 = begin(row0);
  TRow &row1 = row(align,1);
  TRowIterator it1 = begin(row1);
  int i = 0, dif = 0;
  while ( i++ < align_start ) {  it0++; it1++; }
  while ( i++ <= align_end) {
    if ( (*it0) != (*it1) ) dif++;     ////if(isGap(it1))
    it0++; it1++;
  }
  sim_rate = 1 - dif / (float) align_len;
  return  sim_rate >= sim_th;
}

int align_alu_cons_call(string & ref_fa, AluconsHandler *alucons_fh, float & sim_rate, float sim_th){
  for (int k = 1; k <= 4; k++) {
    if (align_alu_cons(ref_fa, alucons_fh->fetch_alucons(k), sim_rate, sim_th))
      return k;
  }
  return 0;
}

int align_clip_to_consRef(string shortSeq, string longSeq, int & refBegin, int & refEnd,  int clipLen){
  const int max_diff_cnt = 6;
  refBegin = 0, refEnd = 0; // initial 
  int shortLen = shortSeq.size();
  TAlign align;
  seqan::Score<int> scoringScheme(0, -2, -2, -4); 
  resize(rows(align), 2);
  assignSource(row(align,0), shortSeq);  // 2,3
  assignSource(row(align,1), longSeq);   // 1,4, free gap at end
  globalAlignment(align, scoringScheme, seqan::AlignConfig<true, false, false, true>());
  int long0 = toViewPosition(row(align, 1), 0);
  int long1 = toViewPosition(row(align, 1), longSeq.size());
  int short0 = toViewPosition(row(align, 0), 0);
  int short1 = toViewPosition(row(align, 0), shortLen);
  if ( (short0 < long0) or (short1 > long1) or ( short1 - short0 - shortLen >= max_diff_cnt ) )
    return 0;

  TRow &row0 = row(align,0);
  TRowIterator it0 = begin(row0);
  TRow &row1 = row(align,1);
  TRowIterator it1 = begin(row1);
  int i = 0, dif = 0;
  while ( i++ < short0 ) {  it0++; it1++; }
  while ( i++ <= short1 ) {
    if ( (*it0) != (*it1) ) dif++;     ////if(isGap(it1))
    it0++; it1++;
  }
  if ( dif <= max_diff_cnt) {
    if ( clipLen > 0 )  
      refEnd = short0 + clipLen;
    if ( clipLen < 0 )  
      refBegin = short1 + clipLen;
    ///cout << refBegin << " " << refEnd << " " << clipLen << " " << short0 << " " << short1 << endl;
  }    
  return 0;
}

void filter_outlier_pn(string path_input, string fn_suffix, map<int, string> &ID_pn, string chrn, string file_pn_used_output, float percentage_pn_used) {
  string line;
  int ni;
  map < string, int > pn_lineCnt;
  set <int> lineCnt;  
  for (map<int, string>::iterator pi = ID_pn.begin(); pi != ID_pn.end(); pi++) {
    string file_st = path_input + pi->second + "." + fn_suffix + "." + chrn;
    ifstream fin(file_st.c_str());
    ni = 0;
    while (fin >> line) ni++;
    fin.close();
    lineCnt.insert(ni);
    pn_lineCnt[pi->second] = ni;
  }  
  set <int>::iterator li = lineCnt.begin();
  int cnt_th = 0;
  ni = 0;
  while ( ni++ < percentage_pn_used * lineCnt.size())
    if ( ++li != lineCnt.end() ) cnt_th = *li;
  ofstream fout(file_pn_used_output.c_str());
  for (map < string, int >::iterator pi = pn_lineCnt.begin(); pi != pn_lineCnt.end(); pi++)
    if ( pi->second <= cnt_th) fout << pi->first << endl;
  fout.close();
}
