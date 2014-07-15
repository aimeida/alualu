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

// get inferred split positions
bool read_first2col(string fn, vector < pair<int, int> > & insert_pos, bool has_header) {
  insert_pos.clear();
  ifstream fin(fn.c_str());
  if (!fin) return false;
  stringstream ss;
  string line, tmpv;
  if (has_header)  getline(fin, line); 
  int beginPos, endPos;
  int ni = 0;
  std::set < pair <int, int> > lefts_rights;
  while (getline(fin, line)) {
    ni ++;
    ss.clear(); ss.str( line );  
    ss >> beginPos >> endPos;
    if ( abs (beginPos - endPos ) < MAX_POS_DIF ) {
      lefts_rights.insert( get_valid_pair(beginPos, endPos) ); 
    } else if (beginPos and endPos) {
      lefts_rights.insert( make_pair(beginPos, beginPos) ); 
      lefts_rights.insert( make_pair(endPos, endPos) ); 
    }   
  } 
  fin.close();
  if (lefts_rights.empty()) return false;
  std::set < pair <int, int> >::iterator si = lefts_rights.begin();
  int beginPre = (*si).first;
  int endPre =  (*si).second;
  si++;
  for ( ; si != lefts_rights.end(); si++) {
    beginPos = (*si).first;
    endPos = (*si).second;
    if ( abs(beginPre - beginPos) < MAX_POS_DIF ) {
      beginPre = (beginPre + beginPos) / 2 ;
      endPre = (endPre + endPos) / 2 ;
    } else {
      insert_pos.push_back( make_pair(beginPre, endPre) );
      beginPre = beginPos;
      endPre = endPos;
    }
  }
  insert_pos.push_back(make_pair(beginPre, endPre) );
  return !insert_pos.empty(); 
}

bool parseline_del_tmp0(string line0, string & output_line, map <int, EmpiricalPdf *> & pdf_rg, int estimatedAluLen, string line1){
  float *log10_gp = new float[3];
  stringstream ss, ss_out;
  string chrn, pos, posr, token;
  int idx, insert_len;
  int clipCnt = 0, unknowCnt = 0, midCnt = 0;  
  vector < pair <int, int> > unknowInfo;
  if ( !line0.empty() ) { 
    ss.str(line0);
    ss >> chrn >> pos >> posr >> clipCnt >> unknowCnt ;
    for (int i = 0; i < unknowCnt; i++) {
      getline(ss, token, ':');
      seqan::lexicalCast2(idx, token);
      getline(ss, token, ' ');
      seqan::lexicalCast2(insert_len, token);      
      unknowInfo.push_back( make_pair(idx, insert_len) );
    }
  }
  int _aluCnt, _clipCnt, _ukCnt;
  if ( !line1.empty() ) { 
    ss.clear(); ss.str(line1);
    ss >> chrn >> pos >> _aluCnt >> _clipCnt >> _ukCnt;
    midCnt = _aluCnt + _clipCnt;
    unknowCnt += _ukCnt;
    for (int i = 0; i < _ukCnt; i++) {
      getline(ss, token, ':');
      seqan::lexicalCast2(idx, token);
      getline(ss, token, ' ');
      seqan::lexicalCast2(insert_len, token);
      unknowInfo.push_back( make_pair(idx, insert_len) );
    }   
  }

  float prob_ub = pow(10, -LOG10_RATIO_UB);
  float prob_known = (midCnt+clipCnt)/(float)(midCnt + clipCnt + unknowCnt);
  for (int i = 0; i < 3; i++) log10_gp[i] = 0;
  if (prob_known) {
    log10_gp[0] = clipCnt * log10 ( prob_known * prob_ub ) + midCnt * log10 ( prob_known * (1 - prob_ub) ); // no deletion (insertion)
    log10_gp[1] = (midCnt + clipCnt) * log10 (prob_known * 0.5) ; 
    log10_gp[2] = midCnt * log10 ( prob_known * prob_ub ) + clipCnt * log10 ( prob_known * (1 - prob_ub) );
  }

  for (vector < pair <int, int> >::iterator ui = unknowInfo.begin(); ui != unknowInfo.end(); ui++ ) {
    int idx = (*ui).first;
    int insert_len = (*ui).second;
    float p_y = pdf_rg[idx]->pdf_obs(insert_len + estimatedAluLen);   // no deletion (insertion)
    float p_z = pdf_rg[idx]->pdf_obs(insert_len);
    //float freq0 = 0.67;  // high FP ratio      
    float freq0 = ( midCnt + 1 )/(float)(midCnt + clipCnt + 2); // 1 and 2 are psudo count
    log10_gp[0] += log10 (p_y * (1 - prob_known));
    log10_gp[1] += log10 ((freq0 * p_y + (1 - freq0) * p_z) * (1 - prob_known) ) ;
    log10_gp[2] += log10 (p_z * (1 - prob_known));
  }
  
  bool use_this_line = false;
  if ( !p11_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    ss_out << chrn << " " << pos << " " << estimatedAluLen << " " << midCnt << " " << clipCnt << " " << unknowCnt ;
    float *gp = new float[3];
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
    ss_out << " " << setprecision(6) << gp[2] << " " << gp[1] << " " << gp[0];  // NB: switch 00 and 11, unlike alu_delete
    delete gp;    
    output_line = ss_out.str();
    use_this_line = true;
  }
  delete log10_gp;
  return use_this_line;
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
  if (global_align_insert( hasFlagRC(record), read_clip_fa, ref_fa, score, ALIGN_END_CUT, 0.7) ) {
    fout << header << left_right << " " << adj_clipPos << " " << record.qName << " " << clipPos << " " <<  get_cigar(record) << endl;
    return true;
  }
  return false;
}

int align_clip_to_LongConsRef(string shortSeq, string longSeq, int & refBegin, int & refEnd,  int clipLen){
  const int max_diff_cnt = 6;
  refBegin = 0, refEnd = 0; // initial 
  int shortLen = shortSeq.size();
  TAlign align;
  seqan::Score<int> scoringScheme(1, -1, -2, -3); 
  resize(rows(align), 2);
  assignSource(row(align,0), shortSeq);  // 2,3
  assignSource(row(align,1), longSeq);   // 1,4, free gap at end
  globalAlignment(align, scoringScheme, seqan::AlignConfig<true, false, false, true>());
  int l0 = toViewPosition(row(align, 1), 0);
  int l1 = toViewPosition(row(align, 1), longSeq.size());
  int s0 = toViewPosition(row(align, 0), 0);
  int s1 = toViewPosition(row(align, 0), shortLen);
  if ( (s0 < l0) or (s1 > l1) or ( s1 - s0 - shortLen >= max_diff_cnt ) )
    return 0;
  TRow &row0 = row(align,0);
  TRowIterator it0 = begin(row0);
  TRow &row1 = row(align,1);
  TRowIterator it1 = begin(row1);
  int i = 0, dif = 0;
  while ( i++ < s0 ) {  it0++; it1++; }
  while ( i++ <= s1 ) {
    if ( (*it0) != (*it1) ) dif++;     ////if(isGap(it1))
    it0++; it1++;
  }
  if ( dif <= max_diff_cnt) {
    if ( clipLen > 0 )  refEnd = s0 + clipLen;
    else if ( clipLen < 0 )  refBegin = s1 + clipLen;
  }    
  return 0;
}

void align_clip_to_consRef(string shortSeq, string longSeq, int & refBegin, int & refEnd,  int clipLen){
  const float dif_th = 0.2;
  int shortLen = shortSeq.size();
  TAlign align;
  seqan::Score<int> scoringScheme(1, -1, -2, -3); 
  resize(rows(align), 2);
  assignSource(row(align,0), shortSeq);  // 2,3
  assignSource(row(align,1), longSeq);   // 1,4, free gap at end
  globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>());
  int l0 = toViewPosition(row(align, 1), 0);
  int l1 = toViewPosition(row(align, 1), longSeq.size());
  int s0 = toViewPosition(row(align, 0), 0);
  int s1 = toViewPosition(row(align, 0), shortLen);
  int align_len = min(s1, l1) - max(s0, l0);
  int align_len_shortSeq = shortSeq.size();

  if ( abs(align_len_shortSeq - align_len) <= dif_th * align_len_shortSeq 
       and align_len >= CLIP_BP ){
    if ( clipLen > 0 ) // ??? 
      refEnd = s0 + clipLen;  // might > longSeq.size() 
    else if ( clipLen < 0 ) 
      refBegin = s1 + clipLen;  // might be negative
  }
}

bool align_alu_to_consRef(const string & shortSeq, const string & longSeq, float dif_th, string loginfo) {
  TAlign align;
  seqan::Score<int> scoringScheme(1, -1, -2, -3); 
  resize(rows(align), 2);
  assignSource(row(align,0), shortSeq);  // 2,3
  assignSource(row(align,1), longSeq);   // 1,4, free gap at end
  globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>());
  int s0 = toViewPosition(row(align, 0), 0);
  int s1 = toViewPosition(row(align, 0), shortSeq.size());
  int l0 = toViewPosition(row(align, 1), 0); 
  int l1 = toViewPosition(row(align, 1), longSeq.size());  
  int align_len = min(s1, l1) - max(s0, l0);
  int align_len_shortSeq = shortSeq.size();
  if ( l0 > s0 )
    align_len_shortSeq -= (l0 - s0);
  else if ( s1 > l1 )
    align_len_shortSeq -= (s1 - l1);
  bool align_pass =  align_len_shortSeq >= CLIP_BP and abs(align_len_shortSeq - align_len) <= dif_th * align_len_shortSeq ;

  if ( !align_pass and !loginfo.empty())  {
    cout << loginfo << ": align_len " << align_len << endl;
    cout << shortSeq << endl;
    cout << longSeq << endl;
    cout << align << endl;
  }
    
  return align_pass;
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
