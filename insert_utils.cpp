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
    } else {
      if (beginPos) lefts_rights.insert( make_pair(beginPos, beginPos) ); 
      if (endPos) lefts_rights.insert( make_pair(endPos, endPos) ); 
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

bool parseline_del_tmp0(string line0, string & output_line, map <int, EmpiricalPdf *> & pdf_rg, float logPE, int estimatedAluLen) {
  float *log10_gp = new float[3];
  stringstream ss, ss_out;
  string chrn, insertMid, debugInfo, token;
  int idx, insert_len, midCnt, clipCnt, unknowCnt;
  vector < pair <int, int> > unknowInfo;
  ss.str(line0); 
  ss >> chrn >> insertMid >> debugInfo >> midCnt >> clipCnt >> unknowCnt;

  float ph0 = 0.5;
  logPE = - abs(logPE);
  float logPM = log10 ( 1 - pow(10, logPE) );

  float *gp = new float[3];
  for (int i = 0; i < 3; i++) log10_gp[i] = 0;
  if (clipCnt + midCnt > 0 ) {
    log10_gp[0] = clipCnt * logPE + midCnt * logPM;
    log10_gp[1] = midCnt * log10 (ph0) + clipCnt * log10 (1 - ph0) + (midCnt + clipCnt) * logPM  ;
    log10_gp[2] = midCnt * logPE + clipCnt * logPM;
  }

  for (int i = 0; i < unknowCnt; i++) {
    getline(ss, token, ':');
    seqan::lexicalCast2(idx, token);
    getline(ss, token, ' ');
    seqan::lexicalCast2(insert_len, token);      
    float p_y, p_z;
    pdf_rg[idx]->ratio_obs(insert_len + estimatedAluLen, insert_len, abs(logPE) - 0.5, p_y, p_z);
    log10_gp[0] += log10 (p_y);
    log10_gp[1] += log10 (ph0 * p_y + (1 - ph0) * p_z) ;
    log10_gp[2] += log10 (p_z);
  }
  
  bool use_this_line = false;
  if ( !p11_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    string exact_left, exact_right;
    split_by_sep(debugInfo, exact_left, exact_right, ',');
    if (exact_left=="0") exact_left = exact_right;
    if (exact_right=="0") exact_right = exact_left;
    ss_out << chrn << " " << exact_left << " " << exact_right << " " << estimatedAluLen << " " << midCnt << " " << clipCnt << " " << unknowCnt ;
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
    ss_out << " " << setprecision(6) << gp[2] << " " << gp[1] << " " << gp[0];  // NB: switch 00 and 11, unlike alu_delete
    use_this_line = true;
    output_line = ss_out.str();
  }

  delete gp;    
  delete log10_gp;
  return use_this_line;
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
