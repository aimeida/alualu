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

void normalize_llh(float *loggp, float th){
  th = - abs(th);
  float maxv = max(max(loggp[0], loggp[1]), loggp[2]);
  float minv = min(min(loggp[0], loggp[1]), loggp[2]);
  int id_max = -1, id_mid = -1, id_min = -1;
  for (int i = 0; i < 3; i ++ ) {
    if ( abs(loggp[i] - maxv) < 0.00001 and id_max < 0) id_max = i;
    else if ( abs(loggp[i] - minv) < 0.00001 and id_min < 0) id_min = i;
    else id_mid = i;
  }
  //if ( ! (id_max > -1 and id_mid > -1 and id_min > -1) ) cout << loggp[0] << " " <<  loggp[1] << " " << loggp[2] << endl;
  assert(id_max > -1 and id_mid > -1 and id_min > -1);
  loggp[id_mid] += abs(loggp[id_max]);
  loggp[id_min] += abs(loggp[id_max]);
  loggp[id_max] = 0;
  if ( loggp[id_min] < th ) {
    loggp[id_min] = th;
    float th_mid = loggp[id_max] - (loggp[id_max] - loggp[id_mid]) * 0.8;
    if ( loggp[id_mid] < th_mid ) loggp[id_mid] = th_mid;
  }
}

int parseline_ins(string line0, ostream & fout, map <int, EmpiricalPdf *> & pdf_rg, float logPE, int estimatedAluLen, int errCode, bool test_print, float *test_gp) {
  //errCode1: 0/0 or .. ==> 0/1, reduce ph0 and useLen
  //errCode2: 1/1 ==> 0/1
  //errCode3: 0/1 or .. ==> 1/1, reduce ph0 ?? (0.2)
  const float ratioMax = 6.; 
  const float maxLh = 40.;
  float ph0 = 0.3;
  if (errCode == 1) ph0 = 1. / (1. + ratioMax); 
  if (test_print) cout << "ph0 is " << ph0 << endl;

  float *log10_gp = new float[3];
  float *log10_gpu = new float[3];
  float *gp = new float[3];
  for (int i = 0; i < 3; i++) { log10_gp[i] = 0; log10_gpu[i] = 0; gp[i] = 0; }
  stringstream ss;
  string chrn, insertMid, debugInfo, token;
  int idx, insert_len, midCnt, clipCnt, unknowCnt;
  vector < pair <int, int> > unknowInfo;
  ss.str(line0); 
  ss >> chrn >> insertMid >> debugInfo >> midCnt;
  string exact_left, exact_right;
  split_by_sep(debugInfo, exact_left, exact_right, ',');  
  //bool both_side = (exact_left != "0") and (exact_right != "0");
  if ( exact_left  == "0" ) exact_left = exact_right;
  if ( exact_right == "0" ) exact_right = exact_left;
  if ( midCnt < 0 ) {
    fout << chrn << " " << exact_left << " " << debugInfo << " " << midCnt << endl;
    return 0;
  } 
  ss >> clipCnt >> unknowCnt;
  int cntCnt = midCnt + clipCnt;
  int covCnt = cntCnt + unknowCnt;
  if ( covCnt < 3 or covCnt > 1000 ) return 0; // considered as missing  
  if ( min(midCnt, clipCnt) >= MID_COV_CNT ) {
    //if ( midCnt * ratioMax < clipCnt) 
    fout << chrn << " " << exact_left << " " << debugInfo << " " << midCnt << " " << clipCnt << " " << unknowCnt 
	 << " 0 1 0 " << estimatedAluLen << endl;
    return 0;
  }
  logPE = - abs(logPE);
  float down_weight = abs(logPE) - 1.; // use less 
  float maxLhCnt = min( (float) ( abs(logPE) * 0.8 * cntCnt), maxLh);
  float maxLhLen = min( (float) ( down_weight * 0.8 * unknowCnt), maxLh);  // generally use less 
  if (cntCnt > 0 ) {
    float logPM = log10 ( 1 - pow(10, logPE) );
    log10_gp[0] = clipCnt * logPE + midCnt * logPM;
    log10_gp[1] = midCnt * log10 (ph0) + clipCnt * log10 (1 - ph0) + (midCnt + clipCnt) * logPM  ;
    log10_gp[2] = midCnt * logPE + clipCnt * logPM;
    if (test_print)   cout << maxLhCnt << " ## " <<  log10_gp[2] << " " << log10_gp[1] << " " << log10_gp[0] << " , ";
    normalize_llh(log10_gp, maxLhCnt); 
    if (test_print)   cout << log10_gp[2] << " " << log10_gp[1] << " " << log10_gp[0] << endl;
  }
  if ( unknowCnt > 0) {
    if (errCode == 4 ) estimatedAluLen = 180; // reduce to minume
    for (int i = 0; i < unknowCnt; i++) {
      getline(ss, token, ':');
      seqan::lexicalCast2(idx, token);
      if ( pdf_rg.find(idx) == pdf_rg.end()) idx = 0; // for debug, if idx not exists
      getline(ss, token, ' ');
      seqan::lexicalCast2(insert_len, token);      
      float p_y, p_z;
      //if (test_print) cout << insert_len << " " << pdf_rg[idx]->pdf_obs(insert_len + estimatedAluLen) << " " << pdf_rg[idx]->pdf_obs(insert_len) << endl;   
      pdf_rg[idx]->ratio_obs(insert_len + estimatedAluLen, insert_len, down_weight, p_y, p_z);
      log10_gpu[0] += log10 (p_y);
      log10_gpu[1] += log10 (ph0 * p_y + (1 - ph0) * p_z) ;
      log10_gpu[2] += log10 (p_z);
    }
    if (test_print)   cout <<  maxLhLen << " ## " <<  log10_gpu[2] << " " << log10_gpu[1] << " " << log10_gpu[0] << " , ";
    normalize_llh(log10_gpu, maxLhLen);
    if (test_print)   cout << log10_gpu[2] << " " << log10_gpu[1] << " " << log10_gpu[0] << endl;
  }
  
  for (int i = 0; i < 3; i++) log10_gp[i] += log10_gpu[i];
  if (test_print)   cout << "logAll " <<  log10_gp[2] << " " << log10_gp[1] << " " << log10_gp[0] << endl;
  
  if (midCnt >= 3 and log10_gp[2] > max (log10_gp[1], log10_gp[0]) ) { // call 0/0 when midCnt obs                     
    if (errCode == 0 ) return 1;
    if (errCode == -1) cerr << "##errCode1\n" ;
    else cerr << "##errCode1 " << line0 << endl;
  }
  if (clipCnt >= 3 and log10_gp[0] > max (log10_gp[1], log10_gp[2]) ) {
    ///////if (errCode == 0) return 2;
    if (errCode == -1) cerr << "##errCode2\n" ;
    else cerr << "##errCode2 " << line0 << endl;
  }
  if ( clipCnt <= 2 and midCnt > 3 * max(clipCnt,2) and log10_gp[0] < max (log10_gp[1], log10_gp[2]) ) {
    ///////if (errCode == 0 ) return 3;
    if (errCode == -1) cerr << "##errCode3\n" ;
    else cerr << "##errCode3 " << line0 << endl;
  }
  if ( log10_gp[0] > max(log10_gp[1], log10_gp[2]) and log10_gpu[2] > max(log10_gpu[0], log10_gpu[1]) ) {
    if (errCode == 0 ) return 4;
    if (errCode == -1) cerr << "##errCode4\n" ;
    else cerr  << "##errCode4 " << line0 << endl;
  }

  if ( !p11_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1  
    fout << chrn << " " << exact_left << " " << debugInfo << " " << midCnt << " " << clipCnt << " " << unknowCnt
	 << " " << setprecision(6) << gp[2] << " " << gp[1] << " " << gp[0] << " " << estimatedAluLen << endl;  // NB: switch 00 and 11
  } else 
    fout << chrn << " " << exact_left << " " << debugInfo << " " << -(midCnt + clipCnt + unknowCnt) << endl;

  if ( test_print ) cout << gp[2] << " " << gp[1] << " " << gp[0] << endl;
  if ( errCode == -1 )  for (int i = 0; i < 3; i++ ) test_gp[i] = gp[i];

  delete log10_gp;
  delete log10_gpu;
  delete gp;    
  return 0;
}

// used for debugging 
int parseline_cnt(string line0) {
  stringstream ss;
  string chrn, insertMid, debugInfo, token;
  int midCnt, clipCnt, unknowCnt;
  vector < pair <int, int> > unknowInfo;
  ss.str(line0); 
  ss >> chrn >> insertMid >> debugInfo >> midCnt;
  string exact_left, exact_right;
  split_by_sep(debugInfo, exact_left, exact_right, ',');  
  if ( exact_left  == "0" ) exact_left = exact_right;
  if ( exact_right == "0" ) exact_right = exact_left;
  if ( midCnt < 0 ) return 0;
  ss >> clipCnt >> unknowCnt;
  int c1 = midCnt + clipCnt + unknowCnt;
  if ( c1 > 300 ) 
    cout << line0 << endl;
  if (midCnt == 0 or clipCnt == 0 )
    cout << c1 << " # " << midCnt << " " << clipCnt << endl;
  else 
    cout << c1 << " & " << (midCnt+0.1)/(clipCnt + 0.1) << endl;
  return 0;
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

bool align_alu_cons(seqan::CharString &ref_fa, seqan::CharString alucons, float & sim_rate,float sim_th){
  TAlign align;
  seqan::Score<int> scoringScheme(0, -1, -1, -2); 
  resize(rows(align), 2);
  assignSource(row(align,0), ref_fa);  // 2,3
  assignSource(row(align,1), alucons);   // 1,4, free gap at end
  globalAlignment(align, scoringScheme, seqan::AlignConfig<true, false, false, true>());
  int align_start = max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  int align_end = min(toViewPosition(row(align, 0), length(ref_fa)), toViewPosition(row(align, 1), length(alucons)));
  sim_rate = 0;
  int align_len = align_end - align_start;
  if ( align_len <= CLIP_BP or align_len <= sim_th * length(ref_fa))
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
  ///if (sim_rate >= sim_th) cout << "ok " << align << endl;
  return  sim_rate >= sim_th;
}

string align_alu_cons_call(seqan::CharString & ref_fa, AluconsHandler *alucons_fh, float & sim_rate, float sim_th){
  for ( vector<string>::iterator si = (alucons_fh->seq_names).begin(); si != (alucons_fh->seq_names).end(); si++) {
    alucons_fh->update_seq_name(*si);
    for (int k = 1; k <= 4; k++) 
      if (align_alu_cons(ref_fa, alucons_fh->fetch_alucons(k), sim_rate, sim_th))
	return *si;
  }  
  return "";
}

bool covered_reads(BamFileHandler * bam_fh, string chrn, int p1, int p2, int minCnt){
  if (! bam_fh->jump_to_region(chrn, p1, p2) ) return false;
  seqan::BamAlignmentRecord record;
  int cnt = 0;
  while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, p1, p2, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_insert_read(record)) continue; 
      if (record.rID != record.rNextId) continue;
      if (abs(record.tLen) > DISCORDANT_LEN ) {
	cnt++;
      } else {
	int r1 = min ( record.beginPos, record.beginPos + abs(record.tLen) );
	int r2 = max ( record.beginPos, record.beginPos + abs(record.tLen) );
	if ( min(p2, r2) - max (r1, p1) > 0 )
	  cnt++;
      }
      if ( cnt >= minCnt) return true;
  }
  return cnt >= minCnt;
}
