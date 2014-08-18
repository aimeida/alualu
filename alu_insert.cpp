#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include "delete_utils.h"

struct PAR_DELETE0 {
  int CLIP_Q;
};

struct REALIGNS {
  int rid, pos_key, aBegin, aEnd;
  string qName;
};

void alu_mate_flag( BamFileHandler * bam_fh, string fn_output, string &header, map <string, int> & rg_to_idx){
  ofstream fout (fn_output.c_str() );  
  fout << header ;
  seqan::BamAlignmentRecord record;
  while ( bam_fh -> fetch_a_read(record) ) { 
    if ( ! QC_insert_read(record) ) continue;
    if ( ! key_exists(bam_fh->rID_chrn, record.rID) ) continue;
    int len_seq = length(record.seq);
    if ( count_non_match(record) >= CLIP_BP )  continue; // want very good quality reads 
    ///// allow pair to be mapped to decoy,  key_exists(bam_fh->rID_chrn, record.rNextId)
    if ( record.rID < record.rNextId or
	 (record.rID == record.rNextId and record.beginPos < record.pNext and abs(record.tLen) > DISCORDANT_LEN) ) {
      fout << record.qName << " " << record.rID << " " << record.rNextId <<  " " << record.beginPos << " " << record.pNext
	   << " " << hasFlagRC(record) << " " << hasFlagNextRC(record) << " " << len_seq << " " << get_rgIdx(rg_to_idx, record) << endl;      
    }
  }
  fout.close();
}

int parse_line3( string &line, int & pos) {
  int lr_num, rr_num;
  stringstream ss;
  ss.str(line);
  ss >> pos >> lr_num >> rr_num;  
  return lr_num + rr_num;
}

void parse_line1(string &line, int pos, int offset, vector <int> & clipls, vector <int> & cliprs)  {
  stringstream ss;
  ss.str(line);
  string pn, tmpv;
  int clipLen, p;
  ss >> pn >> tmpv;
  if ( tmpv.substr(0,3) != "Alu" ) {
    ss >> clipLen >> p;
    if ( abs( p - pos) <= offset ) {
      if (clipLen < 0 ) clipls.push_back(p);
      else  cliprs.push_back(p);
    }
  }
}

void get_scan_pos(string fn, list <int> & pos_to_scan, int combine_bp){ 
  pos_to_scan.clear();
  stringstream ss;
  string line, qname;
  int this_rID, this_pos, len_read, previous_end;  
  bool this_isRC;
  ifstream fin( fn.c_str());
  assert(fin); // init, 700 * 30 /100 = 210 
  getline(fin, line);   
  while (pos_to_scan.empty()) {
    getline(fin, line);
    ss.clear(); ss.str(line);
    ss >> qname >> this_rID >> this_pos >> this_isRC >> len_read;
    if ( this_pos > SCAN_WIN_LEN ) pos_to_scan.push_back(this_pos);
  }
  while (getline(fin, line)) {
    previous_end = this_pos + len_read;    
    ss.clear(); ss.str(line);
    ss >> qname >> this_rID >> this_pos >> this_isRC >> len_read;
    if (this_pos - previous_end > 3 * combine_bp) pos_to_scan.push_back(previous_end + combine_bp);
    if (this_pos - pos_to_scan.back() > 2 * combine_bp) pos_to_scan.push_back(this_pos - combine_bp);
  }
}                                     

// need smarter way to combine positions ==> minimum number of regions to cover all reads info
void join_location(string file1, string file2, int pos_dif, int max_region_len){ 
  ifstream fin( file1.c_str());
  assert(fin);    
  string line;
  getline(fin, line); // skip header
  getline(fin, line);
  int pos, pos_pre;
  int nowCount = parse_line3(line, pos_pre);
  //cout << "nowC " << nowCount << endl;
  int maxCount = nowCount;
  vector<int>  maxCount_pos;
  maxCount_pos.push_back(pos_pre);
  stringstream sout; 
  while (getline(fin, line)) {
    nowCount = parse_line3(line, pos);
    if ( pos - pos_pre <= pos_dif and pos - (*maxCount_pos.begin()) <= max_region_len) {// same block
      if ( nowCount > maxCount) {
	maxCount = nowCount;
	maxCount_pos.clear();
	maxCount_pos.push_back(pos);
      } else if ( nowCount == maxCount) {
	maxCount_pos.push_back(pos);
	}
    } else { // create new block      
      sout << maxCount << " " << maxCount_pos.front() << " " << maxCount_pos.back() <<  endl;
      maxCount_pos.clear();
      maxCount_pos.push_back(pos);      
      maxCount = nowCount;
    }
    pos_pre = pos;
  }
  fin.close();
  ofstream fout( file2.c_str());
  fout << sout.str();
  fout.close();
}

void alumate_counts_filter(string & path0, string pn, string fn_suffix_input, string fn_suffix_output, string fn_suffix_output2, vector<string>  &chrns, int th_cnt){
  int lr_num, rr_num, this_rID, this_pos, len_read;
  string line, qname;
  bool this_isRC;
  stringstream ss;
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    string chrn = *ci;
    string f_input = path0 + chrn + "/" + pn + fn_suffix_input;
    string f_output = path0 + chrn + "/" + pn + fn_suffix_output;
    string f_output2 = path0 + chrn + "/" + pn + fn_suffix_output2;
    list <int> pos_to_scan;    
    get_scan_pos(f_input, pos_to_scan, 8); // no big diff when change 8 to 20
    // cout << f_input << " " << pos_to_scan.size() << endl;
    list< READ_INFO *> lr_reads;    
    ifstream fin ( f_input.c_str());  
    assert(fin); 
    ofstream fout( f_output.c_str());
    fout << "check_pos lr_num rr_num\n";	
    getline(fin, line); // skip header
    getline(fin, line); 
    ss.clear(); ss.str(line);
    ss >> qname >> this_rID >> this_pos >> this_isRC >> len_read;
    lr_reads.push_back(new READ_INFO(this_pos, len_read, !this_isRC, "")); 
    bool readable = true;
    while ( readable and pos_to_scan.size()) {
      int check_pos = pos_to_scan.front();
      pos_to_scan.pop_front();
      while (readable and (lr_reads.empty() or lr_reads.back()-> endPos < check_pos + SCAN_WIN_LEN) ) {
	if (getline(fin, line)) {
	  ss.clear(); ss.str(line);
	  ss >> qname >> this_rID >> this_pos >> this_isRC >> len_read;
	  if ( len_read > 50)
	    lr_reads.push_back(new READ_INFO(this_pos, len_read, !this_isRC, ""));
	} else 
	  readable = false;
      }
      if (lr_reads.empty()) continue;      
      lr_num = 0; rr_num = 0;      	
      list< READ_INFO *>::iterator ri = lr_reads.begin();
      while (ri != lr_reads.end() ) {
	if ( (*ri)->beginPos < check_pos - SCAN_WIN_LEN ) {
	  delete *ri;
	  ri = lr_reads.erase(ri); // <==> lr_reads.erase(ri++)
	  continue;
	}
	if ( (*ri)->endPos < check_pos and (*ri)->should_be_left) lr_num++;
	else if ( (*ri)->beginPos > check_pos and !(*ri)->should_be_left ) rr_num++;
	ri++;
	if ( (*ri)->endPos >= check_pos + SCAN_WIN_LEN ) break;
      }
      if (lr_num + rr_num >= th_cnt) fout << check_pos << " " << lr_num  << " " << rr_num << endl;	
    }
    fin.close();
    fout.close();
    join_location(f_output, f_output2, MAX_POS_DIF, MAX_LEN_REGION);    		       
    cout << "done " << f_output << endl;
  }
}

void filter_location_rep(string path0, string pn, string file1_suffix, string file2_suffix, vector<string> &chrns, RepMaskPos *repmaskPos){
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    string chrn = *ci;
    ifstream fin( (path0 + chrn + "/" + pn + file1_suffix).c_str());
    assert(fin);
    ofstream fout( (path0 + chrn + "/" + pn + file2_suffix).c_str());
    string line;
    stringstream ss;
    int readCnt, pos_left, pos_right;
    vector<int>::iterator bi = repmaskPos->beginP[chrn].begin();
    vector<int>::iterator ei = repmaskPos->endP[chrn].begin();
    int bei = 0;
    int be_size = repmaskPos->beginP[chrn].size();
    while ( getline(fin, line) ) { 
      ss.clear(); ss.str( line );
      ss >> readCnt >> pos_left >> pos_right;
      if (pos_right <= 0) continue;
      while ( pos_left >= (*ei) and bei < be_size) { bi++; ei++; bei++; }
      if ( min(*ei, pos_right) - max(*bi, pos_left) < 0) // NB: < 0, not <= 0
	fout << line << endl;  // not in Alu
    }  
    fin.close();
    fout.close();
  }
}

int read_sort_by_col(string fn, int coln, bool has_header, list< IntString> &rows_list) {
  string line, tmpv;
  int pos;
  stringstream ss;
  ifstream fin( fn.c_str());
  assert(fin); 
  rows_list.clear();
  if (has_header) getline(fin, line);
  size_t rown = 0;
  while (getline(fin, line)) {
    rown++;
    ss.clear(); ss.str( line );
    for (int i = 0; i < coln-1; i++) ss >> tmpv;
    ss >> pos;
    rows_list.push_back( make_pair(pos, line) );
  }
  fin.close();
  rows_list.sort(compare_IntString);
  if (rows_list.size() != rown ) cerr << "##### ERROR #### " << fn << endl;
  return rows_list.size();
}

// qname A_chr B_chr A_beginPos B_beginPos A_isRC B_isRC len_read
void read_match_rid(string fn, int col_val, list< AlumateINFO *> & alumate_list){
  string line, qname;
  int rid1, rid2, pos1, pos2, len_read, rgIdx;
  bool rc1, rc2;
  stringstream ss;
  ifstream fin( fn.c_str());
  assert(fin); 
  getline(fin, line); 
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> qname >> rid1 >> rid2 >> pos1 >> pos2 >> rc1 >> rc2 >> len_read >> rgIdx; 
    if (rid1 == col_val)
      alumate_list.push_back(new AlumateINFO(qname, len_read, rid2, rid1, pos2, pos1, rgIdx, rc2, rc1));
    if (rid2 == col_val) 
      alumate_list.push_back(new AlumateINFO(qname, len_read, rid1, rid2, pos1, pos2, rgIdx, rc1, rc2));
  }
}

void keep_alu_mate(BamFileHandler * bam_fh, int alu_rid, string file1, MapFO & fileMap, string file_alu){
  const int ALU_MIN_LEN = 50;  // min overlap in alu region
  list< AlumateINFO *> alumate_list;
  read_match_rid(file1, alu_rid, alumate_list);  
  alumate_list.sort(AlumateINFO::sort_pos2);  
  cout << "alu file " << file_alu << endl;
  AluRefPos *alurefpos = new AluRefPos(file_alu, 100, 100);   // combine nearby positions, considered as Alu             
  int bei = 0, n_alu = 0;
  int be_size = alurefpos->db_size - 1; // -1 ! otherwise seg fault
  int left_pos;
  for (list< AlumateINFO *>::iterator ri = alumate_list.begin(); ri != alumate_list.end();  ri++) {
    assert ( (*ri)->rid2 == alu_rid );
    left_pos = (*ri)->pos2;   
    while ( ( left_pos >= alurefpos->get_endP() ) and bei < be_size) { alurefpos->nextdb(); bei++; }
    if (  min(  alurefpos->get_endP(), left_pos + (*ri)->len_read ) - max( alurefpos->get_beginP(), left_pos) > ALU_MIN_LEN and
	  key_exists(bam_fh->rID_chrn, (*ri)->rid1) and n_alu++) 
      *(fileMap[(*ri)->rid1]) << (*ri)->qname << " " << (*ri)->rid1 << " " << (*ri)->pos1 << " " << (*ri)->RC1 << " " << (*ri)->len_read 
			      << " " << (*ri)->rid2 << " "  << left_pos << " " << (*ri)->RC2 << " " << (*ri)->rgIdx << " " << alurefpos->get_type() << endl ;          
  }
  ////// about 18% err counts are mapped to Alu region(normal chrns)
  cout << file1 << " err counts: " << alumate_list.size() << ", alu mate: " << n_alu << endl; 
  AlumateINFO::delete_list(alumate_list);
  delete alurefpos;      
}

bool write_tmpFile_all_loci( vector<string> &fn_inputs, vector <string> & fn_idx, string fn_output) {
  stringstream ss;
  string line;
  ofstream fout;
  vector<string>::iterator di;
  if (!fn_idx.empty()) di = fn_idx.begin();
  //cout << "####write to " << fn_output << " " << fn_inputs.size() << " " << fn_idx.size() << endl;
  fout.open(fn_output.c_str());
  for (vector<string>::iterator fi = fn_inputs.begin(); fi != fn_inputs.end(); fi ++ ) {
    ifstream fin( (*fi).c_str());
    if (!fin) {
      cerr << "#### !!! " << *fi << " does not exists! skip it\n";
      continue;
    }
    if ( fn_idx.empty() )  { while (getline(fin, line)) fout << line << endl; }
    else  { while (getline(fin, line)) fout << line << " " << *di << endl; di++;}
    fin.close();
  }
  fout.close();
  // sorting 
  list< IntString> rows_list;
  read_sort_by_col(fn_output, 2, false, rows_list);      
  fout.open(fn_output.c_str());
  for (list< IntString>::iterator ri = rows_list.begin(); ri!=rows_list.end(); ri++) 
    fout << (*ri).second << endl;
  fout.close();
  rows_list.clear();
  return true;
}

void join_all_loci(string f_in, string f_out, int block_dist, int block_len) {
  string line, pn;
  int cnt, sum_cnt, pa, pb, pa_block, pb_block;
  set <string> pns;
  ifstream fin(f_in.c_str());
  assert(fin);
  ofstream fout(f_out.c_str());
  stringstream ss;
  getline(fin, line);
  ss.str(line);
  ss >> cnt >> pa_block >> pb_block;
  while (ss >> pn)  pns.insert(pn);
  sum_cnt = cnt;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> cnt >> pa >> pb;
    if ( pa <= pb_block + block_dist and // same block
	 ( pa <= pa_block + block_len or pa <= pb_block) ) {
      pb_block = min( pa_block + block_len, max(pb_block, pb) );
      sum_cnt += cnt;
    } else {  // new block
      if (pb_block > 0) {
	fout << sum_cnt << " " << pa_block << " " << pb_block;
	for (set<string>::iterator si = pns.begin(); si != pns.end(); si++ ) fout << " " << *si;
	fout << endl;
      }
      pns.clear();
      pa_block = pa;
      pb_block = pb;
      sum_cnt = cnt;
    }
    while (ss >> pn)  pns.insert(pn);
  }
  fout << sum_cnt << " " << pa_block << " " << pb_block;
  for (set<string>::iterator si = pns.begin(); si != pns.end(); si++ ) fout << " " << *si;
  fout << endl;
  fin.close();
  fout.close();
}

void combine_fn(vector<string> & fn_inputs, string fn_output, vector<string> & fn_idx){
  write_tmpFile_all_loci(fn_inputs, fn_idx, fn_output + ".tmp");
  join_all_loci(fn_output + ".tmp", fn_output, MAX_POS_DIF, 2 * MAX_LEN_REGION);
  system( ("rm " + fn_output + ".tmp").c_str() );  
}

bool combine_pn_pos(string chrn, map <int, string> &id_pn_map, string tmpn,  string output_prefix, int group_size, string path0) {
  string fn_prefix = output_prefix + chrn ; // some random names
  // first iteration
  map <string, vector<string> > fnOut_fnInput;   
  map <string, vector<string> >::iterator fni;
  map <string, vector<string> > fnOut_pns;
  string fnOut,  fnInput_prefix, pn;
  int i, j, n_size, n_group;
  n_size = id_pn_map.size();
  n_group = (int)(n_size / group_size);
  for (i = 0; i < n_group; i++) {
    fnOut = (n_size <= group_size) ? fn_prefix : (fn_prefix + "." + int_to_string(i) + ".tmp0");
    for ( j = 0; j < group_size; j++) {
      pn = id_pn_map[ i*group_size+j ];
      fnOut_fnInput[fnOut].push_back(path0 + chrn  + "/" + pn + tmpn);
      fnOut_pns[fnOut].push_back(pn);
    }
  }
  fnOut = (n_size <= group_size) ? fn_prefix : (fn_prefix + "." + int_to_string(n_group) + ".tmp0");
  for ( j = n_group * group_size; j < n_size; j++) {
    pn = id_pn_map[ j ];
    fnOut_fnInput[fnOut].push_back(path0 + chrn + "/" + pn + tmpn);
    fnOut_pns[fnOut].push_back(pn);
  }
  for (fni = fnOut_fnInput.begin(); fni != fnOut_fnInput.end(); fni ++) 
    combine_fn(fni->second, fni->first, fnOut_pns[fni->first]);		  
  fnOut_pns.clear();
  if ( fnOut_fnInput.size() <= 1)  return true;

  vector<string> empty_vec;
  int fn_idx = 0;    
  while ( true ) {
      n_size =  fnOut_fnInput.size();
      n_group = (int)( n_size / group_size);
      fnOut_fnInput.clear();
      for (i = 0; i < n_group; i++) { 
	fnOut = (n_size <= group_size) ? fn_prefix  : (fn_prefix + "."+ int_to_string(i) + ".tmp"+ int_to_string(fn_idx+1) );
	for ( j = 0; j < group_size; j++) fnOut_fnInput[fnOut].push_back( fn_prefix + "." + int_to_string(i*group_size+j) + ".tmp"+ int_to_string(fn_idx) );
      }
      fnOut = (n_size <= group_size) ? fn_prefix : (fn_prefix+ "." + int_to_string(n_group) + ".tmp"+ int_to_string(fn_idx+1) );
      for ( j = n_group * group_size; j < n_size; j++)   fnOut_fnInput[fnOut].push_back( fn_prefix + "." + int_to_string(j) +  ".tmp" + int_to_string(fn_idx) );	
      for (fni = fnOut_fnInput.begin(); fni != fnOut_fnInput.end(); fni ++) 
	combine_fn(fni->second, fni->first, empty_vec);
      fn_idx++;
      if ( fnOut_fnInput.size() <= 1) 	break;
  }
  return true;
}

bool clipreads_at_insertPos(string pn, string chrn, BamFileHandler *bam_fh, FastaFileHandler *fasta_fh,  string fin_pos, string fout_reads, int CLIP_Q) {
  const int MAX_REGION_LEN = 400;
  ifstream fin( (fin_pos).c_str());
  if (!fin) return false;
  string line, numAluRead, _pn;
  ofstream fout(fout_reads.c_str());
  fout << "pn region_begin region_end left_right adj_clipPos qName clipPos cigar\n" ;
  stringstream ss, ss_header;
  int region_begin, region_end, refBegin, refEnd;
  int region_begin0, region_end0;
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> numAluRead >> region_begin >> region_end;
    int _tmpwin = ( MAX_REGION_LEN - region_end + region_begin) / 2;
    region_begin0 = region_begin - _tmpwin;
    region_end0 = region_end + _tmpwin;
    bool pn_has_alumate = false;
    while ( ss >> _pn ) 
      if ( _pn == pn) { pn_has_alumate=true; break;}
    if (!pn_has_alumate) continue;
    refBegin = region_begin - SCAN_WIN_LEN;
    refEnd = region_end + SCAN_WIN_LEN;
    if (! bam_fh->jump_to_region(chrn, refBegin, refEnd) ) continue;
    ss_header.clear(); ss_header.str("");
    ss_header << pn << " " << region_begin << " " << region_end << " " ;    
    seqan::BamAlignmentRecord record;
    while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, refBegin, refEnd, record);
      if (read_status == "stop" ) break;
      // allow one pair is clipped, one pair is alu_read. More reads considered
      if (read_status == "skip" or !QC_insert_read(record)) continue; 
      int beginPos = record.beginPos;
      int endPos = beginPos + getAlignmentLengthInRef(record);
      int adj_clipPos, align_len;
      seqan::CharString ref_fa;
      list <char> cigar_opts;
      list <int> cigar_cnts;
      
      int _cliplen;
      seqan::CharString _clipSeq;            
      
      bool s1 = has_soft_last(record, CLIP_BP);
      bool s2 = has_soft_first(record, CLIP_BP); 
      if (s1 and s2) continue;
      if (!s1 and !s2) continue;
      if (!aluclip_RC_match(record) ) continue;
      if ( s1 and endPos >= region_begin0 and endPos < region_end0 ) {
	adj_clipPos = endPos;
	if (length(ref_fa) == 0) {
          fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
          parse_cigar(get_cigar(record), cigar_opts, cigar_cnts);
        }
	if (clipLeft_move_right( record.seq, ref_fa, cigar_cnts, refBegin, adj_clipPos, align_len) and 
	    trim_clip_soft_last(record, _cliplen, _clipSeq, CLIP_Q) )
	  fout << ss_header.str() << "L " << adj_clipPos << " " << record.qName << " " << endPos << " " <<  get_cigar(record) << endl;	
      } 
      else if ( s2 and beginPos >= region_begin0 and beginPos < region_end0 ) {
	adj_clipPos = beginPos;
	if (length(ref_fa) == 0) {
          fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
          parse_cigar(get_cigar(record), cigar_opts, cigar_cnts);
        }
	if (clipRight_move_left( record.seq, ref_fa, cigar_cnts, refBegin, adj_clipPos, align_len) and 
	    trim_clip_soft_first(record, _cliplen, _clipSeq, CLIP_Q))
	  fout << ss_header.str() << "R " << adj_clipPos << " " << record.qName << " " << beginPos << " " <<  get_cigar(record) << endl;
      }
    }
  }
  fin.close();
  fout.close();
  return true;
}

void filter_noise ( string fn, int r, int minCov) {  // round by resolution
  ///cout << round_by_resolution(98, 10) << endl; 
  ifstream fin( fn.c_str());
  assert(fin);
  string header, line, tmpv;
  vector < pair<int, string> > pr_vec;
  map < int, int > pr_cnt;
  getline(fin, header);   
  while (getline(fin, line)) {
    stringstream ss;
    int pos, posr;
    ss.clear(); ss.str(line);
    for (int i = 0; i < 4; i++) ss >> tmpv;
    ss >> pos;
    posr = round_by_resolution(pos,r);
    if ( pr_cnt.find(posr) != pr_cnt.end() ) {
      pr_cnt[posr] += 1;
      pr_vec.push_back( make_pair(posr, line));      
    } else if ( pr_cnt.find(posr-1) != pr_cnt.end() ) {
      pr_cnt[posr-1] += 1;
      pr_vec.push_back( make_pair(posr-1, line));      
    } else if ( pr_cnt.find(posr+1) != pr_cnt.end() ) {
      pr_cnt[posr+1] += 1;
      pr_vec.push_back( make_pair(posr+1, line));      
    } else {
      pr_cnt[posr] = 1;
      pr_vec.push_back( make_pair(posr, line));      
    }
  }
  fin.close();
  ofstream fout( fn.c_str());
  fout << header << endl;
  for (vector < pair<int, string> >::iterator pi = pr_vec.begin(); pi != pr_vec.end(); pi++ ) 
    if (pr_cnt[(*pi).first] >= minCov)  // at least 2 reads !!
      fout << (*pi).second << endl;
  fout.close();
}

// for each insertion region, combine reads from different pns 
void combine_clipreads_by_pos(std::set<string> &pns_used, string path0, string path1) {
  map< pair<string, string>, vector<string> > regionPos_lines;
  for (std::set<string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) {
    string file_input_clip = path0 + *pi;
    ifstream fin(file_input_clip.c_str());
    if (!fin) {
      cout << "error " << file_input_clip <<  " missing \n";
      continue;
    }
    stringstream ss;
    string line, pn, region_begin, region_end;
    getline(fin, line);  // read header 
    while ( getline(fin, line) ) {
      ss.clear(); ss.str(line);
      ss >> pn >> region_begin >> region_end;
      regionPos_lines[make_pair(region_begin, region_end)].push_back(line); 
    }
    fin.close();
  } 
  map< pair<string, string>, vector<string> >::iterator ri;
  for (ri = regionPos_lines.begin(); ri != regionPos_lines.end(); ri ++) {
    string file_output_clip = path1 + (ri->first).first + "_" + (ri->first).second;	
    ofstream fout(file_output_clip.c_str() );
    for (vector<string>::iterator ii = (ri->second).begin(); ii != (ri->second).end(); ii++) 
      fout << *ii << endl;
    fout.close();
  }
}


string pos_to_fn( pair<int, int> ps) {
  return int_to_string(ps.first) + "_" + int_to_string(ps.second);
}

bool nonempty_files_sorted(string & path1, std::set < pair<int, int> > & regionPos_set) {
  DIR *d;
  struct dirent *dir;
  d = opendir(path1.c_str());
  string fn;
  int regionBegin, regionEnd;
  string s_regionBegin, s_regionEnd;
  if (d) {
    while ( (dir = readdir(d))!=NULL )  {
      fn = dir->d_name;      
      if (fn[0] == '.' or check_file_size( path1 + fn)<=0 ) continue;
      split_by_sep(fn, s_regionBegin, s_regionEnd, '_');
      seqan::lexicalCast2(regionBegin, s_regionBegin);
      seqan::lexicalCast2(regionEnd, s_regionEnd);
      regionPos_set.insert( make_pair(regionBegin, regionEnd));
    }
    closedir(d);
    return true;
  }
  return false;
}

void addcnt_if_pos_match( map <int, float> & clipLeft_cnt, int pos, int match_offset, float inc_val) {
  int match_key = 0;
  for ( map <int, float>::iterator mi = clipLeft_cnt.begin(); mi != clipLeft_cnt.end(); mi ++) {
    if ( abs(mi->first - pos) <= match_offset ) {
      match_offset = abs(mi->first - pos);
      match_key = mi->first;
    }
  }
  if ( match_key == 0) {
    addKey(clipLeft_cnt, pos, inc_val);
  } else {
    if (match_key > pos) { // previously add clipRight
      float cnt = clipLeft_cnt[match_key] + inc_val;
      clipLeft_cnt.erase(match_key);
      clipLeft_cnt[pos] = cnt;
      // cout << "offset " << match_key << " " << pos << endl;
    } else {
      addKey(clipLeft_cnt, match_key, inc_val);
    }
  }
}

void region_pos_vote( string fn, int & cl1, int & cl2, int & cr1, int & cr2 ) {
  const float MIN_VOTE_FREQ = 0.4;
  const int MIN_PN_CLIP_CNT = 2; // otherwise ignore this pn 
  cl1 = 0; cl2 = 0; cr1 = 0; cr2 = 0;
  ifstream fin(fn.c_str());
  assert(fin);
  stringstream ss;
  string line, pn;
  char left_right;
  int region_begin, region_end, clipPos;
  map <string, vector <pair<char, int> > > pn_splitori_pos;
  while ( getline(fin, line)) {  
    ss.clear(); ss.str( line );
    ss >> pn >> region_begin >> region_end >> left_right >> clipPos;
    pn_splitori_pos[pn].push_back( make_pair( left_right, clipPos) );
  }
  fin.close();  
  vector < int > cls, crs;
  map <string, vector <pair<char, int> > >::iterator psi;
  vector <pair<char, int> >::iterator si; 
  for ( psi = pn_splitori_pos.begin(); psi != pn_splitori_pos.end(); psi++) {
    int pl1 = 0, pl2 = 0, pr1 = 0, pr2 = 0;
    vector<int> pls, prs;    
    for (si = (psi->second).begin(); si != (psi->second).end(); si ++ ) {
      if ( (*si).first == 'L') pls.push_back( (*si).second ); 
      else prs.push_back( (*si).second ); 
    }
    int fl1, fl2;
    if ( (int)pls.size() >= MIN_PN_CLIP_CNT) 
      major_two_keys(pls, pl1, pl2, fl1, fl2, CLIP_BP_LEFT, MIN_VOTE_FREQ - 0.1);  // each pn 2 votes 
    
    int fr1, fr2;
    if ( (int)prs.size() >= MIN_PN_CLIP_CNT) 
      major_two_keys(prs, pr1, pr2, fr1, fr2, CLIP_BP_LEFT, MIN_VOTE_FREQ - 0.1);     
    
    if ( fl1 <= 1 and fr1 <= 1 and abs(pl1 - pr1) > MAX_POS_DIF ) {
      pl1 = 0; pr1 = 0;
    }
    if (pl1 > 0) {
      cls.push_back(pl1);
      if (pl2 > 0) cls.push_back(pl2);
      else cls.push_back(pl1);
    }
    if (pr1 > 0) {
      crs.push_back(pr1);
      if (pr2 > 0) crs.push_back(pr2);
      else crs.push_back(pr1);
    }
  }  
  
  int fc1, fc2;
  major_two_keys(cls, cl1, cl2, fc1, fc2, CLIP_BP_LEFT, MIN_VOTE_FREQ);
  ////if (region_begin == 26354188) cout << "cl1 " << cl1 << " " << fc1 << endl;
  if ( cl1 and fc1 <= 2 and fc1 / (float) cls.size() < 0.99 ) cl1 = 0;
  if ( cl2 and fc2 <= 2) cl2 = 0;
  
  major_two_keys(crs, cr1, cr2, fc1, fc2, CLIP_BP_LEFT, MIN_VOTE_FREQ);
  if ( cr1 and fc1 <= 2 and fc1 / (float) crs.size() < 0.99 ) cr1 = 0;
  if ( fc2 <= 2 ) cr2 = 0;

 //if (!cls.empty() and !crs.empty()) 
 //  cout << "##1 " << region_begin << " " << cls.size() << " " << crs.size() << " " << cl1 << " " << cr1 << endl;  
}

string get_pos_pair(int cl1, int cl2, int cr1, int cr2){
  if ( !cl1 and !cr1)  return "";
  std::set <int> cms;
  cms.insert(0);
  cms.insert(cl1);
  cms.insert(cl2);
  cms.insert(cr1);
  cms.insert(cr2);
  cms.erase(0);
  std::set <int>::iterator ci = cms.begin();
  stringstream ss;
  int pre_pos = *ci;
  ci++;
  for (; ci != cms.end(); ci++) {
    if ( abs(*ci - pre_pos) < MAX_POS_DIF ) {
      pre_pos = (*ci + pre_pos)/2;
    } else {
      ss << " " << pre_pos;
      pre_pos = *ci;
    }
  }
  ss << " " << pre_pos;
  return ss.str(); 
}

void regions_pos_vote( string path1, std::set < pair<int, int> > & regionPos_set, string fn_output) {
  ofstream fout(fn_output.c_str());
  fout << "region_begin region_end clipMid\n";
  for (std::set < pair <int, int> >::iterator ri = regionPos_set.begin(); ri != regionPos_set.end(); ri++ ) {
    int cl1, cl2, cr1, cr2;
    region_pos_vote(path1 + pos_to_fn(*ri), cl1, cl2, cr1, cr2) ; // at least two reads to vote for this position
//    if ( cl1 + cl2 + cr1 + cr2 >0 ) 
//      cout << pos_to_fn(*ri) << " " << cl1 << " " << cl2 << " " << cr1 << " " << cr2 << endl;
    string cs = get_pos_pair(cl1, cl2, cr1, cr2);   
    if ( cs.size() > 3 ) 
      fout << (*ri).first << " " << (*ri).second << cs << endl;
  }  
}

void neighbor_pos_filename(std::set < pair<int,int> > & regionPos_set, std::set <pair <int, int> > & fns, int clipLeft) {
  const int FLANK_REGION_LEN = 80;
  std::set < pair<int,int> >::iterator ri;
  ri = regionPos_set.find( * (fns.begin()) );
  while ( ri != regionPos_set.begin() and (*ri).first + FLANK_REGION_LEN >= clipLeft) 
    fns.insert( *ri-- );
  ri = regionPos_set.find( * (fns.rbegin()) );
  while ( ri !=regionPos_set.end() and (*ri).second - FLANK_REGION_LEN <= clipLeft ) 
    fns.insert( *ri++ );
}

int approximate_pos_pns(string path1, string fn_input, string fn_output, std::set < pair<int,int> > & regionPos_set, int pos_dif){
  ifstream fin;
  fin.open(fn_input.c_str());
  assert(fin);
  map< int, std::set <pair<int, int> > > clipLeft_fn0;
  map< int, std::set <pair<int, int> > > clipLeft_fn;
  stringstream ss;
  string line;
  int region_begin, region_end, pos2, pos1;
  getline(fin,line);
  while (getline(fin,line)) {
    ss.clear(); ss.str( line );
    ss >> region_begin >> region_end;
    while ( ss >> pos2)   // in each region, several potential breakpoints
      clipLeft_fn0[pos2].insert( make_pair(region_begin, region_end)) ;
  }
  fin.close();
  if (clipLeft_fn0.empty() ) return 0;
  map< int, std::set <pair<int, int> > >::iterator ci = clipLeft_fn0.begin();  
  int pos0 = ci->first;
  clipLeft_fn[pos0].swap(ci->second);
  ci ++;
  for (; ci != clipLeft_fn0.end(); ci++) {
    pos2 = ci->first;
    if ( pos2 - pos0 <= pos_dif ) {
      pos1 = (pos2 + pos0) / 2;
      clipLeft_fn[pos1].swap(clipLeft_fn[pos0]);
      clipLeft_fn.erase(pos0);
      pos2 = pos1;
      for ( std::set <pair<int, int> >::iterator fi = (ci->second).begin(); fi != (ci->second).end(); fi++ ) 
	clipLeft_fn[pos2].insert(*fi);
    } else {
      clipLeft_fn[pos2].swap(ci->second);
    }
    pos0 = pos2;
  }
  //cout << "size0 " << clipLeft_fn0.size() << " " << clipLeft_fn.size() << endl;
  clipLeft_fn0.clear();
  ofstream fout(fn_output.c_str());
  fout << "clipMid pn(count)\n";
  for (map< int, std::set < pair<int, int> > >::iterator pi = clipLeft_fn.begin(); pi != clipLeft_fn.end(); pi++ ) {
    int clipLeft = pi->first;
    map <string, std::set<string> > pn_qName;
    std::set < pair<int, int> > fns;
    size_t fns_size = (pi->second).size() ;  // num of regions voted for this clip pos
    fns.swap(pi->second);  
    neighbor_pos_filename(regionPos_set, fns, clipLeft);  // add more regions to fns
    assert ( fns.size() >= fns_size );
    for ( std::set < pair<int, int> >::iterator fi = fns.begin(); fi != fns.end(); fi++ ) { 
      fin.open( (path1 + pos_to_fn(*fi) ).c_str());
      string line, pn, qName, tmp1, tmp2, tmp3;
      while (getline(fin, line)) {
	ss.clear(); ss.str( line );
	int _p;
	ss >> pn >> tmp1 >> tmp2 >> tmp3 >> _p >> qName;
	if ( abs( _p - clipLeft ) <= MAX_POS_DIF) 
	  pn_qName[pn].insert(qName);  // some reads appear more than once 
      }
      fin.close();
    }
    if ( pn_qName.empty() )  
      continue;      

    fout << clipLeft;  // clipRight unknown for now
    for (map <string, std::set<string> >::iterator pi = pn_qName.begin(); pi != pn_qName.end(); pi++) 
      fout << " " << pi->first << "," << (pi->second).size() ;
    fout << endl;
  }
  fout.close();
  return 1;
}

void vote_clipPos(int clip0, int & clipl, int & clipr, string path0, map <string, int> & pn_cnt, float freq_th){
  clipl = clipr = 0;
  vector <int> clipls, cliprs;
  for (map <string, int>::iterator pi = pn_cnt.begin(); pi != pn_cnt.end(); pi++ ) {
    ifstream fin( (path0 + pi->first).c_str() );
    assert(fin);
    string line, tmpv0; 
    char left_right;
    int adj_clipPos, tmpv1, tmpv2;
    stringstream ss;
    getline(fin, line);   
    while ( getline(fin, line) ) {
      ss.clear(); ss.str( line );
      ss >> tmpv0 >> tmpv1 >> tmpv2 >> left_right >> adj_clipPos;
      //cout << clip0 << " " << tmpv2 << " " << adj_clipPos << endl;
      if ( clip0 + 300 < tmpv2  ) break;
      if ( abs(clip0 - adj_clipPos) <= MAX_POS_DIF ) {
	pi->second += 1;
	if ( left_right == 'L' )  clipls.push_back(adj_clipPos); 
	else if ( left_right == 'R' ) cliprs.push_back(adj_clipPos); 
      } 
    }
    fin.close(); 
  }
  //cout << "vote " << clip0 << " " << clipls.size() << " " << cliprs.size() << endl;
  ////if ( major_key_freq(clipls, clipl, MAX_POS_DIF, freq_th, 3) >= 2  // eg: freq_th = 0.7, at least 3 clip read exact the same   
  major_key_freq(clipls, clipl, MAX_POS_DIF, freq_th, 3);  // eg: freq_th = 0.7, at least 3 clip read exact the same 
  major_key_freq(cliprs, clipr, MAX_POS_DIF, freq_th, 3); 
}

void exact_pos_pns(string path0, string fn_input, string fn_output, float freq_th ) {
  ifstream fin(fn_input.c_str());
  assert(fin);
  ofstream fout(fn_output.c_str());
  stringstream ss;
  string line, pn, cnt, tmpv;
  int clip0, clipl, clipr;
  getline(fin,line);
  fout << "clipLeft clipRight pn(count)\n";
  while (getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> clip0;
    map <string, int> pn_cnt;
    while ( ss >> tmpv ) {
      split_by_sep (tmpv, pn, cnt, ',');
      pn_cnt[pn] = 0;
    }

    int tmpcnt;
    seqan::lexicalCast2(tmpcnt, cnt);    
    if (pn_cnt.size() <= 1 and tmpcnt <= 2 ) continue;  // at least 2 reads to support 
    
    vote_clipPos(clip0, clipl, clipr, path0, pn_cnt, freq_th);
    //cout << clip0 << " " << pn_cnt.size() << " " << clipl << " " << clipr << endl;
    if (!clipl and !clipr ) continue;
    fout << clipl << " " << clipr;
    for (map <string, int>::iterator pi = pn_cnt.begin(); pi != pn_cnt.end(); pi++ ) 
      fout << " " << pi->first;
    fout << endl;
  }
  fin.close();
  fout.close();
}

bool is_clipread(int & clipLen, seqan::CharString & clipSeq, seqan::BamAlignmentRecord & record, int thisEnd, int clipLeft, int clipRight, PAR_DELETE0 & pars) {
  if ( has_soft_first(record, CLIP_BP) and abs( record.beginPos - clipRight) < MAX_POS_DIF
       and trim_clip_soft_first(record, clipLen, clipSeq, pars.CLIP_Q) )
    return true;
  if ( has_soft_last(record, CLIP_BP) and abs( thisEnd - clipLeft ) < MAX_POS_DIF
       and trim_clip_soft_last(record, clipLen, clipSeq, pars.CLIP_Q) )
    return true;
  return false;
}

string classify_read( seqan::BamAlignmentRecord & record, int clipLeft, int clipRight, PAR_DELETE0 & pars, int & clipLen, seqan::CharString & clipSeq) {
  string iread = "useless";    
  if (!aluclip_RC_match(record)) 
    return iread;
  
  int thisEnd = record.beginPos + (int)getAlignmentLengthInRef(record);   
  if (record.rID != record.rNextId or abs(record.tLen) >= DISCORDANT_LEN) {
    if ( ( thisEnd < clipLeft + 5 and !hasFlagRC(record) ) or ( record.beginPos > clipRight - 5 and hasFlagRC(record) ) )
      iread = "alu_read";  
    if ( is_clipread(clipLen, clipSeq, record, thisEnd, clipLeft, clipRight, pars) ) {
      if (iread == "alu_read")
	iread = "alu_clip";  // is clip read + mate is Alu
      else 
	iread = "clip_read";  
    }
    return iread;    
  }
  bool read_is_left = left_read(record);
  // otherwise only look at proper mapped reads 
  if ( hasFlagRC(record) == hasFlagNextRC(record) or read_is_left == hasFlagRC(record) ) 
    return "useless";  
  if ( is_clipread(clipLen, clipSeq, record, thisEnd, clipLeft, clipRight, pars) ) 
    return "clip_read";
  int trimb, trime;
  get_trim_length(record, trimb, trime, pars.CLIP_Q);
//  cout << "###trim " << trimb << " " << trime << " " << record.beginPos << " "  << thisEnd << endl;
//  debug_print_read(record);
  if (trimb + trime < (int) length(record.seq) - CLIP_BP and count_non_match(record) <= 5 and
      record.beginPos + trimb <= min(clipLeft, clipRight) - CLIP_BP and thisEnd - trime >= max(clipLeft, clipRight) + CLIP_BP ) 
    return "skip_read";  // after trimming, fewer reads are considered as skip_read, thus less skip_clip conflict, more clip read
  int pair_begin = read_is_left ? record.beginPos : record.pNext;
  int pair_end = pair_begin + abs(record.tLen);
  if ( pair_begin < clipLeft - CLIP_BP and pair_end > clipLeft + CLIP_BP) {
    if ( (read_is_left and has_soft_first(record, CLIP_BP) ) or (!read_is_left and has_soft_last(record, CLIP_BP) ) )
      return "useless";
    return "unknow_read";
  }
  return "useless";
}

int tmp0_reads(string chrn, int clipLeft, int clipRight, BamFileHandler *bam_fh, ostream & fout, ostream & fout2,  map < int, vector<ALUREAD_INFO> > & rid_alureads, map<string, int> & rg_to_idx, PAR_DELETE0 & pars, bool test_print = false, int exact_left = 0, int exact_right = 0) {
  int region_begin = clipLeft - ALU_FLANK; // no need for larger flank region, if i only want to build consensus
  int region_end = clipRight + ALU_FLANK;	
  if (! bam_fh->jump_to_region(chrn, region_begin, region_end)) 
    return 0;
  map < seqan::CharString, string > rg_str;
  map < seqan::CharString, string > qname_iread;
  map < seqan::CharString, string > qname_clipInfo;
  map < seqan::CharString, string > qname_skipInfo;
  while ( true ) {
    seqan::BamAlignmentRecord record;
    string read_status = bam_fh->fetch_a_read(chrn, region_begin, region_end, record);
    if (read_status == "stop" ) break;
    if (read_status == "skip" or !QC_insert_read(record)) continue;      

    string iread2 = "";   
    get_mapVal(qname_iread, record.qName, iread2);  	
    if ( !iread2.empty() and iread2.substr(0,3) == "alu")  continue;
      
    int clipLen;
    seqan::CharString clipSeq;
    string iread = classify_read(record, clipLeft, clipRight, pars, clipLen, clipSeq);                  
    
///    string tmpstr1 = "FCC20N3ACXX:2:2301:9955:98649#CTTGATAG";
///    string tmpstr2 = toCString(record.qName);
///    if (tmpstr1 == tmpstr2)
///      cout << "tmpstr " << iread << endl;

    if ( iread == "alu_read" or iread == "alu_clip") {// mate is mapped to Alu
      int this_pos;
      if (hasFlagRC(record)) this_pos = -record.beginPos;  // alu is left read
      else this_pos = record.beginPos;	   // alu is right read
      ALUREAD_INFO one_aluRead = ALUREAD_INFO(record.qName, this_pos, clipLeft, record.pNext, length(record.seq), hasFlagRC(record)==hasFlagNextRC(record) );
	rid_alureads[record.rNextId].push_back(one_aluRead);
    }     
    
    bool conflict = false;
    if  ( (iread == "clip_read" and iread2 == "skip_read") or (iread == "skip_read" and iread2 == "clip_read") )
      conflict = true;	
    
    bool write_unknow = ( iread == "unknow_read") ? true : false;
    if ( !iread2.empty() and iread2 != "unknow_read" ) 
      write_unknow = false;
    
    if ( conflict or write_unknow) {
      qname_iread[record.qName] = "unknow_read";
      int unknow_exact = 0;  // unknow for now, may belong to clip reads afterwards
      if (has_soft_first(record, CLIP_BP)) unknow_exact = record.beginPos;
      else if (has_soft_last(record, CLIP_BP)) unknow_exact = - ( record.beginPos + (int)getAlignmentLengthInRef(record)); //left of read
      if ( rg_str.find(record.qName) == rg_str.end() or unknow_exact != 0) { // if exists, no replace 
	stringstream rg_ss;
	int rgIdx = get_rgIdx(rg_to_idx, record);
	rg_ss << rgIdx << ":" << abs(record.tLen) << ":" << min(record.beginPos, record.pNext) << ":" << unknow_exact;
	rg_str[record.qName] = rg_ss.str();
      }
      continue;
    } 
    if (iread == "clip_read" or iread == "alu_clip") {
      qname_iread[record.qName] = iread;
      stringstream ssf;
      ssf << clipSeq  << " " << clipLen << " ";
      if (record.cigar[0].operation == 'S') ssf << record.beginPos; 	    
      else ssf << record.beginPos + (int)getAlignmentLengthInRef(record);	    
      int rgIdx = get_rgIdx(rg_to_idx, record);
      ssf << " " << get_cigar(record) << " " << rgIdx << " " << record.beginPos << " " ;
      if ( record.rID == record.rNextId )  ssf << record.beginPos + record.tLen;
      else  ssf << "-1";
      qname_clipInfo[record.qName] = ssf.str();
      continue;
    }		
    if (iread == "skip_read") {      
      qname_iread[record.qName] = iread;	
      stringstream ssf;
      ssf << record.beginPos << " " << record.beginPos + (int)getAlignmentLengthInRef(record); 
      qname_skipInfo[record.qName] = ssf.str();
      continue;
    }
  }   
  
  if (!test_print) {
    vector<string>  str_unknowreads;   
    bool flag = false;
     for (map < seqan::CharString, string >::iterator qi = qname_iread.begin(); qi != qname_iread.end(); qi++ ) {
       if ( qi->second == "skip_read") {
         fout2 << clipLeft << " " << qname_skipInfo[ qi->first] << " " << qi->first << endl;  
       } else if ( qi->second == "clip_read" or  qi->second == "alu_clip") {
         fout << clipLeft << " " << qname_clipInfo[ qi->first] << " " << qi->first << endl;
         flag = true;   ///// at least one clip_read, otherwise ignore, if (!rid_alureads.empty()) flag = true;
       } else if ( qi->second == "unknow_read") {
         str_unknowreads.push_back( rg_str[qi->first] );	
       }
     }   
     
     if (flag and !str_unknowreads.empty()) {  
       fout2 << clipLeft << " UR " << str_unknowreads.size() ;
       for  ( vector< string > ::iterator si = str_unknowreads.begin(); si != str_unknowreads.end(); si++ )
         fout2 << " " << *si ;
       fout2 << endl;
     }
  } else {  // print unknown reads
    std::set <seqan::CharString> uds;
    for (map < seqan::CharString, string >::iterator qi = qname_iread.begin(); qi != qname_iread.end(); qi++ ) 
      if ( qi->second == "unknow_read") 
	uds.insert(qi->first);
    cout << "unknown reads " << uds.size() << endl;
    bam_fh->jump_to_region(chrn, region_begin, region_end);
    int tmpmax = max (exact_left, exact_right);
    int tmpmin = (!exact_left or !exact_right) ? tmpmax : min(exact_left, exact_right);
    while ( true ) {
      seqan::BamAlignmentRecord record;
      string read_status = bam_fh->fetch_a_read(chrn, region_begin, region_end, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_insert_read(record)) continue;      
      if ( uds.find(record.qName) != uds.end() ) {
	int pair_begin = min(record.beginPos, record.pNext);
	if ( pair_begin < tmpmin - CLIP_BP_MID and pair_begin + abs(record.tLen)> tmpmax + CLIP_BP_MID) 
	  debug_print_read(record);
      }
    }
  }
  return 1;
}

// too slow, not used any more
void alumate_reads1(BamFileHandler *bam_fh, string file_db, ofstream & fout, int query_rid, vector<ALUREAD_INFO> & aluinfos) {
  map <string, string> qnames;
  if (!file_db.empty()) {
    ifstream fin(file_db.c_str());
    assert(fin);
    stringstream ss;
    string line, qname, tmpv, alu_type;
    getline(fin,line);
    while (getline(fin, line)) {
      ss.clear(); ss.str(line); 
      ss >> qname ;
      for ( int i = 0; i < 8; i++) ss >> tmpv;
      ss >> alu_type;
      qnames[qname] = alu_type;
    }    
    fin.close();
  }  
  for ( vector<ALUREAD_INFO>::iterator ai = aluinfos.begin(); ai != aluinfos.end(); ai++ ) {
    if ( !bam_fh->jump_to_region( query_rid, (*ai).pos - 20, (*ai).pos + 20) ) continue;
    string mapped_db = "na";
    string _qn =  toCString( (*ai).qName);
    get_mapVal(qnames, _qn, mapped_db);
    if ( mapped_db == "na") continue;
    seqan::BamAlignmentRecord record;
    while (bam_fh->fetch_a_read(record) and record.beginPos <= (*ai).pos + 3) {
      if (record.qName == (*ai).qName) {
	string record_seq = toCString(record.seq);
	if ((*ai).sameRC)  seqan::reverseComplement(record_seq); 	
	fout << (*ai).clipLeft << " " << record_seq << " " << (*ai).sameRC << " " << mapped_db << endl;
	break;
      }
    }
  }
}

void alumate_reads2(string file_db, ofstream & fout, int rid,  vector<ALUREAD_INFO> & aluinfos, string alu_type_default) {
  map <string, string> qnames;
  if ( !file_db.empty()) {
    ifstream fin(file_db.c_str());
    assert(fin);
    stringstream ss;
    string line, qname, tmpv, alu_type;
    getline(fin,line);
    while (getline(fin, line)) {
      ss.clear(); ss.str(line); 
      ss >> qname ;
      for ( int i = 0; i < 8; i++) ss >> tmpv;
      ss >> alu_type;
      qnames[qname] = alu_type;
    }
  fin.close();
  }
  for ( vector<ALUREAD_INFO>::iterator ai = aluinfos.begin(); ai != aluinfos.end(); ai++ )  {
    string atype = alu_type_default;
    string qx = toCString( (*ai).qName);
    get_mapVal(qnames, qx, atype);
    fout << (*ai).clipLeft << " " << atype << " " << rid << " " << (*ai).pos << " " << (*ai).sameRC << " "
	 << (*ai).readLen << " " <<  (*ai).qName << " " << (*ai).this_pos << endl;
  }
}

void count_clipr(string fn_clip, map < int, pair <int, int> > & exact_pos, map <int, set<string> > & aluclip_cnts, map <int, vector< string > > & unknow_str) {
  // if not in exact_pos, write into unknow_str
  ifstream fin ( fn_clip.c_str()) ;
  assert(fin);
  stringstream ss;
  string line, tmpv1, tmpv2, rgIdx, qName;
  int pos, clipc, clipLen, p1, p2;
  getline(fin, line); // read header
  while (getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> pos ; 
    map < int, pair <int, int> >::iterator ei = exact_pos.find(pos); 
    if ( ei == exact_pos.end() )   // check positions that failed !!!
      continue;
    int exact_left = (ei->second).first;
    int exact_right = (ei->second).second;    
    ss >> tmpv1 >> clipLen >> clipc >> tmpv2 >> rgIdx >> p1 >> p2 >> qName;   // p1, p2 is pair_begin and pair_end
    bool check_unknow = false;
    if (clipLen < 0) {
      if ( ( exact_left and abs(clipc - exact_left) <= CLIP_BP_LEFT ) or 
	   ( !exact_left and abs(clipc - exact_right) <= MAX_POS_DIF ) )
	aluclip_cnts[pos].insert(qName);
      else check_unknow = false;
    } else if (clipLen > 0 ) {
      if ( ( exact_right and abs(clipc - exact_right) <= CLIP_BP_LEFT ) or 
	   ( !exact_right and abs(clipc - exact_left) <= MAX_POS_DIF ) )
	aluclip_cnts[pos].insert(qName);
      else check_unknow = false;
    } 
    if ( check_unknow and p2 > 0 and abs(p2 - p1) < DISCORDANT_LEN ) {
      check_unknow = false;
      if ( exact_left and exact_right) {
	if ( min(p1, p2) <= min(exact_left, exact_right) - CLIP_BP and max(p1, p2) >= max(exact_left, exact_right) + CLIP_BP )  
	  check_unknow = true;
      } else {
	if ( min(p1, p2) <= max(exact_left, exact_right) - MAX_POS_DIF and max(p1, p2) >= max(exact_left, exact_right) + MAX_POS_DIF)  
	  check_unknow = true;	  
      }
      if (check_unknow) {
	string tmps = rgIdx + ":" + int_to_string( abs(p2 - p1) );
	unknow_str[pos].push_back( tmps );
      }
    }
  } 
}

int count_alur(string file_input, map < int, pair <int, int> > & exact_pos, map <int, set <string> > & aluclip_cnts, map < int, string > & rid_chrn, string file_fa, AluconsHandler *alucons_fh, vector <string> & chrns, bool alumate_align) {
  ifstream fin(file_input.c_str());
  assert(fin);
  stringstream ss;
  string line;
  getline(fin, line); // read header
  vector < REALIGNS > to_realign;  // rid_pos_alignBegin_alignEnd
  while (getline(fin, line)) {
    ss.clear(); ss.str(line); 
    string alu_type, qName;
    int pos, rid, aluBegin, sameRC, read_len, thisBegin;    
    ss >> pos;
    map < int, pair <int, int> >::iterator ei = exact_pos.find(pos); 
    if ( ei == exact_pos.end() )   // check positions that failed !!!
      continue;
    ss >> alu_type >> rid >> aluBegin >> sameRC >> read_len >> qName >> thisBegin;
    map <int, set <string> >::iterator aci = aluclip_cnts.find(pos);
    if ( aci != aluclip_cnts.end() and (aci->second).find(qName) != (aci->second).end() )
      continue;
    int exact_left = (ei->second).first;
    int exact_right = (ei->second).second;    
    int tmpmax = max (exact_left, exact_right);
    int tmpmin = (!exact_left or !exact_right) ? tmpmax : min(exact_left, exact_right);
    if ( (thisBegin < 0 and abs(thisBegin) < tmpmax - CLIP_BP) or  // RC should match
	 (thisBegin > 0 and abs(thisBegin) > tmpmin + CLIP_BP) )
      continue;    
    if ( alu_type.substr(0,3) == "Alu" ) {  // match one Alu in Database
      aluclip_cnts[pos].insert(qName); 
    } else if ( aci !=aluclip_cnts.end() ) {  // don't bother to realign if clip read does not exists
      REALIGNS ra = {rid, pos, aluBegin, aluBegin + read_len, qName};
      to_realign.push_back(ra);
    } 
  }  
  fin.close();
  
  if (!alumate_align) 
    return 0;

  FastaFileHandler * fasta_fh = NULL;
  string chrn_pre = "chr0";
  map <int, int > pre_count;
  map <int, int > toalign_count;
  map <int, int > added_count;

  for (vector < REALIGNS >::iterator ti = to_realign.begin(); ti != to_realign.end(); ti ++ ) {
    string chrn = rid_chrn[(*ti).rid];
    if ( std::find( chrns.begin(), chrns.end(), chrn) == chrns.end())
      continue;
    if ( chrn_pre != chrn ) {
      if ( chrn_pre != "chr0") delete fasta_fh;
      chrn_pre = chrn;
      fasta_fh = new FastaFileHandler(replace_str0_str(file_fa, chrn, "chr0"), chrn);
    }
    seqan::CharString ref_fa;    
    float sim_rate;
    fasta_fh->fetch_fasta_upper((*ti).aBegin, (*ti).aEnd, ref_fa);
    string alu_type = align_alu_cons_call(ref_fa, alucons_fh, sim_rate, 0.8);
    // keep track of counts
    if ( pre_count.find((*ti).pos_key) == pre_count.end() ) {
      pre_count[(*ti).pos_key] = aluclip_cnts[(*ti).pos_key].size();
    }
    
    addKey(toalign_count, (*ti).pos_key, 1);
    if (!alu_type.empty()) {
      aluclip_cnts[(*ti).pos_key].insert((*ti).qName);
      //cout << "aligned " << chrn << " " << alu_type << " " << (*ti).pos_key << " " << (*ti).qName << endl;
      addKey(added_count, (*ti).pos_key, 1);
    }
  }
  if (chrn_pre != "chr0" ) delete fasta_fh;
  // print out counts
  cout << "debug print: pos, preCnt, toalignCnt, addedCnt\n";
  for (map <int, int >::iterator pi = pre_count.begin(); pi != pre_count.end(); pi ++ ) {
    int addcnt = 0;
    get_mapVal( added_count, pi->first, addcnt);
    cout << pi->first << " " << pi->second << " " << toalign_count[pi->first] << " " << addcnt << endl;
  }
  return 0;
}


void count_skipr(string file_input, map < int, pair <int, int> > & exact_pos, map <int, set <string> > & aluclip_cnts, map <int, int > & skip_cnts, map <int, vector <string> > & unknow_str){
  ifstream fin(file_input.c_str());
  assert(fin);
  stringstream ss;
  string line, col2, ustr;
  int pos;
  getline(fin, line); // read header
  while (getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> pos >> col2;
    if ( aluclip_cnts.find(pos) == aluclip_cnts.end()) // not interesting if neither clip or alu reads exist
      continue;
    map < int, pair <int, int> >::iterator ei = exact_pos.find(pos); 
    int exact_left = (ei->second).first;
    int exact_right = (ei->second).second;    
    int tmpmax = max (exact_left, exact_right);
    int tmpmin = (!exact_left or !exact_right) ? tmpmax : min(exact_left, exact_right);
    if ( col2 == "UR" ) {
      int unknowCnt;
      int aluC = 0;
      ss >> unknowCnt;
      string token, idx, s_insert_len;
      int pair_begin, insert_len, unknow_exact;
      for (int i = 0; i < unknowCnt; i++) {
	getline(ss, idx, ':');
	getline(ss, s_insert_len, ':');
	seqan::lexicalCast2(insert_len, s_insert_len);      
	getline(ss, token, ':');
	seqan::lexicalCast2(pair_begin, token);      
	getline(ss, token, ' '); 
	seqan::lexicalCast2(unknow_exact, token);  // check if they are real unknow reads or mis_classified clip reads 
	if ( unknow_exact < 0 ) {
	  if ( (exact_left and abs(-unknow_exact - exact_left) < CLIP_BP ) or 
	       (!exact_left and abs(-unknow_exact - exact_right) <= CLIP_BP_MID ) )
	    aluclip_cnts[pos].insert("UR"+int_to_string(aluC++));
	} else if  ( unknow_exact > 0 ) {
	  if ( (exact_right and abs(unknow_exact - exact_right) < CLIP_BP ) or 
	       (!exact_right and abs(unknow_exact - exact_left) <= CLIP_BP_MID ) ) 
	    aluclip_cnts[pos].insert("UR"+int_to_string(aluC++));
	} else if ( pair_begin < tmpmin - CLIP_BP_MID and pair_begin + insert_len > tmpmax + CLIP_BP_MID) {
	  unknow_str[pos].push_back(idx + ":" + s_insert_len);
	}	  
      }
    } else {
      int thisBegin, thisEnd;
      seqan::lexicalCast2( thisBegin, col2);
      ss >> thisEnd;
      if (thisBegin <= tmpmin - CLIP_BP and thisEnd >= tmpmax + CLIP_BP)  // prefer under call      
	addKey( skip_cnts, pos, 1);
    }
  }
}

void read_seq(string pn, string fn, map <int, vector<string> > & pos_seqs, string skip_mk = ""){
  ifstream fin( fn.c_str() );
  assert (fin);
  string line,tmpv;
  int pos;
  stringstream ss;
  getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pos >> tmpv;
    if (tmpv == skip_mk) continue;
    pos_seqs[pos].push_back( replace_str0_str(line, pn, int_to_string(pos) ));
  }
  fin.close();
}

void read_clip_pass(string fn_clip_pass, map < int, pair <int, int> > & exact_pos) {
  string line;
  stringstream ss;
  int pos, clipl, clipr;
  ifstream fin ( fn_clip_pass.c_str()) ;
  assert(fin);
  getline(fin, line); // read header
  while (getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> pos >> clipl >> clipr;
    exact_pos[pos] = make_pair(clipl, clipr);
  }
  fin.close();
}

void write_tmp1(vector <string> & chrns, string file_tmp1, string fn_clip_pass0, string fn_clip0, string fn_alu0, string fn_su0, map < int, string > & rid_chrn, string file_fa, AluconsHandler *alucons_fh) {
  ofstream fout(file_tmp1.c_str());
  fout << "chr insertMid debugInfo aluclipCnt skipCnt unknowCnt unknowStr\n";  
  for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
    string chrn = *ci;
    string fn_clip_pass = replace_str0_str(fn_clip_pass0, chrn, "chr0");
    map < int, pair <int, int> > exact_pos;
    read_clip_pass(fn_clip_pass, exact_pos);
    map <int, set <string> > aluclip_cnts;
    map <int, vector <string> > unknow_str;
    map <int, int > skip_cnts;
    string fn_clip = replace_str0_str(fn_clip0, chrn, "chr0");
    string fn_alu = replace_str0_str(fn_alu0, chrn, "chr0");
    string fn_su = replace_str0_str(fn_su0, chrn, "chr0");    

    count_clipr(fn_clip, exact_pos, aluclip_cnts, unknow_str);
    count_alur(fn_alu, exact_pos, aluclip_cnts, rid_chrn, file_fa, alucons_fh, chrns, false);  
    //count_alur(fn_alu, exact_pos, aluclip_cnts, rid_chrn, file_fa, alucons_fh, chrns, true);  //re_align, might be slow 
    count_skipr(fn_su, exact_pos, aluclip_cnts, skip_cnts, unknow_str);
    for (map <int, set <string> >::iterator ai = aluclip_cnts.begin(); ai != aluclip_cnts.end(); ai++ ) {
      int exact_left = (exact_pos[ai->first]).first;
      int exact_right = (exact_pos[ai->first]).second; 
      int exact_mid = (!exact_left or !exact_right) ? max(exact_left, exact_right) : (exact_left + exact_right)/2;
      int skipCnt = 0;
      get_mapVal( skip_cnts, ai->first, skipCnt);
      fout << chrn << " " << exact_mid << " " << exact_left << "," << exact_right << " " << (ai->second).size()
	   << " " << skipCnt << " ";
      map <int, vector <string> >::iterator ui = unknow_str.find(ai->first);
      if ( ui == unknow_str.end() ) {
	fout << "0";
      } else {
	fout << (ui->second).size();
	for ( vector <string>::iterator it = (ui->second).begin(); it != (ui->second).end(); it++) 
	  fout << " " << *it;
      }
      fout << endl;
    }
  }
  fout.close();
}

void write_tmp2(string fn_tmp1, string fn_tmp2, map <int, EmpiricalPdf *> & pdf_rg, float log10RatioUb, int fixed_len){
  ifstream fin0(fn_tmp1.c_str());
  assert(fin0);
  ofstream fout(fn_tmp2.c_str());
  fout << "chr insertBegin insertEnd insertLen midCnt clipCnt unknowCnt 00 01 11\n";  
  stringstream ss;
  string line, output_line;
  getline(fin0, line); // read header
  map<string, map<int, string> > tmp1_info;
  while (getline(fin0, line)) 
    if (parseline_del_tmp0(line, output_line, pdf_rg, log10RatioUb, fixed_len, 2))  // at least 2 counts support clip 
      fout << output_line << endl;
  fin0.close();
  fout.close();
}

int main( int argc, char* argv[] )
{
  string config_file = argv[1];
  string opt =  argv[2];
  if (argc < 3) exit(1);
  boost::timer clocki;    
  clocki.restart();
  ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
  
  int alucons_len = seqan::lexicalCast<int> (cf_fh.get_conf("alucons_len"));
  float consensus_freq = seqan::lexicalCast<float> (cf_fh.get_conf("consensus_freq"));
  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  chrns.push_back("chrX");
  chrns.push_back("chrY");
  
  string path0 = cf_fh.get_conf( "file_alu_insert0") ;    
  check_folder_exists(path0);
  string path1 = cf_fh.get_conf( "file_alu_insert1") ;    
  check_folder_exists(path1);
  string pathClip = path1 + "clip/";
  check_folder_exists(pathClip);     
  string pathCons = path1 + "cons/";
  check_folder_exists(pathCons);
  string pathDel0 = path1 + "fixed_delete0/";
  check_folder_exists(pathDel0);    

  map <int, string> ID_pn;
  get_pn(cf_fh.get_conf( "file_pn"), ID_pn);
  string file_dist_prefix = cf_fh.get_conf("file_dist_prefix");
  string file_pn_used = cf_fh.get_conf( "file_pn_used");
  std::set <string> pns_used;
  
  if ( opt != "preprocess" and opt != "write_tmps_pn" ) read_file_pn_used( file_pn_used, pns_used); // some pn is ignored, due to too many reads

  if ( opt == "preprocess" ) {
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++ ) {
      check_folder_exists(path0 + *ci +  "/" );
      check_folder_exists(pathClip + *ci +  "/");
      check_folder_exists(pathClip + *ci + "_pos/");
      check_folder_exists(pathCons + *ci);
      check_folder_exists(pathCons + *ci + "_pos/");
    } 

    string bam_input = cf_fh.get_conf( "file_bam_prefix") + ID_pn[0] + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bam_input+".bai");    
    string bam_rid_chrn = cf_fh.get_conf("bam_rid_chrn");
    if ( check_file_size(bam_rid_chrn) <= 0 ) {
      ofstream fmap(bam_rid_chrn.c_str());
      for ( map<int, string>::iterator bi = bam_fh->rID_chrn.begin(); bi !=  bam_fh->rID_chrn.end(); bi++) 
	fmap << bi->first << " " << bi->second << endl;
      fmap.close();
    }
    delete bam_fh;
    
  } else if (opt == "write_tmps_pn") {     // write tmp0 files 
    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    string file1 = path0 + pn + ".tmp1";
    map<int, string>::iterator rc ;
    //// write pn.tmp1.
    map <string, int> rg_to_idx;    
    parse_reading_group( get_name_rg(file_dist_prefix, pn), rg_to_idx );    

//    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input);
//    string header1 = "qname A_chr B_chr A_beginPos B_beginPos A_isRC B_isRC len_read rgIdx\n";
//    alu_mate_flag(bam_fh, file1, header1, rg_to_idx);
//    //  cout << "done with reading this bam file, time used " << clocki.elapsed() << endl; // about 2 hours    
//
//    //// write chr*/pn.tmp1  use database to filter alu mate
//    string header2 = "qname this_rID this_pos this_isRC len_read alu_rID alu_pos alu_isRC rgIdx alu_type\n";
//    MapFO fileMap;    
//    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) {
//      fileMap[rc->first] = new ofstream( ( path0 + rc->second + "/" + pn + ".tmp1").c_str() );
//      assert(fileMap[rc->first]);
//      *(fileMap[rc->first]) << header2 ;
//    }     
//    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) 
//      keep_alu_mate(bam_fh, rc->first,  file1, fileMap, cf_fh.get_conf( "file_alupos_filter") + "alu_" + rc->second);           
//    close_fhs(fileMap, bam_fh->rID_chrn);    
//   
//    int col_idx = 0;
//    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) {
//      string file_tmp1 =  path0 + rc->second +"/" + pn + ".tmp1";
//      if (!col_idx) col_idx = get_col_idx(file_tmp1, "this_pos"); 
//      sort_file_by_col<int> (file_tmp1, col_idx, true);
//    }
//    move_files(path0+"tmp1s/", file1);    
//   
//    alumate_counts_filter(path0, pn, ".tmp1", ".tmp2", ".tmp3", chrns, LEFT_PLUS_RIGHT); 
//     delete bam_fh;
    
    //// filter out rep regions,  combine rep regions(dist < REP_MASK_JOIN bp)
    const int REP_MASK_JOIN = 200;
    RepMaskPos *repmaskPos;
    std::set < string > skip_type;
    repmaskPos = new RepMaskPos(cf_fh.get_conf( "file_repeatMask"), chrns, REP_MASK_JOIN); 
    filter_location_rep( path0, pn, ".tmp3", ".tmp3st", chrns, repmaskPos);    
    delete repmaskPos;
    
  }  else if (opt == "combine_pos_pns") { // combine positions from multiple individuals

    const int NUM_PN_JOIN = 10;  // max number of pns to sort and join, only matters how fast it runs
    ifstream fin(file_pn_used.c_str());
    if (!fin) {
      cerr << "which individuals should be used for further analysis?\n";
      cerr << "please write this information should be in file " << file_pn_used << endl;
      return 1;
    }
    map<int, string> id_pn_map;   
    int i = 0;
    string _pn;
    while (fin >> _pn) id_pn_map[i++] = _pn;
    fin.close();              
    cout << id_pn_map.size() << " individuals out of " << ID_pn.size() << " individuals are used, change "
	 << file_pn_used << " if you want to modify\n";

    string output_prefix = path1 + "insert_pos.";
    /// output "sum_LR_alumate regionBegin regionEnd pns"    
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      ////combine_pn_pos(*ci, id_pn_map, ".tmp3", output_prefix, NUM_PN_JOIN, path0);  // used for simulation data 
      combine_pn_pos(*ci, id_pn_map, ".tmp3st", output_prefix, NUM_PN_JOIN, path0);  // takes about 12 mins to run !   
      system( ("rm " + output_prefix + *ci + "*tmp?").c_str()  );      
    }        
    
  } else if (opt == "clipReads_pn") { // clip reads for each pn
    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];    
    if ( pns_used.find(pn) == pns_used.end() ){
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bam_input+".bai");
    string file_fa = cf_fh.get_conf("file_fa_prefix");    
    int CLIP_Q = seqan::lexicalCast<int> (cf_fh.get_conf("CLIP_Q"));
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {      
      string fin_pos = path1 + "insert_pos." + *ci; // eg: 10181 positions for chr1      
      string fout_reads = pathClip + *ci + "/" + pn;
      FastaFileHandler * fasta_fh = new FastaFileHandler(file_fa + *ci + ".fa", *ci);
      // move clip pos such that matched length is maximized 
      clipreads_at_insertPos(pn, *ci, bam_fh, fasta_fh, fin_pos, fout_reads, CLIP_Q);  
      filter_noise(fout_reads, MAX_POS_DIF, 2);
      delete fasta_fh;
    }
    delete bam_fh;
    cout << "output to " << pathClip << "chr*/" << pn  << endl;

  } else if ( opt == "clipReads_pns" ) { 

    assert (argc == 4);
    string chrn = argv[3];    
    string tmp_path0 = pathClip + chrn + "/";
    string tmp_path1 = pathClip + chrn + "_pos/";
    string file_clipRegion = pathClip + chrn + ".clip_region";
    string file_clipPos = pathClip + chrn + ".clip_pn";
    string tmp_path2 = pathCons + chrn + "/";  // to be used in next step 

    system( ("rm "+ tmp_path1 + "*").c_str() );      
    system( ("rm "+ tmp_path2 + "*").c_str() );      
    combine_clipreads_by_pos(pns_used, tmp_path0, tmp_path1);
    std::set < pair<int, int> > regionPos_set;    
    assert( nonempty_files_sorted( tmp_path1, regionPos_set) );
    if (regionPos_set.empty()) return 0;
    regions_pos_vote(tmp_path1, regionPos_set, file_clipRegion);   
    if ( approximate_pos_pns(tmp_path1, file_clipRegion, file_clipPos+".tmp", regionPos_set, MAX_POS_DIF))
      exact_pos_pns(tmp_path0, file_clipPos+".tmp", file_clipPos, consensus_freq);
    
  } else if ( opt == "write_tmp0_pn" ) {   // write alu and clip reads for each PN

    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];
    if ( pns_used.find(pn) == pns_used.end() ){ 
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bam_input+".bai");    
    map<string, int> rg_to_idx;
    parse_reading_group( get_name_rg(file_dist_prefix, pn), rg_to_idx );    
    PAR_DELETE0 pars = { seqan::lexicalCast<int> (cf_fh.get_conf("CLIP_Q")) };
    
    for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){      
      vector< pair<int, int> > insert_pos;      
      string file_clipPos = pathClip + *ci + ".clip_pn";      
      // write down clip reads 
      string file_clip2 = pathCons + *ci + "/" + pn + ".clip";
      ofstream fclip ( file_clip2.c_str());  
      fclip << "clipPos seq clipLen clipPos_by_cigar cigar rgIdx beginPos pairEnd qName\n";   
      // use database, faster, also discard half of the reads
      string file_alu2 = pathCons + *ci + "/" + pn + ".aludb";
      ofstream falu ( file_alu2.c_str());  
      falu << "clipPos alu_type rid aluBegin sameRC read_len qName thisBegin\n";      
      // mapping span insertion breakpoint
      string file_skip = pathCons + *ci + "/" + pn + ".skip";  
      ofstream fskip ( file_skip.c_str());  
      fskip << "clipPos thisBegin(UR) thisEnd(Cnt) qName(Str)\n";

      if ( !read_first2col( file_clipPos , insert_pos, true) ) {
	fclip.close(); falu.close(); fskip.close(); continue;	
      }
      map < int, vector<ALUREAD_INFO>  > rid_alureads;       // info of alu reads 
      for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) 
	tmp0_reads(*ci, (*pi).first, (*pi).second, bam_fh, fclip, fskip, rid_alureads, rg_to_idx, pars);            
      
      fclip.close();      
      fskip.close();
      
      if ( !rid_alureads.empty() ) {
	for ( map < int, vector<ALUREAD_INFO> >::iterator ri = rid_alureads.begin(); ri != rid_alureads.end(); ri++) {
	  string chrn, file_db;
	  if ( bam_fh->get_chrn(ri->first, chrn) ) file_db =  path0 + chrn +"/" + pn + ".tmp1";	    
	  else file_db = "";
	  alumate_reads2(file_db, falu, ri->first,  ri->second, "na");
	}
	falu.close();
	sort_file_by_col<int> (file_alu2, 1, true);  // sort by clip pos      
      } else 
	falu.close();
    }   
    delete bam_fh;

  } else if (opt == "clipPos_pns" ) {  // make sure each insertion has exact breakpoint 

    assert (argc == 4);
    string chrn = argv[3];
    string pathCons1 = pathCons + chrn + "/" ;
    map <int, vector<string> > pos_seqs;
    for (std::set <string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) 
      read_seq(*pi, pathCons1 + *pi + ".clip", pos_seqs);

    string file_clip_pass = pathCons + chrn + "_clip_pass" ;
    ofstream fout2( file_clip_pass.c_str() );
    fout2 << "file_name clipLeft clipRight\n";
    if (pos_seqs.empty()) {
      fout2.close(); 
      return 0;
    }
    //// rewrite reads to another folder      
    string pathCons2 = pathCons + chrn + "_pos/" ;
    system( ("rm "+ pathCons2 + "*").c_str() );      
    for (map <int,  vector<string> >::iterator pi = pos_seqs.begin(); pi != pos_seqs.end() ; pi++ ) {
      string file_cons1 = pathCons2 + int_to_string(pi->first) + ".clip" ;  // rewrite clip reads 
      ofstream fout1(file_cons1.c_str());
      for ( vector< string > ::iterator si = (pi->second).begin(); si != (pi->second).end(); si++ )
	fout1 << *si << endl;
      sort_file_by_col<string> (file_cons1, 1, true); // sort by pn
      fout1.close();		
    } 
    
    for (map <int,  vector<string> >::iterator pi = pos_seqs.begin(); pi != pos_seqs.end() ; pi++ ) {
      int posl = pi->first;
      vector <int> clipls, cliprs;
      for ( vector<string>::iterator si = (pi->second).begin(); si != (pi->second).end(); si++ ) 
	parse_line1(*si, posl, CLIP_BP_MID, clipls, cliprs);
      int clipl, clipr;
      int n_clipl = major_key_freq(clipls, clipl, CLIP_BP_LEFT, consensus_freq, 2);   
      int n_clipr = major_key_freq(cliprs, clipr, CLIP_BP_LEFT, consensus_freq, 2);
      if (!clipl and !clipr ) continue;
      int pos_dif =  (!clipl or !clipr) ? 0 : abs(clipl - clipr);
      ///cout << "clipl " << clipl << " " << n_clipl <<  " clipr " << clipr << " " << n_clipr << endl;
      if ( pos_dif > MAX_POS_DIF ){ 
	if ( n_clipl <= 1 and n_clipr <= 1 ) 
	  continue;
	int nl = clipls.size();
	int nr = cliprs.size ();
	if ( nl > nr and nr <= 5)
	  fout2 << pi->first << " " << clipl << " 0 " <<  endl;	  	
	else if ( nr > nl and nl <= 5)
	  fout2 << pi->first << " 0 " << clipr << " " <<  endl;	  	
	else
	  cout << "err? " << pi->first << " " << clipl << " " << clipr << " " <<  endl;	  
      } else {	
	fout2 << pi->first << " " << clipl << " " << clipr << " " <<  endl;	  
      }      
    }
    fout2.close() ;
    
  } else if (opt == "write_tmp2_pn") { 
    
     assert (argc == 4);
     string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];
     if ( pns_used.find(pn) == pns_used.end() ){ 
       cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
       return 0;
     }
     
     string file_clip_pass = pathCons + "chr0_clip_pass" ;
     string file_clip = pathCons + "chr0/" + pn + ".clip"; 
     string file_alu = pathCons + "chr0/" + pn + ".aludb"; 
     string file_su = pathCons + "chr0/" + pn + ".skip"; // include skip and unknown reads
     string file_tmp1 = pathDel0 + pn + ".tmp1";
     map < int, string > rid_chrn;
     get_chrn(cf_fh.get_conf("bam_rid_chrn"), rid_chrn);
     string file_fa = cf_fh.get_conf("file_fa_prefix") + "chr0.fa";
     get_chrn(cf_fh.get_conf("bam_rid_chrn"), rid_chrn);
     AluconsHandler *alucons_fh = new AluconsHandler(cf_fh.get_conf("file_alu_cons"));
     write_tmp1(chrns, file_tmp1, file_clip_pass, file_clip, file_alu, file_su, rid_chrn, file_fa, alucons_fh);	    
     delete alucons_fh;
     string file_tmp2 = pathDel0 + pn + ".tmp2";
     map <int, EmpiricalPdf *> pdf_rg;    
     string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
     read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
     float log10RatioUb = seqan::lexicalCast<float> (cf_fh.get_conf("LOG10_RATIO_UB"));
     write_tmp2(file_tmp1, file_tmp2, pdf_rg, log10RatioUb, alucons_len);     
     EmpiricalPdf::delete_map(pdf_rg);
     move_files(pathDel0 + "tmp1s/", file_tmp1);
     move_files(pathDel0 + "tmp2s/", file_tmp2);

   } else if (opt == "fixed_vcf_pns") {
    vector <string> pns;
    read_file_pn_used(file_pn_used, pns); 
    string path_input = pathDel0 + "tmp2s/";
    string tmp_file_pn = path_input + *(pns.begin()) + ".tmp2";
    int col_idx =  get_col_idx(tmp_file_pn, "00");
    
    cout << "col_idx " << col_idx << endl;
    string vcf_prefix = pathDel0 + int_to_string( pns.size());

    combine_pns_llh(path_input, ".tmp2", vcf_prefix + ".pos.tmp", pns, chrns, col_idx);
    combine_pns_vcf(path_input, ".tmp2", vcf_prefix + ".vcf.tmp", pns, chrns, col_idx);  
    
    map <string, std::set<int> > chrn_aluBegin;
    string fn_rm = vcf_prefix + ".pos.torm"; 
    write_rm2(vcf_prefix + ".pos.tmp", fn_rm, chrn_aluBegin, 1.92, false); //qchisq(0.95, df=1)  [1] 3.841459 
    filtered_vcf(vcf_prefix + ".vcf.tmp", vcf_prefix + ".vcf", -1, chrn_aluBegin);
  
  } else if (opt == "debug0") { 

    string chrn = "chr4";
    int refBegin = 77948411 - 200 ; 
    int refEnd = 77948411;

    string file_fa = cf_fh.get_conf("file_fa_prefix");    
    FastaFileHandler * fasta_fh = new FastaFileHandler(file_fa + chrn + ".fa", chrn);
    seqan::CharString ref_fa;
    fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
    cout << chrn << " " << refBegin << " " << refEnd << endl;
    cout << ref_fa << endl;

  } else if (opt == "debug1") {  // check a region 
    
    string pn = "1006-02";
    string chrn = "chr1";
    int clipLeft = 162628642;
    int clipRight = 162628619;
    
    map<string, int> rg_to_idx;
    parse_reading_group( get_name_rg(file_dist_prefix, pn), rg_to_idx );    
    PAR_DELETE0 pars = { seqan::lexicalCast<int> (cf_fh.get_conf("CLIP_Q")) };
    map < int, vector<ALUREAD_INFO>  > rid_alureads;       // info of alu reads 

    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bai_input);
    tmp0_reads(chrn, clipLeft, clipRight, bam_fh, cout, cout, rid_alureads, rg_to_idx, pars, true, clipLeft, clipRight);    
    delete bam_fh;

  } else if (opt == "debug2") { 

    string pn = "1006-01";
    map <int, EmpiricalPdf *> pdf_rg;    
    string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
    read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
    string line, output_line;
    //line = "chr1 162628630 162628642,162628619 2 6 3 0:483 0:517 0:473";
    line = "chr1 162628630 162628642,162628619 3 0 1 1:163";
    if (parseline_del_tmp0(line, output_line, pdf_rg, 3, alucons_len, 2, true))
      cout << output_line << endl;

    line = "chr1 162628630 162628642,162628619 3 0 0";
    if (parseline_del_tmp0(line, output_line, pdf_rg, 3, alucons_len, 2, true))
      cout << output_line << endl;


  } else {

    cout << "unknown option ! \n";

  }

  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
