#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"

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

void alumate_counts_filter(string fn_input, string fn_output, string fn_output2, vector<string>  &chrns, int th_cnt){
  int lr_num, rr_num, this_rID, this_pos, len_read;
  string line, qname;
  bool this_isRC;
  stringstream ss;
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    string chrn = *ci;
    string f_input = fn_input + "." + chrn;
    string f_output = fn_output + "." + chrn;
    string f_output2 = fn_output2 + "." + chrn;
    list <int> pos_to_scan;    
    get_scan_pos(f_input, pos_to_scan, 8); // no big diff when change 8 to 20
    // cout << f_input << " " << pos_to_scan.size() << endl;
    list< READ_INFO *> lr_reads;    
    ifstream fin ( f_input.c_str());  // init, 700 * 30 /100 = 210 
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

void filter_location_rep(string file1, string file2, vector<string> &chrns, RepMaskPos *repmaskPos){
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    string chrn = *ci;
    ifstream fin( (file1 + "." + chrn).c_str());
    assert(fin);
    ofstream fout( (file2 + "." + chrn).c_str());
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

void sort_file_by_col(string fn, string col_name) {
  ifstream fin( fn.c_str());
  assert(fin); 
  string header, line, tmp1;
  int colv, coln = 0;
  bool col_name_exist = false;
  getline(fin, header);
  stringstream ss;
  ss.clear(); ss.str( header );
  while ( ss >> tmp1) { 
    coln++;
    if (tmp1 == col_name) {
      col_name_exist = true;
      break;
    }
  }
  assert (col_name_exist);
  //cout << "readline col " << coln << endl;
  list< IntString> rows_list;  
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    for (int i = 0; i < coln-1; i++) ss >> tmp1;
    ss >> colv;
    rows_list.push_back( make_pair(colv, line) );
  }
  fin.close();
  rows_list.sort(compare_IntString);
  ofstream fout( fn.c_str() );
  fout << header << endl;
  for (list< IntString>::iterator ri = rows_list.begin(); ri != rows_list.end(); ri++) 
    fout << ri->second << endl;
  fout.close();
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
  list< AlumateINFO *> alumate_list;
  read_match_rid(file1, alu_rid, alumate_list);  
  alumate_list.sort(AlumateINFO::sort_pos2);  
  AluRefPos *alurefpos = new AluRefPos(file_alu);                
  //////alurefpos->debug_print();
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
      fnOut_fnInput[fnOut].push_back(path0 + pn + "." + tmpn + "." + chrn);
      fnOut_pns[fnOut].push_back(pn);
    }
  }
  fnOut = (n_size <= group_size) ? fn_prefix : (fn_prefix + "." + int_to_string(n_group) + ".tmp0");
  for ( j = n_group * group_size; j < n_size; j++) {
    pn = id_pn_map[ j ];
    fnOut_fnInput[fnOut].push_back(path0 + pn + "." + tmpn + "." + chrn);	
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

void filter_pos_freq(string fn_input, string fn_output, float freq_min, float freq_max, int pn_cnt) {
  ifstream fin(fn_input.c_str());
  assert(fin);
  ofstream fout(fn_output.c_str());
  fout << "sum_LR_alumate regionBegin regionEnd pns\n";
  stringstream ss;
  string line, n_reads, refBegin, refEnd, pn;
  while (getline(fin, line)) {
    int n_pn = 0;
    ss.clear(); ss.str( line );
    ss >> n_reads >> refBegin >> refEnd;
    while ( ss >> pn) n_pn ++ ;
    float n_freq = (float)n_pn / pn_cnt; 
    if ( n_freq <= freq_max and n_freq >= freq_min) 
      fout << line << endl;
  }
  fin.close();
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

  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  chrns.push_back("chrX");
  chrns.push_back("chrY");

  string path0 = cf_fh.get_conf( "file_alu_insert0") ;    
  string path1 = cf_fh.get_conf( "file_alu_insert1") ;    
  string path_move;
  check_folder_exists(path0);
  check_folder_exists(path1);

  map<int, string> ID_pn;
  get_pn(cf_fh.get_conf( "file_pn"), ID_pn);
  
  if (opt == "write_tmps") {     

    int idx_pn = seqan::lexicalCast<int> (argv[3]);
    string pn = ID_pn[idx_pn];
    cout << "reading pn: " << idx_pn << " " << pn << "..................\n";
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, "");
    string file1_prefix = path0 + pn + ".tmp1";

//    //// step1: write pn.tmp1.
//    map <string, int> rg_to_idx;    
//    parse_reading_group( get_name_rg(cf_fh.get_conf("file_dist_prefix"), pn), rg_to_idx );    
//    string header1 = "qname A_chr B_chr A_beginPos B_beginPos A_isRC B_isRC len_read rgIdx\n";
//    alu_mate_flag(bam_fh, file1_prefix, header1, rg_to_idx);
//    cout << "done with reading this bam file, time used " << clocki.elapsed() << endl; // about 2 hours    
//    //// step2: write pn.tmp1.chr*.  use database to filter alu mate
//    string header2 = "qname this_rID this_pos this_isRC len_read alu_rID alu_pos alu_isRC rgIdx alu_type\n";
//    MapFO fileMap;    
//    map<int, string>::iterator rc ;
//    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) {
//      fileMap[rc->first] = new ofstream( (file1_prefix + "." + rc->second).c_str() );  // write down reads whose pair is mapped to Alu
//      assert(fileMap[rc->first]);
//      *(fileMap[rc->first]) << header2 ;
//    }     
//    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) 
//      keep_alu_mate(bam_fh, rc->first,  file1_prefix, fileMap, cf_fh.get_conf( "file_alupos_prefix") + rc->second);           
//    close_fhs(fileMap, bam_fh->rID_chrn);    
//    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) 
//      sort_file_by_col(file1_prefix + "." + rc->second, "this_pos");
//    move_files(path0+"tmp1s/", file1_prefix);    
    
#ifdef DEBUG_MODE   // only look at chr1 for now 
    cerr << "!!!! DEBUG_MODE: only chr1 tested\n";
    chrns.clear();
    chrns.push_back("chr1");
#endif       

    string file2_prefix =  path0 + pn + ".tmp2";    // left and right counts at each check position
    string file3_prefix =  path0 + pn + ".tmp3";    // combine *tmp2 files into regions
    alumate_counts_filter(file1_prefix, file2_prefix, file3_prefix, chrns, LEFT_PLUS_RIGHT); 
    move_files(path0+"tmp1s/", path0 + pn + ".tmp1.chr*") ;
    move_files(path0+"tmp2s/", path0 + pn + ".tmp2.chr*") ;
    //// filter out rep regions,  combine rep regions(dist < REP_MASK_JOIN bp)
    RepMaskPos *repmaskPos;
    repmaskPos = new RepMaskPos(cf_fh.get_conf( "file_repeatMask"), chrns, REP_MASK_JOIN); 
    string file3st_prefix = path0 + pn + ".tmp3st";  // subset of *tmp3
    filter_location_rep(file3_prefix, file3st_prefix, chrns, repmaskPos);    
    delete bam_fh;
    delete repmaskPos;
    
    }  else if (opt == "combine_pos") { // combine positions from multiple individuals

#ifdef DEBUG_MODE   // only look at chr1 for now 
    chrns.clear();
    chrns.push_back("chr1");
#endif               

    string file_pn_used = cf_fh.get_conf( "file_pn_used");
    ifstream fin(file_pn_used.c_str());
    if (!fin) {
      cerr << "which individuals should be used for further analysis?\n";
      cerr << "this information should be in file " << file_pn_used << endl;
      return 1;
    }

    assert(fin);
    map<int, string> id_pn_map;   
    int i = 0;
    string pn;
    while (fin >> pn) id_pn_map[i++] = pn;
    fin.close();              
    cout << id_pn_map.size() << " individuals out of " << ID_pn.size() << " individuals are used, change "
	 << file_pn_used << " if you want to modify\n";

    string output_prefix = path1 + "insert_pos.";
    float freq_min = seqan::lexicalCast<float> (cf_fh.get_conf("freq_min"));
    float freq_max = seqan::lexicalCast<float> (cf_fh.get_conf("freq_max"));
    string fn_suffix = get_name_suffix(freq_min, freq_max);   
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      combine_pn_pos(*ci, id_pn_map, "tmp3st", output_prefix, 10, path0);  // takes about 12 mins to run !   
      system( ("rm " + output_prefix + *ci + "*tmp?").c_str()  );      
      string fn_input = path1 + "insert_pos."+ *ci;
      filter_pos_freq( fn_input, fn_input + fn_suffix, freq_min, freq_max, id_pn_map.size() );
    }        
  }
  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
