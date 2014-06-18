#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include "delete_utils.h"

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

bool clipreads_at_insertPos(string pn, string chrn, BamFileHandler *bam_fh, FastaFileHandler *fasta_fh,  string fin_pos, string fout_reads) {
  ifstream fin( (fin_pos).c_str());
  if (!fin) return false;
  string line, numAluRead, _pn;
  getline(fin, line); // skip header
  ofstream fout(fout_reads.c_str());
  fout << "pn region_begin region_end left_right adj_clipPos qName clipPos cigar tLen seq\n" ;
  stringstream ss, ss_header;
  int region_begin, region_end, refBegin, refEnd;
  
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> numAluRead >> region_begin >> region_end;
    bool pn_has_alumate = false;
    while ( ss >> _pn ) 
      if ( _pn == pn) { pn_has_alumate=true; break;}
    if (!pn_has_alumate) continue;

    refBegin = region_begin - 2 *  DEFAULT_READ_LEN;
    refEnd = region_end + 2 * DEFAULT_READ_LEN;    
    if (! bam_fh->jump_to_region(chrn, refBegin, refEnd) ) continue;
    ss_header.clear(); ss_header.str("");
    ss_header << pn << " " << region_begin << " " << region_end << " " ;    

    seqan::BamAlignmentRecord record;
    while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, refBegin, refEnd, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_delete_read(record)) continue;  // search for clip reads 
      int beginPos = record.beginPos;
      int endPos = beginPos + getAlignmentLengthInRef(record);
      int adj_clipPos, align_len;
      seqan::CharString ref_fa;
      list <char> cigar_opts;
      list <int> cigar_cnts;
      
      if (has_soft_last(record, CLIP_BP) and endPos >= region_begin - FLANK_REGION_LEN and endPos < region_end + FLANK_REGION_LEN) {
	adj_clipPos = endPos;
	if (length(ref_fa) == 0) {
	  fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
	  parse_cigar(get_cigar(record), cigar_opts, cigar_cnts);
	}
	if (clipLeft_move_right( record.seq, ref_fa, cigar_cnts, refBegin, adj_clipPos, align_len) and 
	    align_clip_to_ref('L', adj_clipPos, endPos, align_len, record, fasta_fh, fout, ss_header.str()) )
	  continue;
      }
      
      if ( has_soft_first(record, CLIP_BP) and beginPos >= region_begin - FLANK_REGION_LEN and beginPos < region_end + FLANK_REGION_LEN) {
	adj_clipPos = beginPos;
	if (length(ref_fa) == 0) {
	  fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
	  parse_cigar(get_cigar(record), cigar_opts, cigar_cnts);
	}
	if (clipRight_move_left( record.seq, ref_fa, cigar_cnts, refBegin, adj_clipPos, align_len) and 
	    align_clip_to_ref('R', adj_clipPos, beginPos, align_len, record, fasta_fh, fout, ss_header.str()) )
	  continue;
      }
    }
  }
  fin.close();
  fout.close();
  return true;
}

// for each insertion region, combine reads from different pns 
void combine_clipreads_by_pos(std::set<string> &pns_used, ConfigFileHandler & cf_fh, string fn_suffix, string chrn) {
  map< pair<string, string>, vector<string> > regionPos_lines;
  for (std::set<string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) {
    string file_input_clip = cf_fh.get_conf( "file_clip_reads") + chrn + "/" + *pi + fn_suffix;
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
    string file_output_clip = cf_fh.get_conf( "file_clip_reads") + chrn + "_pos/" + (ri->first).first + "_" + (ri->first).second;	
    ofstream fout(file_output_clip.c_str() );
    for (vector<string>::iterator ii = (ri->second).begin(); ii != (ri->second).end(); ii++) 
      fout << *ii << endl;
    fout.close();
  }
}

void parse_fn(string fn, int & regionBegin, int & regionEnd) {
  stringstream ss;
  ss.str(fn);
  string token;
  getline(ss, token, '_');
  seqan::lexicalCast2(regionBegin, token);
  getline(ss, token, '_');
  seqan::lexicalCast2(regionEnd, token);
}

bool nonempty_files_sorted(string & path1, map < pair<int, int>, string> & regionPos_fn) {
  DIR *d;
  struct dirent *dir;
  d = opendir(path1.c_str());
  string fn;
  int regionBegin, regionEnd;
  if (d) {
    while ( (dir = readdir(d))!=NULL )  {
      fn = dir->d_name;
      if (fn[0] == '.' or !is_nonempty_file(fn)) continue;
      parse_fn(fn, regionBegin, regionEnd);
      regionPos_fn[ make_pair(regionBegin, regionEnd)] = fn;
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

// each pn can have max 2 pos
bool pn_vote_max2pos( vector <pair<char, int> > & pos_of_pn, int & p1, int & p2, float & f1, float &f2){
  p1 = p2 = 0;  
  f1 = f2 = 0;
  int vsize = pos_of_pn.size();
  if ( vsize < 1) return false;
  map <int, float> clipLeft_cnt;
  for (vector <pair<char, int> >::iterator pi = pos_of_pn.begin(); pi != pos_of_pn.end(); pi ++ ) {
    int match_offset =  ((*pi).first == 'L') ? CLIP_BP_LEFT : CLIP_BP_RIGHT ;
    addcnt_if_pos_match(clipLeft_cnt, (*pi).second, match_offset, 1);
  }

  //debugprint_map(clipLeft_cnt);
  multimap<float, int> cnt_clipLeft = flip_map(clipLeft_cnt);
  multimap<float, int>::reverse_iterator li = cnt_clipLeft.rbegin();
  p1 = (li->second);
  f1 = (li->first) / vsize;
  if ( cnt_clipLeft.size() == 1) return true;
  li ++; 
  if ( (f2 = (li->first) / vsize) > 0.3) { // 2nd vote
    p2 = (li->second);
    float f12 = f1 + f2;
    f1 /= f12;
    f2 /= f12;
  } else {
    f1 = 1;
    f2 = 0;
  }
  return true;
}

void region_pos_vote( string fn, int & clipLeft1, int & clipLeft2 ) {
  ifstream fin(fn.c_str());
  assert(fin);
  stringstream ss;
  string line, pn, sleft_right;
  int region_begin, region_end, clipPos;
  map <string, vector <pair<char, int> > > pn_splitori_pos;
  while ( getline(fin, line)) {  
    ss.clear(); ss.str( line );
    ss >> pn >> region_begin >> region_end >> sleft_right >> clipPos;
    pn_splitori_pos[pn].push_back( make_pair( sleft_right[0], clipPos) );
  }
  fin.close();  
  map < int, float> clipLeft_cnt;
  float vsize = 0;  
  for (map <string, vector <pair<char, int> > >::iterator pi = pn_splitori_pos.begin(); pi != pn_splitori_pos.end(); pi++) {
    int p1, p2;
    float f1, f2;
    if (!pn_vote_max2pos(pi->second, p1, p2, f1, f2)) 
      continue;
    //cout << pi->first << " " << p1 << " " << p2 << " " << f1 << " " << f2 << endl;
    assert ( abs (f1 + f2 - 1.) < 1e-5) ;
    // fixme : add weight to the count ?
    if (p1 > 0) vsize += 1;
    if (p1 > 0) addcnt_if_pos_match(clipLeft_cnt, p1, CLIP_BP_RIGHT, f1);
    if (p2 > 0) addcnt_if_pos_match(clipLeft_cnt, p2, CLIP_BP_RIGHT, f2);
  }
  multimap<float, int> cnt_clipLeft = flip_map(clipLeft_cnt);
  multimap<float, int>::reverse_iterator li = cnt_clipLeft.rbegin();
  clipLeft1 = 0, clipLeft2 = 0 ;
  if ( li->first / vsize >= MIN_VOTE_FREQ) clipLeft1 = li->second;
  if ( cnt_clipLeft.size() > 1) {
    li++;
    if ( li->first / vsize >= MIN_VOTE_FREQ) clipLeft2 = li->second;
  }
}

void regions_pos_vote( string path1, string fn_output) {
  map < pair<int,int>,  string> regionPos_fn;
  assert( nonempty_files_sorted( path1, regionPos_fn) );
  ofstream fout(fn_output.c_str());
  fout << "region_begin region_end clipLeft\n";
  for (map < pair <int, int>, string >::iterator ri = regionPos_fn.begin(); ri != regionPos_fn.end(); ri++ ) {
    int clipLeft1, clipLeft2 ;
    //if ( (ri->first).first > 26984) break;
    //cout << path1 + ri->second << endl;
    region_pos_vote( path1 + ri->second, clipLeft1, clipLeft2);
    if (clipLeft1 > 0) {
      fout << (ri->first).first << " " << (ri->first).second << " " << clipLeft1 ;
      if ( clipLeft2 > 0 )
	fout << " " << clipLeft2 ;
      fout << endl;
    }
  }  
}

void exactpos_pns(string fn_input, string path1, string fn_output){
  ifstream fin;
  fin.open(fn_input.c_str());
  assert(fin);
  map< int, vector<string> > clipLeft_fn;
  stringstream ss;
  string line, region_begin, region_end;
  int clipLeft;
  getline(fin,line);
  while (getline(fin,line)) {
    ss.clear(); ss.str( line );
    ss >> region_begin >> region_end;
    while ( ss >> clipLeft) 
      clipLeft_fn[clipLeft].push_back(region_begin+"_"+region_end);
  }
  fin.close();
  cout << "NB: private insertion is ignored !\n";
  ofstream fout(fn_output.c_str());
  fout << "clipLeft clipRight pn(count)\n";
  for (map< int, vector<string> >::iterator pi = clipLeft_fn.begin(); pi != clipLeft_fn.end(); pi++ ) {
    int clipLeft = pi->first;
    map <string, std::set<string> > pn_qName;
    for ( vector<string>::iterator fi = (pi->second).begin(); fi != (pi->second).end(); fi++ ) {
      fin.open( (path1 + *fi).c_str());
      assert(fin);
      string line, pn, tmp1, tmp2, sleft_right, qName;
      stringstream ss;
      int clipPos;
      while (getline(fin, line)) {
	ss.clear(); ss.str( line );
	ss >> pn >> tmp1 >> tmp2 >> sleft_right >> clipPos >> qName;
	int match_offset =  ( sleft_right[0] == 'L') ? CLIP_BP_LEFT : CLIP_BP_RIGHT ;
	if ( abs(clipPos - clipLeft ) <= match_offset )
	  pn_qName[pn].insert(qName);
      }
      fin.close();
    }
    if (pn_qName.size() <= 1) continue;
    fout << clipLeft << " -1" ;  // clipRight unknown for now
    for (map <string, std::set<string> >::iterator pi = pn_qName.begin(); pi != pn_qName.end(); pi++) 
      fout << " " << pi->first << "," << (pi->second).size() ;
    fout << endl;
  }
  fout.close();
}

void clip_skip_unknow_reads(string pn, string chrn, vector< pair<int, int> > & insert_pos, BamFileHandler *bam_fh, string file_output,  map < int, vector<ALUREAD_INFO> > & rid_alureads, map < int, vector < string >  > & insertBegin_unknowreads, map < int, int > & insertBegin_skipreads,  map < int, int > & insertBegin_midreads, map<string, int> & rg_to_idx, int alucons_len) {

  ofstream fout(file_output.c_str());    
  fout << "clipPos seq_RC_correct read_type sameRC qName\n"; // col2, some seqs are RC according to sameRC flag
  seqan::BamAlignmentRecord record;
  string read_status;
  for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
    int clipLeft = (*pi).first;
    int clipRight = (*pi).second;
    if (clipRight < clipLeft) clipRight = clipLeft;    
    int offset_left = CLIP_BP_LEFT;
    int offset_right = ( (*pi).second < 0 ) ? CLIP_BP_RIGHT : CLIP_BP_LEFT;    
    int region_begin = clipLeft - ALU_FLANK; // no need for larger flank region, if i only want to build consensus
    int region_end = clipRight + ALU_FLANK;	
    if (! bam_fh->jump_to_region(chrn, region_begin, region_end)) 
      continue;    

    map < seqan::CharString, string > rg_str;
    std::set < seqan::CharString > clipreads_qname;
    std::set < seqan::CharString > skipreads_qname;
    while ( true ) {
      read_status = bam_fh->fetch_a_read(chrn, region_begin, region_end, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_insert_read(record)) continue;
      int thisEnd = record.beginPos + (int)getAlignmentLengthInRef(record);   
      if (record.rID != record.rNextId or abs(record.tLen) >= DISCORDANT_LEN) {
	if ( ( thisEnd < clipLeft and !hasFlagRC(record) ) or 
	     ( record.beginPos > clipRight and hasFlagRC(record) ) ) {
	  ALUREAD_INFO one_aluRead = ALUREAD_INFO(record.qName, clipLeft, record.pNext, pn, hasFlagRC(record)==hasFlagNextRC(record) );
	  rid_alureads[record.rNextId].push_back(one_aluRead);
	}
	continue;
      }      
      // otherwise only look at proper mapped reads 
      bool read_is_left = left_read(record);
      if ( hasFlagRC(record) == hasFlagNextRC(record) or read_is_left == hasFlagRC(record) ) continue;       
      // clip reads
      if ( ( has_soft_first(record, CLIP_BP) and abs( record.beginPos - clipRight) <=  offset_right) or
	   ( has_soft_last(record, CLIP_BP) and abs( thisEnd - clipLeft ) <= offset_left ) ) {  	    
	fout << clipLeft << " " << record.seq  <<  " clipRead -1 " << record.qName << endl;
	clipreads_qname.insert(record.qName);
	continue;
      }      
      // skip reads, fixme: add align_aluSkip(record.seq, fasta_fh, length(record.seq) ??? in ins_del.cpp
      if ( (record.beginPos) < clipLeft - CLIP_BP and thisEnd > clipLeft + CLIP_BP and count_non_match(record) <= 5 ) {
	skipreads_qname.insert(record.qName);
	continue;
      }
      // unknown read
      if ( (read_is_left and (record.beginPos) < clipLeft and record.pNext > clipLeft ) or 
	   (!read_is_left and (record.pNext) < clipLeft and thisEnd > clipLeft ) ) {
	int rgIdx = get_rgIdx(rg_to_idx, record);
	stringstream rg_ss;
	rg_ss << " " << rgIdx << ":" << abs(record.tLen) + alucons_len;
	rg_str[record.qName] = rg_ss.str();	
      }
    }
    int err_cnt = 0; 
    if ( !clipreads_qname.empty() and !skipreads_qname.empty()) 
      for ( std::set < seqan::CharString >::iterator si = skipreads_qname.begin(); si != skipreads_qname.end(); si++) {
	std::set < seqan::CharString >::iterator ci = clipreads_qname.find( *si );
	if ( ci != clipreads_qname.end()) {
	  clipreads_qname.erase(ci);  // more sure about skipreads
	  err_cnt += 1;
	}
      }
    //if (err_cnt) cout << err_cnt << " err reads are discarded\n"; 
    insertBegin_skipreads[clipLeft] = skipreads_qname.size();
    insertBegin_midreads[clipLeft] = clipreads_qname.size();
    
    for ( map < seqan::CharString, string >::iterator ri = rg_str.begin(); ri != rg_str.end(); ri++) 
      if ( clipreads_qname.find(ri->first) == clipreads_qname.end() and skipreads_qname.find(ri->first) == skipreads_qname.end() ) 
	insertBegin_unknowreads[clipLeft].push_back( ri->second );
    rg_str.clear(); 
    
    //cout << clipLeft  << "  "  <<  clipreads_qname.size() << endl;

  }
  fout.close();
}

void add_aluReads(BamFileHandler *bam_fh, string file_cons,  map < int, vector<ALUREAD_INFO> > & rid_alureads, AluconsHandler *alucons_fh, map < int, int > & insertBegin_midreads) {
  fstream fout;
  map < string, bool > align_alucons_pass;
  fout.open(file_cons.c_str(), fstream::app|fstream::out);
  seqan::BamAlignmentRecord record;
  for ( map < int, vector<ALUREAD_INFO> >::iterator ri = rid_alureads.begin(); ri != rid_alureads.end(); ri++) {
    for ( vector<ALUREAD_INFO>::iterator ai = (ri->second).begin(); ai != (ri->second).end(); ai++ ) {
      if ( !bam_fh->jump_to_region( ri->first, (*ai).pos - 20, (*ai).pos + 20) ) continue;
      while (bam_fh->fetch_a_read(record) and record.beginPos <= (*ai).pos + 3) {
	if (record.qName == (*ai).qName) {
	  string record_seq = toCString(record.seq);
	  // fixme: why ok at alu_insdel "reverseComplement(record_seq)" without seqan:: ?
	  if ((*ai).sameRC)  seqan::reverseComplement(record_seq); 
	  map < string, bool >::iterator api = align_alucons_pass.find(record_seq);
	  if ( api != align_alucons_pass.end() ) {
	    if ( api->second) {
	      fout << (*ai).clipLeft << " " << record_seq << " aluRead " << (*ai).sameRC << " " <<  record.qName <<  endl;	      
	      addKey(insertBegin_midreads, (*ai).clipLeft, 1);
	    }
	  } else {
	    float sim_rate;
	    bool aa_pass = align_alu_cons_call( record_seq, alucons_fh, sim_rate, ALUCONS_SIMILARITY);
	    align_alucons_pass[record_seq] = aa_pass;
	    if (aa_pass) {
	      fout << (*ai).clipLeft << " " << record_seq << " aluRead " << (*ai).sameRC << " " <<  record.qName <<  endl;	      
	      addKey(insertBegin_midreads, (*ai).clipLeft, 1);
	    //cout << "add " << (*ai).clipLeft << endl;
	    }
	  }
	  break;
	}
      }  // this pos, done
    }   
    //cout << ri->first << " " << (ri->second).size() << endl;
    (ri->second).clear();
  }
  fout.close();
}
  
void write_tmp1(string chrn, vector< pair<int, int> > & insert_pos, ofstream & fout, map < int, vector< string> > & insertBegin_unknowreads, map < int, int > & insertBegin_skipreads,  map < int, int > & insertBegin_midreads,  int alucons_len) {
  for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
    int clipLeft = (*pi).first;
    int clipRight = (*pi).second;
    if (clipRight < clipLeft) clipRight = clipLeft;
    int midCnt = 0, clipCnt = 0;
    string unknowStr;
    get_mapVal(insertBegin_midreads, clipLeft, midCnt);
    get_mapVal(insertBegin_skipreads, clipLeft, clipCnt);
    map < int, vector< string> >::iterator ui = insertBegin_unknowreads.find(clipLeft);
    int unknowCnt = ( ui == insertBegin_unknowreads.end() ? 0 : (ui->second).size() );    
    if ( !midCnt and !unknowCnt  ) continue;
    fout << chrn << " " << clipLeft << " " << clipRight << " " << alucons_len  << " " << midCnt << " " << clipCnt << " " << unknowCnt;
    if (unknowCnt) 
      for ( vector< string > ::iterator ui2 = (ui->second).begin(); ui2 != (ui->second).end(); ui2 ++ ) fout << *ui2 ;	
    fout << endl;    
  }
}

bool parseline_del_tmp1(string &line, string & output_line, map <int, EmpiricalPdf *> & pdf_rg){
  float *log10_gp = new float[3];
  stringstream ss, ss_out;
  string chrn;
  int clipLeft, clipRight, estimatedAluLen, midCnt, clipCnt, unknowCnt;
  ss.clear(); ss.str(line); 
  ss >> chrn >> clipLeft >> clipRight >> estimatedAluLen >> midCnt >> clipCnt >> unknowCnt ;
  float prob_ub = pow(10, -LOG10_RATIO_UB);
  float prob_known = (midCnt+clipCnt)/(float)(midCnt + clipCnt + unknowCnt);
  for (int i = 0; i < 3; i++) log10_gp[i] = 0;
  if (prob_known) {
    log10_gp[0] = clipCnt * log10 ( prob_known * prob_ub ) + midCnt * log10 ( prob_known * (1 - prob_ub) );
    log10_gp[1] = (midCnt + clipCnt) * log10 (prob_known * 0.5) ; 
    log10_gp[2] = midCnt * log10 ( prob_known * prob_ub ) + clipCnt * log10 ( prob_known * (1 - prob_ub) );
  }
  if (unknowCnt) { 
    int insert_len, idx;
    string token;
    for (int i = 0; i < unknowCnt; i++) {
      getline(ss, token, ':');
      seqan::lexicalCast2(idx, token);
      getline(ss, token, ' ');
      seqan::lexicalCast2(insert_len, token);      
      float p_y = pdf_rg[idx]->pdf_obs(insert_len);
      float p_z = pdf_rg[idx]->pdf_obs(insert_len - estimatedAluLen);
      //float freq0 = 0.67;  // high FP ratio      
      float freq0 = ( midCnt + 1 )/(float)(midCnt + clipCnt + 2); // 1 and 2 are psudo count
      log10_gp[0] += log10 (p_y * (1 - prob_known));
      log10_gp[1] += log10 ((freq0 * p_y + (1 - freq0) * p_z) * (1 - prob_known) ) ;
      log10_gp[2] += log10 (p_z * (1 - prob_known));
    }
  }
  bool use_this_line = false;
  if ( !p11_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    ss_out << chrn << " " << clipLeft << " " << estimatedAluLen << " " << midCnt << " " << clipCnt << " " << unknowCnt ;
    float *gp = new float[3];
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
    ss_out << " " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2]; 
    delete gp;    
    output_line = ss_out.str();
    use_this_line = true;
  }
  delete log10_gp;
  return use_this_line;
}

void write_tmp2(string fn_tmp1, string fn_tmp2, map <int, EmpiricalPdf *> & pdf_rg){
  ofstream fout(fn_tmp2.c_str());
  fout << "chr insertBegin insertLen midCnt clipCnt unknowCnt 00 01 11\n";  
  string line, output_line;
  ifstream fin(fn_tmp1.c_str());
  assert(fin);
  getline(fin, line); // read header
  while (getline(fin, line))
    if (parseline_del_tmp1(line, output_line, pdf_rg))
      fout << output_line << endl;
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
  
#ifdef DEBUG_MODE   // only look at chr1 for now 
  cerr << "!!!! DEBUG_MODE: only chr1 tested\n";
  chrns.clear();
  chrns.push_back("chr1");
#endif       

  string path0 = cf_fh.get_conf( "file_alu_insert0") ;    
  string path1 = cf_fh.get_conf( "file_alu_insert1") ;    
  string path_move;
  check_folder_exists(path0);
  check_folder_exists(path1);
  map <int, string> ID_pn;
  get_pn(cf_fh.get_conf( "file_pn"), ID_pn);
  
  float freq_min = seqan::lexicalCast<float> (cf_fh.get_conf("freq_min"));
  float freq_max = seqan::lexicalCast<float> (cf_fh.get_conf("freq_max"));
  string fn_suffix = get_name_suffix(freq_min, freq_max);   
  
  if (opt == "write_tmps_pn") {     
    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input);
    string file1_prefix = path0 + pn + ".tmp1";
    map<int, string>::iterator rc ;

    //// step1: write pn.tmp1.
    map <string, int> rg_to_idx;    
    parse_reading_group( get_name_rg(cf_fh.get_conf("file_dist_prefix"), pn), rg_to_idx );    
    string header1 = "qname A_chr B_chr A_beginPos B_beginPos A_isRC B_isRC len_read rgIdx\n";
    alu_mate_flag(bam_fh, file1_prefix, header1, rg_to_idx);
    cout << "done with reading this bam file, time used " << clocki.elapsed() << endl; // about 2 hours    
    //// step2: write pn.tmp1.chr*.  use database to filter alu mate
    string header2 = "qname this_rID this_pos this_isRC len_read alu_rID alu_pos alu_isRC rgIdx alu_type\n";
    MapFO fileMap;    
    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) {
      fileMap[rc->first] = new ofstream( (file1_prefix + "." + rc->second).c_str() );  // write down reads whose pair is mapped to Alu
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << header2 ;
    }     
    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) 
      keep_alu_mate(bam_fh, rc->first,  file1_prefix, fileMap, cf_fh.get_conf( "file_alupos_prefix") + rc->second);           
    close_fhs(fileMap, bam_fh->rID_chrn);    
    
    int col_idx = 0;
    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) {
      string file1 = file1_prefix + "." + rc->second;
      if (!col_idx) col_idx = get_col_idx(file1, "this_pos"); 
      sort_file_by_col<int> (file1, col_idx, true);
    }
    move_files(path0+"tmp1s/", file1_prefix);    
    
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
    
  }  else if (opt == "combine_pos_pns") { // combine positions from multiple individuals
    string file_pn_used = cf_fh.get_conf( "file_pn_used");
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
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      combine_pn_pos(*ci, id_pn_map, "tmp3st", output_prefix, 10, path0);  // takes about 12 mins to run !   
      system( ("rm " + output_prefix + *ci + "*tmp?").c_str()  );      
      string fn_input = path1 + "insert_pos."+ *ci;
      filter_pos_freq( fn_input, fn_input + fn_suffix, freq_min, freq_max, id_pn_map.size() );
    }        
    
  } else if (opt == "clipReads_by_pn") { // clip reads for each pn
    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];
    
    std::set <string> pns_used;
    read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
    if ( pns_used.find(pn) == pns_used.end() ){
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    string fout_path = cf_fh.get_conf( "file_clip_reads");
    check_folder_exists(fout_path);     
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bam_input+".bai");
    string file_fa = cf_fh.get_conf("file_fa_prefix");    
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {      
      string fin_pos = path1 + "insert_pos." + *ci + fn_suffix;      
      string tmp_path = fout_path + *ci + "/";
      check_folder_exists( tmp_path);
      string fout_reads = tmp_path + pn + fn_suffix;
      FastaFileHandler * fasta_fh = new FastaFileHandler(file_fa + *ci + ".fa", *ci);
      clipreads_at_insertPos(pn, *ci, bam_fh, fasta_fh, fin_pos, fout_reads);
      delete fasta_fh;
    }
    delete bam_fh;
    cout << "output to " << fout_path << "chr*/" << pn << fn_suffix << endl;
    
  } else if ( opt == "clipReads_pos_pns" ) { 

    std::set <string> pns_used;
    read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
    string fout_path = cf_fh.get_conf( "file_clip_reads");
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      combine_clipreads_by_pos(pns_used, cf_fh, fn_suffix, *ci);
      string tmp_path = fout_path + *ci + "_pos/";
      check_folder_exists( tmp_path );      
      string file_insert_pos = path1 + "clip/" + *ci + fn_suffix;
      regions_pos_vote(tmp_path, file_insert_pos + ".tmp");
      exactpos_pns(file_insert_pos + ".tmp", tmp_path, file_insert_pos);
    }
  } else if ( opt == "fixed_delete0_pn" ) {   // write fasta reads for consensus, use clip reads first 
    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];

    std::set <string> pns_used;
    read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
    if ( pns_used.find(pn) == pns_used.end() ){ // a bit slow
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bam_input+".bai");
    AluconsHandler *alucons_fh = new AluconsHandler(cf_fh.get_conf("file_alu_cons"), cf_fh.get_conf("type_alu_cons"));    
    map<string, int> rg_to_idx;
    string file_dist_prefix = cf_fh.get_conf("file_dist_prefix");
    parse_reading_group( get_name_rg(file_dist_prefix, pn), rg_to_idx );
    
    string path_cons = cf_fh.get_conf("file_ins_cons");  // for consensus sequence
    check_folder_exists( path_cons);
    string path_del0 = cf_fh.get_conf("file_ins_fixed_del0");  // write tmp0 file
    check_folder_exists( path_del0);    
    string file_tmp1 = path_del0 + pn + ".tmp1";
    string file_tmp2 = path_del0 + pn + ".tmp2";    
    ofstream fout1 ( file_tmp1.c_str());
    fout1 << "chr insertBegin insertEnd estimatedAluLen midCnt clipCnt unknowCnt unknowStr\n";
    for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
      check_folder_exists( path_cons + *ci);
      string file_output = path_cons + *ci + "/" + pn;
      string file_clipPos = path1 + "clip/" + *ci + fn_suffix;
      vector< pair<int, int> > insert_pos;      
      read_first2col( file_clipPos , insert_pos, true);   // debug, read first 10 col!!
      map < int, vector<ALUREAD_INFO>  > rid_alureads;       // info of alu reads 
      map < int, vector<string> > insertBegin_unknowreads;   
      map < int, int > insertBegin_skipreads;   
      map < int, int > insertBegin_midreads;   // midreads (inserting a alu seq) = clip + alu reads 
      clip_skip_unknow_reads(pn, *ci, insert_pos, bam_fh, file_output, rid_alureads, insertBegin_unknowreads, insertBegin_skipreads, insertBegin_midreads, rg_to_idx, alucons_fh->seq_len);      
      add_aluReads(bam_fh, file_output, rid_alureads, alucons_fh, insertBegin_midreads); // get seq from bam instead of genome.fa      
      sort_file_by_col<int> (file_output, 1, true);  // sort by clip pos
      rid_alureads.clear();      
      write_tmp1(*ci, insert_pos, fout1, insertBegin_unknowreads,  insertBegin_skipreads, insertBegin_midreads, alucons_fh->seq_len);
    }   
    fout1.close();
    map <int, EmpiricalPdf *> pdf_rg;    
    string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
    read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
    write_tmp2(file_tmp1, file_tmp2, pdf_rg);
    EmpiricalPdf::delete_map(pdf_rg);
    move_files(path_del0 + "tmp1s/", file_tmp1);
    move_files(path_del0 + "tmp2s/", file_tmp2);

  } else if (opt == "fixed_vcf_pns") {
    string path_del0 = cf_fh.get_conf("file_ins_fixed_del0") ;
    vector <string> pns;
    read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns); 
    string path_input = path_del0 + "tmp2s/";
    
    string tmp_file_pn = path_input + *(pns.begin()) + ".tmp2";
    int col_idx =  get_col_idx(tmp_file_pn, "00");
    assert (col_idx == 7 );

    string fn_pos = path_del0 + int_to_string( pns.size()) + ".pos";
    filter_by_llh_noPrivate(path_input, ".tmp2", fn_pos, pns, chrns, col_idx);
    string fn_vcf = path_del0 + int_to_string( pns.size()) + ".vcf";  
    combine_pns_vcf_noPrivate(path_input, ".tmp2", fn_vcf, pns, chrns, col_idx);  

  } else {
    cout << "unknown option ! \n";
  }

  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
