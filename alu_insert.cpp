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

bool clipreads_at_insertPos(string pn, string chrn, BamFileHandler *bam_fh, FastaFileHandler *fasta_fh,  string fin_pos, string fout_reads) {
  const int FLANK_REGION_LEN = 80;
  ifstream fin( (fin_pos).c_str());
  if (!fin) return false;
  string line, numAluRead, _pn;
  ofstream fout(fout_reads.c_str());
  fout << "pn region_begin region_end left_right adj_clipPos qName clipPos cigar\n" ;
  stringstream ss, ss_header;
  int region_begin, region_end, refBegin, refEnd;
  
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> numAluRead >> region_begin >> region_end;
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

void fn_to_pos(string fn, int & regionBegin, int & regionEnd) {
  stringstream ss;
  ss.str(fn);
  string token;
  getline(ss, token, '_');
  seqan::lexicalCast2(regionBegin, token);
  getline(ss, token, '_');
  seqan::lexicalCast2(regionEnd, token);
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
  if (d) {
    while ( (dir = readdir(d))!=NULL )  {
      fn = dir->d_name;      
      if (fn[0] == '.' or check_file_size( path1 + fn)<=0 ) continue;
      fn_to_pos(fn, regionBegin, regionEnd);
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

int region_pos_vote( string fn, int & clipLeft1, int & clipLeft2 ) {
  const float MIN_VOTE_FREQ = 0.4;
  clipLeft1 = 0, clipLeft2 = 0 ;
  ifstream fin(fn.c_str());
  assert(fin);
  stringstream ss;
  string line, pn, sleft_right;
  int region_begin, region_end, clipPos, cnt = 0;
  map <string, vector <pair<char, int> > > pn_splitori_pos;
  while ( getline(fin, line)) {  
    ss.clear(); ss.str( line );
    ss >> pn >> region_begin >> region_end >> sleft_right >> clipPos;
    pn_splitori_pos[pn].push_back( make_pair( sleft_right[0], clipPos) );
    cnt ++;
  }
  fin.close();  
  map < int, float> clipLeft_cnt;
  float vsize = 0;  
  // each pn vote max two positions,  fixme : add weight to the count ?
  for (map <string, vector <pair<char, int> > >::iterator pi = pn_splitori_pos.begin(); pi != pn_splitori_pos.end(); pi++) {
    int p1, p2;
    float f1, f2;
    if (!pn_vote_max2pos(pi->second, p1, p2, f1, f2)) 
      continue;
    //cout << pi->first << " " << p1 << " " << p2 << " " << f1 << " " << f2 << endl;
    assert ( abs (f1 + f2 - 1.) < 1e-5) ;
    if (p1 > 0) vsize += 1;
    if (p1 > 0) addcnt_if_pos_match(clipLeft_cnt, p1, CLIP_BP_RIGHT, f1);
    if (p2 > 0) addcnt_if_pos_match(clipLeft_cnt, p2, CLIP_BP_RIGHT, f2);
  }
  multimap<float, int> cnt_clipLeft = flip_map(clipLeft_cnt);
  multimap<float, int>::reverse_iterator li = cnt_clipLeft.rbegin();
  if ( li->first / vsize >= MIN_VOTE_FREQ) clipLeft1 = li->second;
  if ( cnt_clipLeft.size() > 1) {
    li++;
    if ( li->first / vsize >= MIN_VOTE_FREQ) clipLeft2 = li->second;
  }
  return 0;
}

void regions_pos_vote( string path1, string fn_output, std::set < pair<int, int> > & regionPos_set) {
  ofstream fout(fn_output.c_str());
  fout << "region_begin region_end clipLeft\n";
  for (std::set < pair <int, int> >::iterator ri = regionPos_set.begin(); ri != regionPos_set.end(); ri++ ) {
    int clipLeft1, clipLeft2 ;
    string tmp_file = path1 + pos_to_fn(*ri);
    region_pos_vote( tmp_file, clipLeft1, clipLeft2);
    if (clipLeft1 > 0) {
      fout << (*ri).first << " " << (*ri).second << " " << clipLeft1 ;
      if ( clipLeft2 > 0 )
	fout << " " << clipLeft2 ;
      fout << endl;
    }
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

void exactpos_pns(string fn_input, string path1, string fn_output, std::set < pair<int,int> > & regionPos_set){
  ifstream fin;
  fin.open(fn_input.c_str());
  assert(fin);
  map< int, std::set <pair<int, int> > > clipLeft_fn;
  stringstream ss;
  string line;
  int clipLeft, region_begin, region_end;
  getline(fin,line);
  while (getline(fin,line)) {
    ss.clear(); ss.str( line );
    ss >> region_begin >> region_end;
    while ( ss >> clipLeft) 
      clipLeft_fn[clipLeft].insert( make_pair(region_begin, region_end)) ;
  }
  fin.close();
  ofstream fout(fn_output.c_str());
  fout << "clipLeft clipRight pn(count)\n";
  for (map< int, std::set < pair<int, int> > >::iterator pi = clipLeft_fn.begin(); pi != clipLeft_fn.end(); pi++ ) {
    int clipLeft = pi->first;
    map <string, std::set<string> > pn_qName;
    std::set < pair<int, int> > fns;
    size_t fns_size = (pi->second).size() ;
    fns.swap(pi->second);
    neighbor_pos_filename(regionPos_set, fns, clipLeft);
    assert ( fns.size() >= fns_size );
    for ( std::set < pair<int, int> >::iterator fi = fns.begin(); fi != fns.end(); fi++ ) {
      fin.open( (path1 + pos_to_fn(*fi) ).c_str());
      string line, pn, qName;
      while (getline(fin, line)) {
	if (read_match_clipLeft(line, clipLeft, pn , qName))
	  pn_qName[pn].insert(qName);  // some reads are appear more than once 
      }
      fin.close();
    }
    fout << clipLeft << " -1" ;  // clipRight unknown for now
    for (map <string, std::set<string> >::iterator pi = pn_qName.begin(); pi != pn_qName.end(); pi++) 
      fout << " " << pi->first << "," << (pi->second).size() ;
    fout << endl;
  }
  fout.close();
}


bool is_clipread(int & clipLen, seqan::CharString & clipSeq, seqan::BamAlignmentRecord & record, int thisEnd, int clipLeft, int clipRight, int offset_left, int offset_right ) {  
  if ( has_soft_first(record, CLIP_BP) and abs( record.beginPos - clipRight) <=  offset_right 
       and trim_clip_soft_first(record, clipLen, clipSeq, CLIP_Q) )
    return true;
  if ( has_soft_last(record, CLIP_BP) and abs( thisEnd - clipLeft ) <= offset_left 
       and trim_clip_soft_last(record, clipLeft, clipSeq, CLIP_Q) )
    return true;
  return false;
}

string classify_read( seqan::BamAlignmentRecord & record, int clipLeft, int clipRight, int offset_left, int offset_right ,int & clipLen, seqan::CharString & clipSeq) {
  int thisEnd = record.beginPos + (int)getAlignmentLengthInRef(record);   
  if (record.rID != record.rNextId or abs(record.tLen) >= DISCORDANT_LEN) {
    if ( is_clipread(clipLen, clipSeq, record, thisEnd, clipLeft, clipRight, offset_left, offset_right) ) 
      return "clip_read";  // prefer to use it as clip_read, if it's also alu_read
    if ( ( thisEnd < clipLeft + 5 and !hasFlagRC(record) ) or ( record.beginPos > clipRight - 5 and hasFlagRC(record) ) ) 
      return "alu_read";
    return "useless";    
  }  
  bool read_is_left = left_read(record);
  // otherwise only look at proper mapped reads 
  if ( hasFlagRC(record) == hasFlagNextRC(record) or read_is_left == hasFlagRC(record) ) 
    return "useless";  
  if ( is_clipread(clipLen, clipSeq, record, thisEnd, clipLeft, clipRight, offset_left, offset_right) ) 
    return "clip_read";
  int trimb, trime;
  get_trim_length(record, trimb, trime, CLIP_Q);
  if (trimb + trime < (int) length(record) - CLIP_BP and count_non_match(record) <= 5 and
      record.beginPos + trimb <= clipLeft - CLIP_BP and thisEnd - trime >= clipLeft + CLIP_BP ) 
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

void clip_skip_unknow_reads(string chrn, ofstream &ftmp1, vector< pair<int, int> > & insert_pos, BamFileHandler *bam_fh, string file_clipreads,  map < int, vector<ALUREAD_INFO> > & rid_alureads, map<string, int> & rg_to_idx, int alucons_len) {
  ofstream fout(file_clipreads.c_str());    
  fout << "clipPos seq clipLen clipPos_by_cigar cigar rgIdx tLen\n";   // seq is after correction 
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
    map < seqan::CharString, string > qname_iread;
    map < seqan::CharString, string > qname_foutStr;
    while ( true ) {
      read_status = bam_fh->fetch_a_read(chrn, region_begin, region_end, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_insert_read(record)) continue;      
      map < seqan::CharString, string >::iterator qi = qname_iread.find(record.qName);
      if ( qi != qname_iread.end() and qi->second != "unknow_read")
	continue;

      int clipLen;
      seqan::CharString clipSeq;
      string iread = classify_read(record, clipLeft, clipRight, offset_left, offset_right, clipLen, clipSeq);            
      
      if ( iread == "alu_read") {
	ALUREAD_INFO one_aluRead = ALUREAD_INFO(record.qName, clipLeft, record.pNext, length(record.seq), hasFlagRC(record)==hasFlagNextRC(record) );
	rid_alureads[record.rNextId].push_back(one_aluRead);
      } else {
	bool conflict = false;
	string ir = "";   
	get_mapVal(qname_iread, record.qName, ir);  	
	if  ( (iread == "clip_read" and ir == "skip_read") or (iread == "skip_read" and ir == "clip_read") ) 
	  conflict = true;	
	if  (iread == "clip_read" and !conflict) {
	  qname_iread[record.qName] = iread;
	  stringstream ssf;
	  ssf << clipLeft << " " << clipSeq  << " " << clipLen << " ";
	  if (record.cigar[0].operation == 'S') ssf << record.beginPos; 	    
	  else ssf << record.beginPos + (int)getAlignmentLengthInRef(record);	    
	  int rgIdx = get_rgIdx(rg_to_idx, record);
	  ssf << " " << get_cigar(record) << " " << rgIdx << " " << record.tLen << endl;
	  qname_foutStr[record.qName] = ssf.str();
	}		
	if (iread == "skip_read" and !conflict)
	  qname_iread[record.qName] = iread;	
	if ( iread == "unknow_read" or conflict) {    
	  qname_iread[record.qName] = "unknow_read";
	  int rgIdx = get_rgIdx(rg_to_idx, record);
	  stringstream rg_ss;
	  rg_ss << " " << rgIdx << ":" << abs(record.tLen) + alucons_len;
	  rg_str[record.qName] = rg_ss.str();
	}
      }
    }   
    vector<string>  str_unknowreads;   
    int cnt_skipreads = 0;   
    int cnt_midreads = 0;   // midreads (inserting a alu seq) = clip + alu reads     
    for (map < seqan::CharString, string >::iterator qi = qname_iread.begin(); qi != qname_iread.end(); qi++ ) {
      if ( qi->second == "skip_read") {
	cnt_skipreads++;
      } else if ( qi->second == "clip_read") {
	cnt_midreads++;
	fout << qname_foutStr[ qi->first ];
      } else if ( qi->second == "unknow_read") {
	str_unknowreads.push_back( rg_str[qi->first] );	
      }
    }    
    if (str_unknowreads.empty() and !cnt_midreads)
      continue;
    ftmp1 <<  chrn << " " << clipLeft << " " << clipRight << " " << alucons_len  << " " << cnt_midreads << " " <<  cnt_skipreads << " " << str_unknowreads.size() ;
    for  ( vector< string > ::iterator si = str_unknowreads.begin(); si != str_unknowreads.end(); si++ )
      ftmp1 << " " << *si ;
    ftmp1 << endl;
    
  }
}

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
    map <string, string>::iterator qi = qnames.find( toCString( (*ai).qName) );  
    if ( qi != qnames.end() ) mapped_db = qi->second;
    seqan::BamAlignmentRecord record;
    while (bam_fh->fetch_a_read(record) and record.beginPos <= (*ai).pos + 3) {
      if (record.qName == (*ai).qName) {
	string record_seq = toCString(record.seq);
	if ((*ai).sameRC)  seqan::reverseComplement(record_seq); 	
	if ( mapped_db!="na")
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
    map <string, string>::iterator qi = qnames.find( toCString( (*ai).qName) );
    fout << (*ai).clipLeft << " " ;
    if ( qi != qnames.end() ) fout << qi->second ;
    else    fout << alu_type_default ;
    fout << " " << rid << " " << (*ai).pos << " " << (*ai).sameRC << " " << (*ai).readLen << " " <<  (*ai).qName << endl;
  }
}

void count_alumate(string file_input, map <int, int> & insertBegin_cnt ) {
  insertBegin_cnt.clear();
  ifstream fin(file_input.c_str());
  assert(fin);
  stringstream ss;
  string line, tmpv1, tmpv2, alu_type;
  int pos;
  getline(fin, line); // read header
  while (getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> pos >> tmpv1 >> tmpv2 >> alu_type;
    if (alu_type == "na") continue;
    addKey(insertBegin_cnt, pos, 1);
  }
  fin.close();
}

void write_tmp2_chrn( list< IntString> & rows_list, ofstream & fout) {
  rows_list.sort(compare_IntString);
  for (list< IntString>::iterator ri = rows_list.begin(); ri!=rows_list.end(); ri++) 
    fout << (*ri).second << endl;
  rows_list.clear();
}

void empty_alumates_only(map <int, EmpiricalPdf *> & pdf_rg, string chrn, map <int, int > & insertBegin_cnt, list< IntString> & rows_list, string len_alucon) {
  for ( map <int, int >::iterator bi = insertBegin_cnt.begin(); bi != insertBegin_cnt.end(); bi++ ) {
    int pos = bi->first;
    string line0 = chrn + " " + int_to_string(pos)  + " " + int_to_string(pos) + " " + len_alucon + " 0 0 0";
    string output_line;
    if (parseline_del_tmp1(line0, output_line, pdf_rg, bi->second ))
      rows_list.push_back( make_pair(pos, output_line) );
  }	
  insertBegin_cnt.clear();
}

void write_tmp2(string fn_tmp1, string fn_tmp2, map <int, EmpiricalPdf *> & pdf_rg, string fn_alu0){
  ofstream fout(fn_tmp2.c_str());
  fout << "chr insertBegin insertLen midCnt clipCnt unknowCnt 00 01 11\n";  
  ifstream fin(fn_tmp1.c_str());
  assert(fin);
  stringstream ss;
  string line, output_line, chrn, fn_alu, tmpv1, len_alucon;
  string pre_chrn = "";
  int pos;
  map <int, int > insertBegin_cnt;
  //std::set <int> insertBegin_left;
  list< IntString> rows_list;
  getline(fin, line); // read header
  while (getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> chrn >> pos >> tmpv1 >> len_alucon;

    if ( chrn != pre_chrn) {
      if ( pre_chrn != "") {
	empty_alumates_only(pdf_rg, pre_chrn, insertBegin_cnt, rows_list, len_alucon);
	write_tmp2_chrn(rows_list, fout);
      }
      pre_chrn = chrn;      
      fn_alu = replace_str0_str(fn_alu0, chrn, "chr0");
      count_alumate(fn_alu, insertBegin_cnt);
    }
    
    int cnt_alumate = 0;
    get_mapVal( insertBegin_cnt, pos, cnt_alumate ); 
    if (cnt_alumate != 0 )
      insertBegin_cnt.erase( pos );
    
    if (parseline_del_tmp1(line, output_line, pdf_rg, cnt_alumate))
      rows_list.push_back( make_pair(pos, output_line) );
  } 
  empty_alumates_only(pdf_rg, chrn, insertBegin_cnt, rows_list, len_alucon);
  write_tmp2_chrn(rows_list, fout);
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
  string file_pn_used = cf_fh.get_conf( "file_pn_used");
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++ ) {
    check_folder_exists(path0 + *ci +  "/" );
    check_folder_exists(pathClip + *ci +  "/");
    check_folder_exists(pathClip + *ci + "_pos/");
    check_folder_exists(pathCons + *ci);
  }
  
  if (opt == "write_tmps_pn") {     
    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input);
    string file1 = path0 + pn + ".tmp1";
    map<int, string>::iterator rc ;

    //// write pn.tmp1.
    map <string, int> rg_to_idx;    
    parse_reading_group( get_name_rg(cf_fh.get_conf("file_dist_prefix"), pn), rg_to_idx );    
    string header1 = "qname A_chr B_chr A_beginPos B_beginPos A_isRC B_isRC len_read rgIdx\n";
    alu_mate_flag(bam_fh, file1, header1, rg_to_idx);
    cout << "done with reading this bam file, time used " << clocki.elapsed() << endl; // about 2 hours    

    //// write chr*/pn.tmp1  use database to filter alu mate
    string header2 = "qname this_rID this_pos this_isRC len_read alu_rID alu_pos alu_isRC rgIdx alu_type\n";
    MapFO fileMap;    
    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) {
      fileMap[rc->first] = new ofstream( ( path0 + rc->second + "/" + pn + ".tmp1").c_str() );
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << header2 ;
    }     
    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) 
      keep_alu_mate(bam_fh, rc->first,  file1, fileMap, cf_fh.get_conf( "file_alupos_prefix") + rc->second);           
    close_fhs(fileMap, bam_fh->rID_chrn);    
   
    int col_idx = 0;
    for (rc = bam_fh->rID_chrn.begin(); rc != bam_fh->rID_chrn.end(); rc++) {
      string file_tmp1 =  path0 + rc->second +"/" + pn + ".tmp1";
      if (!col_idx) col_idx = get_col_idx(file_tmp1, "this_pos"); 
      sort_file_by_col<int> (file_tmp1, col_idx, true);
    }
    move_files(path0+"tmp1s/", file1);    
   
    alumate_counts_filter(path0, pn, ".tmp1", ".tmp2", ".tmp3", chrns, LEFT_PLUS_RIGHT); 

    //// filter out rep regions,  combine rep regions(dist < REP_MASK_JOIN bp)
    const int REP_MASK_JOIN = 200;
    RepMaskPos *repmaskPos;
    repmaskPos = new RepMaskPos(cf_fh.get_conf( "file_repeatMask"), chrns, REP_MASK_JOIN); 
    filter_location_rep( path0, pn, ".tmp3", ".tmp3st", chrns, repmaskPos);    
    delete bam_fh;
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
    
  } else if (opt == "clipReads_by_pn") { // clip reads for each pn
    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];    
    std::set <string> pns_used;
    read_file_pn_used( file_pn_used, pns_used); // some pn is ignored, due to too many reads
    if ( pns_used.find(pn) == pns_used.end() ){
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bam_input+".bai");
    string file_fa = cf_fh.get_conf("file_fa_prefix");    
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {      
      string fin_pos = path1 + "insert_pos." + *ci; // eg: 10181 positions for chr1      
      string fout_reads = pathClip + *ci + "/" + pn;
      FastaFileHandler * fasta_fh = new FastaFileHandler(file_fa + *ci + ".fa", *ci);
      clipreads_at_insertPos(pn, *ci, bam_fh, fasta_fh, fin_pos, fout_reads);
      delete fasta_fh;
    }
    delete bam_fh;
    cout << "output to " << pathClip << "chr*/" << pn  << endl;

  } else if ( opt == "clipReads_pos_pns" ) { 
    std::set <string> pns_used;
    read_file_pn_used( file_pn_used, pns_used); // some pn is ignored, due to too many reads
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      string tmp_path0 = pathClip + *ci + "/";
      string tmp_path1 = pathClip + *ci + "_pos/";
      combine_clipreads_by_pos(pns_used, tmp_path0, tmp_path1);
      string file_clipRegion = pathClip + *ci + ".clip_region";
      string file_clipPos = pathClip + *ci + ".clip_pn";
      std::set < pair<int, int> > regionPos_set;
      assert( nonempty_files_sorted( tmp_path1, regionPos_set) );
      regions_pos_vote(tmp_path1, file_clipRegion, regionPos_set);
      exactpos_pns(file_clipRegion, tmp_path1, file_clipPos, regionPos_set);
    }

  } else if ( opt == "fixed_delete0_pn" ) {   // write fasta reads for consensus, use clip reads first 

    assert (argc == 4);
    string pn = ID_pn[seqan::lexicalCast<int> (argv[3])];
    std::set <string> pns_used;
    read_file_pn_used(file_pn_used, pns_used); // some pn is ignored, due to too many reads
    if ( pns_used.find(pn) == pns_used.end() ){ // a bit slow
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bam_input+".bai");
    
    string bam_rid_chrn = cf_fh.get_conf("bam_rid_chrn");
    if ( check_file_size(bam_rid_chrn) <= 0 ) {
      ofstream fmap(bam_rid_chrn.c_str());
      for ( map<int, string>::iterator bi = bam_fh->rID_chrn.begin(); bi !=  bam_fh->rID_chrn.end(); bi++) 
	fmap << bi->first << " " << bi->second << endl;
      fmap.close();
    }
     
    AluconsHandler *alucons_fh = new AluconsHandler(cf_fh.get_conf("file_alu_cons"), cf_fh.get_conf("type_alu_cons"));    
    map<string, int> rg_to_idx;
    string file_dist_prefix = cf_fh.get_conf("file_dist_prefix");
    parse_reading_group( get_name_rg(file_dist_prefix, pn), rg_to_idx );
    
    int min_pn = 2;
    if (min_pn != 1)
      cout << "NB: private insertions are ignored!\n";

    string file_tmp0 = pathDel0 + pn + ".tmp0";
    ofstream fout0 ( file_tmp0.c_str());
    fout0 << "chr insertBegin insertEnd estimatedAluLen midClipCnt clipCnt unknowCnt unknowStr\n";  // midClipCnt replaced midCnt
    for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
      vector< pair<int, int> > insert_pos;      
      string file_clipPos = pathClip + *ci + ".clip_pn";
      read_first2col( file_clipPos , insert_pos, true, min_pn);    
      string file_clip2 = pathCons + *ci + "/" + pn + ".clip";
      map < int, vector<ALUREAD_INFO>  > rid_alureads;       // info of alu reads 
      clip_skip_unknow_reads(*ci, fout0, insert_pos, bam_fh, file_clip2, rid_alureads, rg_to_idx, alucons_fh->seq_len);            

      string file_alu2 = pathCons + *ci + "/" + pn + ".aludb";
      ofstream falu ( file_alu2.c_str());  // use database, faster, also discard half of the reads
      falu << "clipPos alu_type rid beginPos sameRC read_len qName\n";
      string alu_type_default = "na";
      for ( map < int, vector<ALUREAD_INFO> >::iterator ri = rid_alureads.begin(); ri != rid_alureads.end(); ri++) {
	string chrn, file_db;
	if ( bam_fh->get_chrn(ri->first, chrn) ) file_db =  path0 + chrn +"/" + pn + ".tmp1";	    
	else file_db = "";
	/* too slow: alumate_reads1(bam_fh, file_db, falu, ri->first, ri->second); */
	alumate_reads2(file_db, falu, ri->first,  ri->second, alu_type_default);
      }
      falu.close();
      sort_file_by_col<int> (file_alu2, 1, true);  // sort by clip pos      
    }   
    fout0.close();

    map <int, EmpiricalPdf *> pdf_rg;    
    string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
    read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
    string file_tmp2 = pathDel0 + pn + ".tmp2";
    string file_alu2 = pathCons + "chr0/" + pn + ".aludb"; 
    write_tmp2(file_tmp0, file_tmp2, pdf_rg, file_alu2);	    
    EmpiricalPdf::delete_map(pdf_rg);
    move_files(pathDel0 + "tmp0s/", file_tmp0);
    move_files(pathDel0 + "tmp2s/", file_tmp2);

  } else if (opt == "fixed_vcf_pns") {
    vector <string> pns;
    read_file_pn_used(file_pn_used, pns); 
    string path_input = pathDel0 + "tmp2s/";
    string tmp_file_pn = path_input + *(pns.begin()) + ".tmp2";
    int col_idx =  get_col_idx(tmp_file_pn, "00");
    assert (col_idx == 7 );
    string fn_pos = pathDel0 + int_to_string( pns.size()) + ".pos";
    filter_by_llh_noPrivate(path_input, ".tmp2", fn_pos, pns, chrns, col_idx);
    string fn_vcf = pathDel0 + int_to_string( pns.size()) + ".vcf";  
    combine_pns_vcf_noPrivate(path_input, ".tmp2", fn_vcf, pns, chrns, col_idx);  

    
  } else if (opt == "debug1") { 

    string pn = ID_pn[0];
    /*
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bai_input);
    bam_fh->jump_to_region("chr1", 120000, 120010);    
    int ni = 0;
    seqan::BamAlignmentRecord record;
    while ( ni ++ < 10 ) {
      bam_fh -> fetch_a_read(record);
    }
    delete bam_fh;
    */

   map <int, EmpiricalPdf *> pdf_rg;    
   string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
   string file_dist_prefix = cf_fh.get_conf("file_dist_prefix");
   read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
   string line0, output_line;

   line0 = "chr1 3 3 299 0 0 4  0:830  0:802  0:831  3:1072";
   if (parseline_del_tmp1(line0, output_line, pdf_rg, 0 ))
     cout << output_line << endl;

   line0 = "chr1 3 3 299 1 0 4  0:830  0:802  0:831  3:1072";
   if (parseline_del_tmp1(line0, output_line, pdf_rg, 0 ))
     cout << output_line << endl;

   line0 = "chr1 3 3 386 1 0 4  0:830  0:802  0:831  3:1072";
   if (parseline_del_tmp1(line0, output_line, pdf_rg, 0 ))
     cout << output_line << endl;
   
//   line0 = "chr9 99 99 299 3 0 1 0:507 0:495";
//   if (parseline_del_tmp1(line0, output_line, pdf_rg, 0 ))
//     cout << output_line << endl;


  } else {

    cout << "unknown option ! \n";

  }

  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
