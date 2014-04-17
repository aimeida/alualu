#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

inline string get_name(string path, string fn, string suffix){ return path + fn + suffix;}

void parse_line2(string &line, int &this_pos, int &len_read, string &alu_type) {
  string type_flag, qname;
  int this_chr, bad_chr, bad_pos;
  stringstream ss;
  ss.clear(); ss.str( line );
  ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> len_read >> alu_type;
  assert(len_read > 80); // in case of reading error 
}

int parse_line3( string &line, int & pos) {
  int lr_num, rr_num;
  stringstream ss;
  ss.str(line);
  ss >> pos >> lr_num >> rr_num;  
  return lr_num + rr_num;
}


void get_scan_pos(string fn, list <int> & pos_to_scan, int combine_bp/*=8*/){ 
  string line, type_flag, qname, alu_type;
  int this_chr, this_pos, bad_chr, bad_pos,  len_read, previous_end;  
  stringstream ss;
  ifstream fin( fn.c_str());
  assert(fin); // init, 700 * 30 /100 = 210 
  for (int i = 0; i < 4; i++) getline(fin, line); 
  ss.str(line);
  ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> len_read;
  pos_to_scan.push_back(this_pos - combine_bp);
  while (getline(fin, line)) {
    previous_end = this_pos + len_read;
    parse_line2(line, this_pos, len_read, alu_type);
    /////if (previous_end + combine_bp == 1574154 or this_pos - combine_bp == 1574154)
    if (this_pos - previous_end > 3 * combine_bp) pos_to_scan.push_back(previous_end + combine_bp);
    if (this_pos - pos_to_scan.back() > 2 * combine_bp) pos_to_scan.push_back(this_pos - combine_bp);
  }
}                                     

string major_type( map <string, int> &alu_type_count) {
  map <string, int>::iterator atc ;
  int max_count = 0;
  for (atc = alu_type_count.begin(); atc != alu_type_count.end(); atc++) 
    if (atc->second > max_count) max_count = atc->second;
  for (atc = alu_type_count.begin(); atc != alu_type_count.end(); atc++) 
    if (atc->second == max_count) return atc->first;
  return "err";
}

// need smarter way to combine positions ==> minimum number of regions to cover all reads info
void join_location(string file1, string file2, int pos_dif, int max_region_len){ 
  ifstream fin( file1.c_str());
  assert(fin);    
  string line;
  getline(fin, line);
  int pos, pos_pre;
  int nowCount = parse_line3(line, pos_pre);
  //cout << "nowC " << nowCount << endl;
  int maxCount = nowCount;
  vector<int>  maxCount_pos;
  maxCount_pos.push_back(pos_pre);
  stringstream sout; // output small, save in memory
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

void alumate_counts_filter(string fn_input, string fn_output, vector<string>  &chrns){
  int lr_num, rr_num, this_pos, len_read;
  string line, alu_type;
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    string chrn = *ci;
    string f_input = fn_input + "." + chrn;
    string f_output = fn_output + "." + chrn;
    list <int> pos_to_scan;    
    get_scan_pos(f_input, pos_to_scan, 8);
    ////for (list <int>::iterator pts = pos_to_scan.begin(); pts != pos_to_scan.end(); pts++) cout << *pts << endl;    
    list< READ_INFO *> lr_reads;    
    ifstream fin ( f_input.c_str());
    assert(fin); // init, 700 * 30 /100 = 210 
    ofstream fout( f_output.c_str());
    getline(fin, line); // skip header
    getline(fin, line); 
    parse_line2(line, this_pos, len_read, alu_type);
    lr_reads.push_back(new READ_INFO(this_pos, 100, "AluY")); 
    bool readable = true;
    while ( readable and pos_to_scan.size()) {
      int check_pos = pos_to_scan.front();
      pos_to_scan.pop_front();
      while (readable and (lr_reads.empty() or lr_reads.back()-> endPos < check_pos + SCAN_WIN_LEN) ) {
	if (getline(fin, line)) {
	  parse_line2(line, this_pos, len_read, alu_type);
	  lr_reads.push_back(new READ_INFO(this_pos, len_read, alu_type));
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
	if ( (*ri)->endPos < check_pos) lr_num++;
	else if ( (*ri)->beginPos > check_pos) rr_num++;
	ri++;
	if ( (*ri)->endPos > check_pos + SCAN_WIN_LEN ) break;
      }
      if (lr_num + rr_num >= LEFT_PLUS_RIGHT) fout << check_pos << " " << lr_num  << " " << rr_num << endl;	
    }
    fin.close();
    fout.close();
    join_location(f_output, f_output, MAX_POS_DIF, MAX_LEN_REGION);    		       
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
      while ( pos_left >= (*ei) and bei < be_size) { bi++; ei++; bei++; }
      if ( min(*ei, pos_right) - max(*bi, pos_left) <= 0)  fout << line << endl;  // not in Alu
    }  
    fin.close();
    fout.close();
  }
}

int read_sort_by_col(string fn, int coln, bool has_header, set< RowInfo, compare_row > &rows) {
  string line, tmpv;
  int pos;
  stringstream ss;
  ifstream fin( fn.c_str());
  assert(fin); 
  rows.clear();
  if (has_header) getline(fin, line);
  size_t rown = 0;
  while (getline(fin, line)) {
    rown++;
    ss.clear(); ss.str( line );
    for (int i = 0; i < coln-1; i++) ss >> tmpv;
    ss >> pos;
    rows.insert( make_pair(pos, line) );
  }
  fin.close();
  if (rows.size() != rown ) cerr << "##### ERROR #### " << fn << endl;
  return rows.size();
}


int read_sort_by_col(string fn, int coln, bool has_header, list< RowInfo> &rows_list) {
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
  rows_list.sort(compare_list);
  if (rows_list.size() != rown ) cerr << "##### ERROR #### " << fn << endl;
  return rows_list.size();
}

void keep_alu_mate(string file1, string file2, AluRefPos * alurefpos, string &header){
  set< RowInfo, compare_row > rows;
  read_sort_by_col(file1, 4, !header.empty(), rows);      
  // filter reads and keep only alu_mate
  ofstream fout(file2.c_str());
  if (!header.empty()) fout << header << " alu_type " <<  endl;
  vector<int>::iterator bi = alurefpos->beginV.begin();
  vector<int>::iterator ei = alurefpos->endV.begin();
  vector<string>::iterator ti = alurefpos->typeV.begin();
  int bei = 0;
  int be_size = alurefpos->beginV.size() - 1; // -1 ! otherwise seg fault
  int left_pos, len_overlap;
  int ni = 0;
  for (set< RowInfo, compare_row >::iterator ri = rows.begin(); ri!=rows.end(); ri++, ni++) {
    left_pos = (*ri).first;   
    while ( ( left_pos >= (*ei)) and bei < be_size) { bi++; ei++; ti++; bei++; }
    //cout << ni << " " << rows.size() << " " << (*ri).second << " " << *bi << " " << *ei << " " << *ti << " " <<  bei << " " << be_size <<  endl;
    if ( (len_overlap = min(*ei, left_pos + DEFAULT_READ_LEN) - max(*bi, left_pos)) > ALU_MIN_LEN ) 
      fout << (*ri).second << " " << *ti << endl;    
  }
  fout.close();
}

void reorder_column(string fn, MapFO &fileMap, bool has_header){
  string line, type_flag, qname, alu_type;
  int this_chr, this_pos, bad_chr, bad_pos, len_read;  
  stringstream ss;  
  ifstream fin( fn.c_str() );
  if (has_header) getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> bad_chr >> bad_pos >> this_chr >> this_pos >> len_read >> alu_type;      
    *(fileMap[this_chr]) << type_flag << " " << qname <<" " << this_chr <<" " << this_pos <<" " << bad_chr <<" " << bad_pos << " " << len_read << " " << alu_type << endl;
  }
  fin.close();
}

void write_all_location( vector<string> &fns, vector <int>&idx_pns, vector<string> &chrns, string file_prefix) {
  MapSFO fileMap;
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) 
    fileMap[*ci] = new ofstream( ( file_prefix + *ci ).c_str() );
  ifstream fin;
  stringstream ss;
  string line, chrn;
  vector <int>::iterator ii = idx_pns.begin();
  for (vector<string>::iterator fi = fns.begin(); fi!= fns.end(); fi++, ii++){
    fin.open( (*fi).c_str() );
    while (getline(fin, line)) {
      ss.clear(); ss.str( line );
      ss >> chrn;
      *(fileMap[chrn]) << line << " " << *ii << endl;
    }
    fin.close();
  }
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
    delete fileMap[*ci];
    list< RowInfo> rows_list;
    string fn = file_prefix + *ci;
    read_sort_by_col(fn, 3, false, rows_list);      
    ofstream fout(fn.c_str());
    for (list< RowInfo>::iterator ri = rows_list.begin(); ri!=rows_list.end(); ri++) 
      fout << (*ri).second << endl;
    fout.close();
    rows_list.clear();
  }
}

void join_location2(vector<string> &chrns, string prefix_if, string prefix_of, int pos_dif, size_t min_num_pn){
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    ifstream fin( (prefix_if+*ci).c_str() );
    assert(fin);
    ofstream fout( (prefix_of+*ci).c_str() );
    int idx, pa, pb, num;
    string chrn, alu_type;
    set<int> ids;
    fin >> chrn >> num >> pa >> pb >> alu_type >> idx;
    int pa_block = pa;
    int pb_block = pb;
    int pa_pre = pa;
    int reads_num = num;
    ids.insert(idx);
    while (fin >> chrn >> num >> pa >> pb >> alu_type >> idx) {
      if ( pa <= pb_block + pos_dif ) { // same block
	ids.insert(idx);
	pb_block = max(pb_block, pb);
	reads_num += num;
      } else { // new block
	if ( ids.size() >= min_num_pn) {
	  fout << pa_block << " " << pb_block << " " << ids.size() << " " << reads_num << " " << alu_type;
	  for (set<int>::iterator si = ids.begin(); si != ids.end(); si++ )
	    fout << " " << *si;
	  fout << endl;
	}
	ids.clear();
	pa_block = pa;
	pb_block = pb;
	reads_num = num;
      }
      pa_pre = pa;
    }
    if ( ids.size() >= min_num_pn) // last block
      fout << pa_block << " " << pb_block << " " << ids.size() << " " << reads_num << endl;
    fin.close();
    fout.close();
  }
}

void write_insert_fasta(string bam_input, string fin_pos, MapFO & fileMap, map<int, seqan::CharString> const &rID_chrn){
  
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  int opt, idx_pn;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];

  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  string path0 = read_config(config_file, "file_alu_insert0") ;    
  string path1 = read_config(config_file, "file_alu_insert1") ;    
  string path_move;
  check_folder_exists(path0);
  check_folder_exists(path1);

  boost::timer clocki;    
  clocki.restart();

  map<int, string> ID_pn;
  get_pn(read_config(config_file, "file_pn"), ID_pn);
  map<int, seqan::CharString> rID_chrn;
  
  if (opt == 1) {     // takes about 2.5 hrs 
    chrns.push_back("chrX");
    chrns.push_back("chrY");
    seqan::lexicalCast2(idx_pn, argv[3]);
    string pn = ID_pn[idx_pn];
    cout << "reading pn: " << idx_pn << " " << pn << "..................\n";
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    get_rID_chrn(bam_input, chrns, rID_chrn);    
    MapFO fileMap;
    string file1_prefix = get_name(path0, pn, ".tmp1");
    string file2_prefix = get_name(path0, pn, ".tmp2");
    string file3_prefix = get_name(path0, pn, ".tmp3");    

    // step 1, about 2h
    string header = "flag qname bad_chr bad_pos this_chr this_pos len_read"; // bad_*: potential, need check
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      string chrn = toCString(rc->second);
      fileMap[rc->first] = new ofstream( (file1_prefix + "." + chrn).c_str() );
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << header << endl;
    }
    alu_mate_flag(bam_input, rID_chrn, fileMap);
    cerr << "read records, done\n";    

    // step 2, if mate mapped to alu. 5 min. 
    string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      delete fileMap[rc->first];
      string chrn = toCString(rc->second);
      AluRefPos *alurefpos = new AluRefPos(file_alupos_prefix + chrn);          
      keep_alu_mate(file1_prefix + "." + chrn, file2_prefix + "." + chrn, alurefpos, header);       
      delete alurefpos;      
      } 
      
    path_move = path0+"tmp1s/";
    check_folder_exists(path_move);
    system(("mv " + path0 + pn + ".tmp1.chr* " + path_move).c_str());
    cerr << "filter alu mate, done\n";        

    // write pn.tmp3.chr*, rewrite for counting, switch 2 column and sort by pos
    header = "flag qname this_chr this_pos bad_chr bad_pos len_read alu_type"; 
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      string file3 = file3_prefix + "." + toCString(rc->second);
      fileMap[rc->first] = new ofstream( file3.c_str() );
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << header << endl ;      
    }    
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) 
      reorder_column(file2_prefix + "." + toCString(rc->second), fileMap, true);
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) 
      delete fileMap[rc->first];    

    path_move = path0+"tmp2s/";
    check_folder_exists(path_move);
    system(("mv " + path0 + pn + ".tmp2.chr* " + path_move).c_str());    

    /// sort  pn.tmp3.chrX
    cerr << "sorting starts \n";
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      string file3 = file3_prefix + "." + toCString(rc->second);
      set< RowInfo, compare_row > rows;
      read_sort_by_col(file3, 4, true, rows);      
      ofstream fout(file3.c_str());
      fout << header << endl ;      
      for (set< RowInfo, compare_row >::iterator ri = rows.begin(); ri!=rows.end(); ri++) 
	fout << (*ri).second << endl;
      fout.close();
    }    
    
#ifdef DEBUG_MODE   // only look at chr1 for now 
    chrns.clear();
    chrns.push_back("chr1");
#endif       
   
    string file4_prefix = get_name(path0, pn, ".tmp4");    
    alumate_counts_filter(file3_prefix, file4_prefix, chrns);  //  20 mins
    path_move = path0+"tmp3s/";
    check_folder_exists(path_move);
    system(("mv "+ path0 + pn + ".tmp3.chr* " + path_move).c_str());        
    
    // filter out rep regions,  combine rep regions(dist < 100bp)
    RepMaskPos *repmaskPos;
    repmaskPos = new RepMaskPos(read_config(config_file, "file_repeatMask"), 100); 
    string file4st = get_name(path0, pn, ".tmp4st");  // subset of *tmp4
    filter_location_rep(file4_prefix, file4st, chrns, repmaskPos);
    

  } else if (opt == 2) { // get alu sequence for building consensus sequence of the inserted region 
    chrns.push_back("chrX");
    chrns.push_back("chrY");
    seqan::lexicalCast2(idx_pn, argv[3]);
    string pn = ID_pn[idx_pn];
    cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
    string path_aluseq = read_config(config_file, "file_alu_seq");
    check_folder_exists(path_aluseq);
    path_move = path0 + "tmp2s/";
    string fin_pos = get_name(path_move, pn, ".tmp2");
    string fout_fa = get_name(path_aluseq, pn, ".tmp1");
    MapFO fileMap;
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    get_rID_chrn(bam_input, chrns, rID_chrn);
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      string chrn = toCString(rc->second);
      fileMap[rc->first] = new ofstream( (fout_fa + "." + chrn).c_str() );
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << "pos hasRC pair_hasRC pair_seq\n";
    }
    write_insert_fasta(bam_input, fin_pos, fileMap, rID_chrn);
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) 
      delete fileMap[rc->first];

  }  else if (opt == 3) { // combine positions from multiple individuals
    chrns.push_back("chrX");
    chrns.push_back("chrY");
    vector<int> idx_pns;
    vector<string> fns;
    ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
    while (fin >> idx_pn) {
      idx_pns.push_back(idx_pn);
      fns.push_back(get_name(path0, ID_pn[idx_pn], ".tmp5" ));      
    }
    write_all_location(fns, idx_pns, chrns, path1+"tmp.insert_pos."); 
    join_location2(chrns, path1+"tmp.insert_pos.", path1+"insert_pos.", 40, 3); //at least 3 pn has insertion       
  }
  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
