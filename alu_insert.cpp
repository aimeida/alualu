#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

inline string get_name(string path, string fn, string suffix){ return path + fn + suffix;}

void read_line(ifstream &fin, list< READ_INFO *> &lr_reads, bool &readable){
  string line, type_flag, qname, alu_type;
  int this_chr, this_pos, bad_chr, bad_pos, len_read;  
  stringstream ss;  
  readable = false;
  if (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> len_read >> alu_type;      
    lr_reads.push_back(new READ_INFO(this_pos, len_read, alu_type));
    readable = true;
  } 
}

void read_line(ifstream &fin, list< READ_INFO *> &lr_reads, bool &readable, list<int> &pos_to_scan){
  string line, type_flag, qname, alu_type;
  int this_chr, this_pos, bad_chr, bad_pos,  len_read;  
  stringstream ss;  
  readable = false;
  if (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> len_read >> alu_type;      
    lr_reads.push_back(new READ_INFO(this_pos, len_read, alu_type));
    if ( (!pos_to_scan.size()) or this_pos > pos_to_scan.back() ) 
      pos_to_scan.push_back( this_pos );
    readable = true;
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

void check_this_pos(list <READ_INFO *> &lr_reads, int this_pos, int &lr_num, int &rr_num, map <string, int> &alu_type_count){
  alu_type_count.clear();
  int offset = 5;
  int this_end = this_pos + DEFAULT_READ_LEN;
  list <READ_INFO *>::iterator ri;
  ri = lr_reads.begin();
  while (ri != lr_reads.end() ) {
    if ( (*ri)->beginPos < this_pos - SCAN_WIN_LEN ) {
      delete *ri;
      ri = lr_reads.erase(ri); // <==> lr_reads.erase(ri++)
      continue;
    } 
    if ( (*ri)->endPos > this_end + SCAN_WIN_LEN ) break;
    if ( (*ri)->endPos < this_pos + offset) {
      lr_num++;
      addKey(alu_type_count, (*ri)->alu_type );
    } else if ( (*ri)->beginPos > this_end - offset ){
      rr_num++;
      addKey(alu_type_count, (*ri)->alu_type );
    }
    ri++;
  }
}

void alumate_counts_filter(string fn, vector<string>  &chrns){
  bool readable;
  ofstream fout(fn.c_str());
  int this_pos, next_pos, lr_num, rr_num;
  map <string, int> alu_type_count;
  string line;
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    list< READ_INFO *> lr_reads;
    list <int> pos_to_scan;    
    list <int>::iterator pts;
    ifstream fin( (fn + "." +  *ci).c_str());
    assert(fin); // init, 700 * 30 /100 = 210 
    getline(fin, line); // skip header
    for (int i = 0; i < 10; i++) read_line(fin, lr_reads, readable); 
    if (!readable) break;
    for (int i = 0; i < 10; i++) read_line(fin, lr_reads, readable, pos_to_scan); 
    if (!readable) break;
    //cerr << *ci << " pos_to_scan " << pos_to_scan.size() << endl;
    while ( pos_to_scan.size() ) {
      pts = pos_to_scan.begin();
      this_pos = *pts;
      next_pos = *(++pts); 
      pos_to_scan.pop_front();
      lr_num = 0;
      rr_num = 0;      
      check_this_pos(lr_reads, this_pos, lr_num, rr_num, alu_type_count);
      if (lr_num + rr_num >= LEFT_PLUS_RIGHT) 
	fout << *ci << " " << this_pos  << " " << lr_num  << " " << rr_num << " " <<  major_type(alu_type_count) <<  endl;      
      while (pos_to_scan.size() < 5 or (lr_reads.back())-> beginPos < next_pos + SCAN_WIN_LEN ) {
	read_line(fin, lr_reads, readable, pos_to_scan); 
	if (!readable) break;
      }	
    }
    fin.close();
  }
  fout.close();
}

void join_location(string file1, string file2, int pos_dif, int count_dif){ // AND
  string chrn, chrn_pre, alu_type, maxC_alu_type;
  int maxCount, nowCount, pos, lr_num, rr_num, pos_pre, lr_num_pre, rr_num_pre;
  vector<int>  maxCount_pos;
  ifstream fin(file1.c_str());
  assert(fin);
  ofstream fout(file2.c_str());
  //int ti = 0;
  while ( fin >> chrn >> pos >> lr_num >> rr_num >> alu_type) {
    maxC_alu_type = alu_type;
    if (chrn_pre != chrn) {
      chrn_pre = chrn;
      pos_pre = pos;
      lr_num_pre = lr_num;
      rr_num_pre = rr_num;
      maxCount_pos.clear();
      maxCount_pos.push_back(pos);
      maxCount = lr_num + rr_num;
      continue;
    }     
    if ( pos - pos_pre <= pos_dif and abs(lr_num - lr_num_pre) <= count_dif and abs(rr_num - rr_num_pre) <= count_dif) { // same block
      nowCount = lr_num + rr_num;
      if ( nowCount > maxCount) {
	maxCount = nowCount;
	maxCount_pos.clear();
	maxCount_pos.push_back(pos);
	maxC_alu_type = alu_type;  // a bit random 
      } else if ( nowCount == maxCount) {
	maxCount_pos.push_back(pos);
      }
    } else { // create new block      
      fout << chrn << " " << maxCount << " " << maxCount_pos.front() << " " << maxCount_pos.back() << " " << maxC_alu_type <<  endl;
      maxCount = lr_num + rr_num;
      maxCount_pos.clear();
      maxCount_pos.push_back(pos);      
    }
    pos_pre = pos;
    lr_num_pre = lr_num;
    rr_num_pre = rr_num;
  }
  fin.close();
  fout.close();
}

void filter_location_rep(string file1, string file2, RepMaskPos *repmaskPos){
  ifstream fin(file1.c_str());
  assert(fin);
  ofstream fout(file2.c_str());
  string line, chrn, chrn_pre = "chr0";
  stringstream ss;
  int r_num, pos_left, pos_right;
  vector<int>::iterator bi, ei;
  int be_size, bei;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> chrn >> r_num >> pos_left >> pos_right;
    if (chrn != chrn_pre) {
      bi = repmaskPos->beginP[chrn].begin();
      ei = repmaskPos->endP[chrn].begin();
      bei = 0;
      be_size = repmaskPos->beginP[chrn].size();
    } 
    while ( pos_left >= (*ei) and bei < be_size) { bi++; ei++; bei++; }
    if ( min(*ei, pos_right) - max(*bi, pos_left) <= 0)  fout << line << endl;  // not in Alu
  }  
  fin.close();
  fout.close();
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
  if (!header.empty()) fout << header << " sub_type " <<  endl;
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
  string line, type_flag, qname, sub_type;
  int this_chr, this_pos, bad_chr, bad_pos, len_read;  
  stringstream ss;  
  ifstream fin( fn.c_str() );
  if (has_header) getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> bad_chr >> bad_pos >> this_chr >> this_pos >> len_read >> sub_type;      
    *(fileMap[this_chr]) << type_flag << " " << qname <<" " << this_chr <<" " << this_pos <<" " << bad_chr <<" " << bad_pos << " " << len_read << " " << sub_type << endl;
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

    cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    get_rID_chrn(bam_input, chrns, rID_chrn);    
    MapFO fileMap;
    string file1_prefix = get_name(path0, pn, ".tmp1");
    string file2_prefix = get_name(path0, pn, ".tmp2");
    string file3_prefix = get_name(path0, pn, ".tmp3");    
    string header = "flag qname bad_chr bad_pos this_chr this_pos len_read"; // bad_*: potential, need check
    
    // step 1, about 2h
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
    // move tmp files to sub_folder
    path_move = path0+"tmp1s/";
    check_folder_exists(path_move);
    system(("mv "+path0+"*tmp1.chr* " + path_move).c_str());
    cerr << "filter alu mate, done\n";        
    /* write pn.tmp3.chrX
       rewrite for counting, switch 2 column and sort by pos */
    header = "flag qname this_chr this_pos bad_chr bad_pos len_read sub_type"; 
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
    system(("mv "+path0+"*tmp2.chr* " + path_move).c_str());    

    /* sort  pn.tmp3.chrX */
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

  } else if (opt == 12) { // combine positions from multiple individuals, can be added to opt = 1 in the end 
    //  scan for potential regions
    seqan::lexicalCast2(idx_pn, argv[3]);
    string pn = ID_pn[idx_pn];
    string file3_prefix = get_name(path0, pn, ".tmp3");    
    cout << "reading " << file3_prefix << endl;
    alumate_counts_filter(file3_prefix, chrns);
    //    path_move = path0+"tmp3s/";
//    check_folder_exists(path_move);
//    system(("mv "+path0+"*tmp3.chr* " + path_move).c_str());        
    exit(0);
    string file4 = get_name(path0, pn, ".tmp4" );
    join_location(file3_prefix, file4, 50, 4);
    

    // filter out rep regions,  combine rep regions(dist < 100bp)
    RepMaskPos *repmaskPos;
    repmaskPos = new RepMaskPos(read_config(config_file, "file_repeatMask"), 100); 
    string file5 = get_name(path0, pn, ".tmp5" );          
    filter_location_rep(file4, file5, repmaskPos);
    
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
