#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

typedef map<int, ofstream* > MapFO;
typedef map<string, ofstream* > MapSFO;
typedef pair<int, string > RowInfo;

#define ALU_MIN_LEN 50
#define DEFAULT_READ_LEN 100 // if unknown, use default
#define FLAG_LONG -1
#define FLAG_RC -2
#define DISCORDANT_LEN 1000 
#define SCAN_WIN_LEN 400  // 0.99 quantile = 400
#define LEFT_PLUS_RIGHT 5 // minimum sum of left and right reads

inline string get_name(string path, string fn, string suffix){ return path + fn + suffix;}

//write down alu_mate and flag, 
bool alu_mate_flag( string bam_input, map<int, seqan::CharString> &rID_chrn, MapFO &fileMap){
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );

  int rID_pre = -1, ia = 0;
  map<seqan::CharString, pair<int,int> > same_chr_left; // value: (pos, count)
  map<seqan::CharString, pair<int,int> >::iterator sc;

  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( hasFlagQCNoPass(record) or hasFlagDuplicate(record) or hasFlagUnmapped(record) or hasFlagNextUnmapped(record) or (not hasFlagMultiple(record))) continue;
    //// this FLAG is broken ==> if (hasFlagSecondary(record)) writeRecord(bamStreamOut, record);
    if ( rID_chrn.find(record.rID) == rID_chrn.end() ) continue;
    if ( record.rID != rID_pre) {
      if (same_chr_left.size() ) {  // do not print out left reads 
	for (map<seqan::CharString, pair<int,int> >::iterator mc = same_chr_left.begin(); mc!= same_chr_left.end(); mc++)
	  cerr << "left " << mc->first << " " << rID_pre << " " << (mc->second).first << " -1 -1 " << (mc->second).second << endl;
	same_chr_left.clear();
      }
      cerr << "done with " << rID_chrn[rID_pre] << " " << rID_pre << endl;
      rID_pre = record.rID;	
    }       
    if (ia++ % 1000000 == 0) cerr << ia << " " << same_chr_left.size() << endl;
    if (record.rID != record.rNextId) { // NB: redundant! one pair read twice, thus are written twice as well
      if (rID_chrn.find(record.rNextId) != rID_chrn.end())  // both chrn exist
	*(fileMap[record.rNextId]) << "dif_chr " << record.qName << " " << record.rNextId << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << numOfBestHits(record) <<  endl;      
      continue;
    }
    int hits = numOfBestHits(record);      
    if (record.beginPos < record.pNext) {	// left read
      if ( hits > 1 and not_all_match(record) ) { 
	same_chr_left[record.qName] = make_pair(record.beginPos, hits);
      } else if ( (abs(record.tLen) > DISCORDANT_LEN ) or (hasFlagNextRC(record) == hasFlagRC(record))) {  // here to define iLong_RC
	same_chr_left[record.qName] = make_pair(record.beginPos, FLAG_LONG); // non-redundant
      }	
    } else if (record.beginPos > record.pNext) { // right read
      sc = same_chr_left.find(record.qName);  
      if (sc == same_chr_left.end()) {   // left read is unique
	if (hits > 1 and not_all_match(record) )
	  *(fileMap[record.rID]) << "mul_map " << record.qName << " " << record.rID <<  " " << record.beginPos << " " << record.rID << " " << record.pNext << " " << hits << endl;
      } else {                   // left read not multi_mapped
	if ( (sc->second).second == FLAG_LONG) {  // NB, redundant! written twice when right reads are read
	  *(fileMap[record.rID]) << "ilong_RC " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " 1" << hits << endl;
	  *(fileMap[record.rID]) << "ilong_RC " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rID << " " << record.pNext << " " << hits << endl;
	} else if (hits == 1) {  // at least one pair is ok
	  *(fileMap[record.rID]) << "mul_map " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << (sc->second).second << endl;
	}
	same_chr_left.erase(sc); // rm left read
      }	      
    }    
  } 
  seqan::close(inStream);     
  return 0;
}

void left_or_right(list< pair<int, int> > &left_reads, list< pair<int, int> > &right_reads, string &type_flag, int pos, int hits, bool mapped_left) {
  if (type_flag == "dif_chr") {
    left_reads.push_back( make_pair(pos, hits) );
    right_reads.push_back( make_pair(pos, hits) );
  } else if ( type_flag == "mul_map" or type_flag == "ilong_RC") {
    if (mapped_left) left_reads.push_back( make_pair(pos, hits));
    else right_reads.push_back( make_pair(pos, hits));
  }
}

void read_line(ifstream &fin, list< pair<int, int> > &left_reads, list< pair<int, int> > &right_reads, bool &readable){
  string line, type_flag, qname;
  int this_chr, this_pos, bad_chr, bad_pos, hits;  
  stringstream ss;  
  if (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> hits;      
    left_or_right(left_reads, right_reads, type_flag, this_pos, hits, (this_pos < bad_pos));  
    readable = true;
  } else {
    readable = false;
  }
}

// update pos_to_scan along the way 
void read_line(ifstream &fin, list< pair<int, int> > &left_reads, list< pair<int, int> > &right_reads, bool &readable, list<int> &pos_to_scan){
  string line, type_flag, qname;
  int this_chr, this_pos, bad_chr, bad_pos, hits;  
  stringstream ss;  
  if (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> hits;      
    left_or_right(left_reads, right_reads, type_flag, this_pos, hits, (this_pos < bad_pos));  
    if ( (!pos_to_scan.size()) or this_pos > pos_to_scan.back() ) pos_to_scan.push_back( this_pos );
    readable = true;
  } else {
    readable = false;
  }
}

void check_this_pos(list <pair<int, int> > &left_reads, list <pair<int, int> > &right_reads, int this_pos, int &lr_num, int &rr_num){
  list <pair<int, int> >::iterator ri;
  ri = left_reads.begin();
  while (ri != left_reads.end() ) {
    //cerr << (*ri).first << " " << this_pos << endl;
    if ( (*ri).first >= this_pos) break;
    if ( (*ri).first < this_pos - SCAN_WIN_LEN ) { left_reads.erase(ri++); }
    else {lr_num++;  ri++;}
  }  
  ri = right_reads.begin(); 
  while (ri != right_reads.end() ) {
    //cerr << (*ri).first << " " << this_pos << endl;
    if ( (*ri).first > this_pos + SCAN_WIN_LEN - DEFAULT_READ_LEN) { ri++; continue; }
    if ( (*ri).first <= this_pos) {right_reads.erase(ri++); }
    else {rr_num++;  ri++;}
  }
}

void read_counts_location(string file1, string file2, vector<string> &chrns){
  bool readable;
  ofstream fout(file2.c_str());
  int this_pos, next_pos, lr_num, rr_num;
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    list <pair<int, int> > left_reads, right_reads;
    list <int> pos_to_scan;    
    list <int>::iterator pts;
    ifstream fin( (file1 + "." +  *ci).c_str());
    assert(fin); // init, 700 * 30 /100 = 210 
    for (int i = 0; i < 10; i++) read_line(fin, left_reads, right_reads, readable); 
    if (!readable) break;
    for (int i = 0; i < 10; i++) read_line(fin, left_reads, right_reads, readable, pos_to_scan); 
    if (!readable) break;
    //cerr << *ci << " pos_to_scan " << pos_to_scan.size() << endl;
    while ( pos_to_scan.size() ) {
      pts = pos_to_scan.begin();
      this_pos = *pts;
      next_pos = *(++pts); 
      //cerr << "check point0 " << pos_to_scan.size() << " " << this_pos << endl;
      pos_to_scan.pop_front();
      lr_num = 0;
      rr_num = 0;      
      check_this_pos(left_reads, right_reads, this_pos, lr_num, rr_num);
      //cerr << this_pos << " 2 " << lr_num << " " << rr_num << " " << left_reads.size() << " " << right_reads.size() << endl;
      if (lr_num + rr_num >= LEFT_PLUS_RIGHT) fout << *ci << " " << this_pos  << " " << lr_num  << " " << rr_num << endl;      
      if (pos_to_scan.size() < 5 or ( (!right_reads.empty()) and right_reads.back().first < next_pos + SCAN_WIN_LEN )) {
	while ( right_reads.size() < 5 or right_reads.back().first < this_pos + SCAN_WIN_LEN ) {
	  read_line(fin, left_reads, right_reads, readable, pos_to_scan); 
	  if (!readable) break;
	}	
      }      
    }
    //cerr << "done " << *ci << endl;
    fin.close();
  }
  //cerr << "output to " << file2 << endl;
  fout.close();
}

void join_location(string file1, string file2, int pos_dif, int count_dif){ // AND
  string chrn, chrn_pre;
  int maxCount, nowCount, pos, lr_num, rr_num, pos_pre, lr_num_pre, rr_num_pre;
  vector<int>  maxCount_pos;
  ifstream fin(file1.c_str());
  assert(fin);
  ofstream fout(file2.c_str());
  //int ti = 0;
  while ( fin >> chrn >> pos >> lr_num >> rr_num) {
    //if ( ti++ > 20) break;
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
      } else if ( nowCount == maxCount) {
	maxCount_pos.push_back(pos);
      }
    } else { // create new block      
      fout << chrn << " " << maxCount << " " << maxCount_pos.front() << " " << maxCount_pos.back() << endl;
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


struct compare_row {
  bool operator()(const RowInfo& a, const RowInfo& b) const {
    return (a.first == b.first) ? (a.second < b.second) :  (a.first < b.first);
  }
};

bool compare_list(const RowInfo& a, const RowInfo& b) {
  return (a.first == b.first) ? (a.second < b.second) :  (a.first < b.first);
}

int read_sort_by_col(string fn, int coln, bool has_header, set< RowInfo, compare_row > &rows) {
  string line, tmpv;
  int pos;
  stringstream ss;
  ifstream fin( fn.c_str());
  assert(fin); 
  rows.clear();
  if (has_header) getline(fin, line);
  int rown = 0;
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
  int rown = 0;
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
  fout << header;
  vector<int>::iterator bi = alurefpos->beginV.begin();
  vector<int>::iterator ei = alurefpos->endV.begin();
  int bei = 0;
  int be_size = alurefpos->beginV.size();
  int left_pos, len_overlap;
  for (set< RowInfo, compare_row >::iterator ri = rows.begin(); ri!=rows.end(); ri++) {
    left_pos = (*ri).first;    // when to start ????
    while ( left_pos >= (*ei) and bei < be_size) { bi++; ei++; bei++; }
    if ( (len_overlap = min(*ei, left_pos + DEFAULT_READ_LEN) - max(*bi, left_pos)) > ALU_MIN_LEN ) 
      fout << (*ri).second << " " << len_overlap << endl;    
  }
  fout.close();
}

void reorder_column(string fn, MapFO &fileMap, bool has_header){
  string line, type_flag, qname;
  int this_chr, this_pos, bad_chr, bad_pos, hits;  
  stringstream ss;  
  ifstream fin( fn.c_str() );
  if (has_header) getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> bad_chr >> bad_pos >> this_chr >> this_pos >> hits;      
    if (type_flag == "left") continue;
    *(fileMap[this_chr]) << type_flag << " " << qname <<" " << this_chr <<" " << this_pos <<" " << bad_chr <<" " << bad_pos <<" " << hits << endl;
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


void join_location2(vector<string> &chrns, string prefix_if, string prefix_of, int pa_dif, int pb_dif, int min_num_pn){
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    ifstream fin( (prefix_if+*ci).c_str() );
    assert(fin);
    ofstream fout( (prefix_of+*ci).c_str() );
    int idx, pa, pb, num;
    string chrn;
    set<int> ids;
    fin >> chrn >> num >> pa >> pb >> idx;
    int pa_block = pa;
    int pb_block = pb;
    int pa_pre = pa;
    int reads_num = num;
    ids.insert(idx);
    while (fin >> chrn >> num >> pa >> pb >> idx) {
      if ( pa - pa_block <= pa_dif or abs(pb - pb_block) <= pb_dif) {  // same block
	ids.insert(idx);
	pb_block = max(pb_block, pb);
	reads_num += num;
      } else { // new block
	if ( ids.size() >= min_num_pn)
	  fout << pa_block << " " << pb_block << " " << ids.size() << " " << reads_num << endl;
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

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  int opt, idx_pn;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];

  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  chrns.push_back("chrX");
  chrns.push_back("chrY");
  string path0 = read_config(config_file, "file_alu_mate0") ;    
  string path1 = read_config(config_file, "file_alu_mate1") ;    

  boost::timer clocki;    
  clocki.restart();

  if (opt == 1) {    
    seqan::lexicalCast2(idx_pn, argv[3]);
    string pn = get_pn(read_config(config_file, "file_pn"), idx_pn);
    cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    map<int, seqan::CharString> rID_chrn;
    get_rID_chrn(bam_input, chrns, rID_chrn);
    
    MapFO fileMap;
    string file1_prefix = get_name(path0, pn, ".tmp1");
    string file2_prefix = get_name(path0, pn, ".tmp2");
    string file3_prefix = get_name(path0, pn, ".tmp3");
    
    string file1, file2, file3;
    string header = "flag qname bad_chr bad_pos this_chr this_pos bad_num_hits\n"; // bad_*: potential, need check
    // step 1, about 2h
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      file1 = file1_prefix + "." + toCString(rc->second);
      fileMap[rc->first] = new ofstream( file1.c_str() );
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << header ;
    }
    alu_mate_flag(bam_input, rID_chrn, fileMap);
    // step 2, 5 min
    string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      delete fileMap[rc->first];
      file1 = file1_prefix + "." + toCString(rc->second);
      file2 = file2_prefix + "." + toCString(rc->second);
      AluRefPos *alurefpos = new AluRefPos(file_alupos_prefix + toCString(rc->second));          
      keep_alu_mate(file1, file2, alurefpos, header);       
      delete alurefpos;      
    } 
    cerr << "searching alu mate, done\n";    
    // rewrite for counting 
    header = "flag qname this_chr this_pos bad_chr bad_pos bad_num_hits\n"; // bad_*: potential, need check 
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      file3 = file3_prefix + "." + toCString(rc->second);
      fileMap[rc->first] = new ofstream( file3.c_str() );
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << header ;      
    }    
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) 
      reorder_column(file2_prefix + "." + toCString(rc->second), fileMap, true);
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) 
      delete fileMap[rc->first];    
    // sort by col
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      file3 = file3_prefix + "." + toCString(rc->second);
      set< RowInfo, compare_row > rows;
      read_sort_by_col(file3, 4, !header.empty(), rows);      
      ofstream fout(file3.c_str());
      for (set< RowInfo, compare_row >::iterator ri = rows.begin(); ri!=rows.end(); ri++) 
	fout << (*ri).second << endl;
      fout.close();
    }
    // first pass for potential locations 
    read_counts_location(file3_prefix, file3_prefix, chrns);
  } else if (opt == 2) { // combine positions 
    // 1. scan for potential regions 2. filter out rep regions 
    bool repmaskPos_read = false;
    RepMaskPos *repmaskPos;
    map<int, string> ID_pn;
    get_pn(read_config(config_file, "file_pn"), ID_pn);
    vector<string> fns;
    vector<int> idx_pns;
    ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
    while (fin >> idx_pn) {
      idx_pns.push_back(idx_pn);
      string pn = ID_pn[idx_pn];
      string file3 = get_name(path0, pn, ".tmp3" );
      string file4 = get_name(path0, pn, ".tmp4" );
      string file5 = get_name(path0, pn, ".tmp5" );      
      /*
      join_location(file3, file4, 50, 4);
      if (!repmaskPos_read) {
	repmaskPos = new RepMaskPos(read_config(config_file, "file_repeatMask"), 100); // combine regions(dist < 100bp)
	repmaskPos_read = true;
      }
      filter_location_rep(file4, file5, repmaskPos);
      */
      
      fns.push_back(file5);
    }
    fin.close();
    //write_all_location(fns, idx_pns, chrns, path1+"tmp.insert_pos.");
    join_location2(chrns, path1+"tmp.insert_pos.", path1+"insert_pos.", 50, 20, 3);
    
  }
  cerr << "time used " << clocki.elapsed() << endl;
  return 0;  
}
