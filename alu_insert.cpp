#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

typedef map<int, ofstream* > Map;
#define ALU_MIN_LEN 50
#define DEFAULT_READ_LEN 100 // if unknown, use default
#define FLAG_LONG -1
#define DISCORDANT_LEN 1000 
#define SCAN_WIN_LEN 350  // 0.99 quantile = 400, 400-50
#define LEFT_PLUS_RIGHT 5 // minimum sum of left and right reads

bool empty_map(map<seqan::CharString, pair<int,int> > &map_chr, ofstream & fout,int rID, seqan::CharString &chrn){
  if (!map_chr.size()) return false;
  cerr << chrn << " size " << map_chr.size() << endl;
  map<seqan::CharString, pair<int,int> >::iterator mc;
  for (mc = map_chr.begin(); mc!= map_chr.end(); mc++)
    fout << "left " << mc->first << " " << rID << " " << (mc->second).first << " -1 -1 " << (mc->second).second << endl;
  map_chr.clear();
  return true;
}

//search alu_mate by flag, write coordiate file
bool alu_mate_flag( string bam_input, map<int, seqan::CharString> &rID_chrx, ofstream &fout1){
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );

  int rID= -1, ia = 0;
  map<seqan::CharString, pair<int,int> > same_chr_left; // value: (pos, count)
  map<seqan::CharString, pair<int,int> >::iterator sc;

  //bool hasAlignments = false;
  //jumpToRegion(inStream, hasAlignments, context, 10, 100000, 100000, baiIndex); // debug chr19
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( hasFlagQCNoPass(record) or hasFlagDuplicate(record) or hasFlagUnmapped(record) or hasFlagNextUnmapped(record) or (not hasFlagMultiple(record))) continue;
    //if (ic++ > 6000) break;
    //// this FLAG is broken ==> if (hasFlagSecondary(record)) writeRecord(bamStreamOut, record);
    ia++;
    if (ia % 1000000 == 0) cerr << ia << " " << same_chr_left.size() << endl;
    if (record.rID != record.rNextId) { // NB: redundant
      fout1 << "dif_chr " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rNextId << " " << record.pNext << " " << numOfBestHits(record) <<  endl;      
      continue;
    } else {
      if ( rID_chrx.find(record.rID) == rID_chrx.end() ) continue;
      if ( record.rID != rID) {
	empty_map(same_chr_left, fout1, rID, rID_chrx[rID]);
	cerr << "done with " << rID_chrx[rID] << " " << rID << endl;
	rID = record.rID;	
      }       
      int hits = numOfBestHits(record);      
      if (record.beginPos < record.pNext) {	
	if ( hits > 1) { 
	  same_chr_left[record.qName] = make_pair(record.beginPos, hits);
	} else if ( (abs(record.tLen) > DISCORDANT_LEN ) or (hasFlagNextRC(record) == hasFlagRC(record))) { 
	  same_chr_left[record.qName] = make_pair(record.beginPos, FLAG_LONG); // non-redundant
	}	
      } else if (record.beginPos > record.pNext) {
	sc = same_chr_left.find(record.qName);
	if (sc == same_chr_left.end()) {
	  if (hits > 1) fout1 << "mul_map " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << hits << endl;
	} else {
	  if ( (sc->second).second == FLAG_LONG) {
	    fout1 << "ilong_RC " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << hits << endl;
	  } else if (hits == 1) { // at least one pair is unique mapped
	    fout1 << "mul_map " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rID << " " << record.pNext << " " << (sc->second).second << endl;
	  }
	  same_chr_left.erase(sc);
	}	      
      }
    }
  }  
  seqan::close(inStream);     
  return 0;
}

void check_type(string inputfile, Map &fileMap, map<int, AluRefPos *> &rID_alurefpos){
  string type_flag, qname;
  int this_chr, this_pos, bad_chr, bad_pos, hits, len_overlap;
  map<string, int > qname_hits_difchr; 
  map<string, int >::iterator qh;
  map<int, AluRefPos *>::iterator ra;
  //cerr << "call " << inputfile << endl;
  ifstream fin(inputfile.c_str());
  assert(fin);
  getline(fin, qname);
  //  while (fin >> type_flag >> qname >> this_chr >> this_pos >> bad_chr_string >> bad_pos >> hits){
  while (fin >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> hits){
    if (type_flag == "left") continue;
    if (type_flag == "dif_chr") {
      if (this_chr < bad_chr) {
	qname_hits_difchr[qname] = hits;
      } else if ( (qh = qname_hits_difchr.find(qname)) != qname_hits_difchr.end()) {
	int pre_hits = qh->second;
	qname_hits_difchr.erase(qh);
	if ( hits > 1 and pre_hits > 1) continue;
	if ( rID_alurefpos.find(this_chr) == rID_alurefpos.end() or rID_alurefpos.find(bad_chr) == rID_alurefpos.end() ) continue;
	if ( rID_alurefpos[this_chr]->insideAlu(this_pos, this_pos + DEFAULT_READ_LEN, ALU_MIN_LEN, len_overlap) and pre_hits <= 1) {
	  *(fileMap[bad_chr]) << "dif_chr " << qname << " " << bad_chr << " " << bad_pos << " " << this_chr << " " << this_pos << " " << hits << " " << len_overlap << endl;      
	} else if ( rID_alurefpos[bad_chr]->insideAlu(bad_pos, bad_pos + DEFAULT_READ_LEN, ALU_MIN_LEN, len_overlap) and hits <= 1) {  // one end is unique
	  *(fileMap[this_chr]) << "dif_chr " << qname << " " <<  this_chr << " " << this_pos << " " << bad_chr << " " << bad_pos << " " << pre_hits << " " << len_overlap << endl;
	}
      }
    } else{
      if ( (ra = rID_alurefpos.find(this_chr)) == rID_alurefpos.end() ) continue;
      if ( type_flag == "mul_map" and (ra->second)->insideAlu(bad_pos, bad_pos + DEFAULT_READ_LEN, ALU_MIN_LEN, len_overlap) ) { // hits > 1
	*(fileMap[this_chr]) << "mul_map " << qname << " " <<  this_chr << " " << this_pos << " " << bad_chr << " " << bad_pos << " " << hits << " " << len_overlap << endl;
      } else if (type_flag == "ilong_RC") {  // hits = 1
	if ( (ra->second)->insideAlu(bad_pos, bad_pos + DEFAULT_READ_LEN, ALU_MIN_LEN, len_overlap) )
	  *(fileMap[this_chr]) << "ilong_RC " << qname << " " << this_chr << " " << this_pos << " " << bad_chr << " " << bad_pos << " " << hits << " " << len_overlap << endl;
	else if ( (ra->second)->insideAlu(this_pos, this_pos + DEFAULT_READ_LEN, ALU_MIN_LEN, len_overlap) )
	  *(fileMap[this_chr]) << "ilong_RC " << qname << " " << bad_chr << " " << bad_pos << " " << this_chr << " " << this_pos << " " << hits << " " << len_overlap << endl;       
      }
    }
  }  
  fin.close();
}

void add_to_list(list< pair<int, int> > &left_reads, list< pair<int, int> > &right_reads, string &type_flag, int pos, int hits, bool mapped_left) {
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
    add_to_list(left_reads, right_reads, type_flag, this_pos, hits, (this_pos < bad_pos));  
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
    add_to_list(left_reads, right_reads, type_flag, this_pos, hits, (this_pos < bad_pos));  
    if ( this_pos > pos_to_scan.back() ) pos_to_scan.push_back( this_pos );
    readable = true;
  } else {
    readable = false;
  }
}

void check_this_pos(list <pair<int, int> > &left_reads, list <pair<int, int> > &right_reads, int this_pos, int lr_num, int rr_num){
  list <pair<int, int> >::iterator ri;
  for ( ri = left_reads.begin(); ri != left_reads.end(); ri++) {
    if ( (*ri).first >= this_pos) break;
    if ( (*ri).first < this_pos - SCAN_WIN_LEN ) left_reads.erase(ri);
    else lr_num++;
  }
  for ( ri = right_reads.begin(); ri != right_reads.end(); ri++) {
    if ( (*ri).first > this_pos + SCAN_WIN_LEN - DEFAULT_READ_LEN) continue;
    if ( (*ri).first <= this_pos) right_reads.erase(ri);
    else rr_num++;
  }
}

void scan_location(string file1, string file2, vector<string> &chrns){
  bool readable;
  ofstream fout(file2.c_str());
  int this_pos, lr_num, rr_num;
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    list <pair<int, int> > left_reads, right_reads;
    list <int> pos_to_scan;    
    list <int>::iterator pts;
    ifstream fin( (file1 + *ci).c_str());
    assert(fin); // init, 700 * 30 /100 = 210 
    for (int i = 0; i < 100; i++) read_line(fin, left_reads, right_reads, readable); 
    if (!readable) break;
    for (int i = 0; i < 100; i++) read_line(fin, left_reads, right_reads, readable, pos_to_scan); 
    if (!readable) break;
    while (pos_to_scan.size()) { // read and scan in turn
      pts = pos_to_scan.begin();
      this_pos = *pts;
      lr_num = 0;
      rr_num = 0;
      check_this_pos(left_reads, right_reads, this_pos, lr_num, rr_num);
      if (lr_num + rr_num >= LEFT_PLUS_RIGHT) fout << *ci << " " << this_pos  << " " << lr_num  << " " << rr_num << endl;
      if (pos_to_scan.size() < 3 or right_reads.back().first < *(++pts) + SCAN_WIN_LEN ) {
	while ( right_reads.back().first < this_pos + SCAN_WIN_LEN ) read_line(fin, left_reads, right_reads, readable, pos_to_scan); 
	if (!readable) break;
      }
    }
    fin.close();
  }
  fout.close();
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);

  int opt, idx_pn;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];
  seqan::lexicalCast2(idx_pn, argv[3]);

//  unsigned coverage_max;
//  seqan::lexicalCast2(coverage_max, (read_config(config_file, "coverage_max")));

  string pn = get_pn(read_config(config_file, "file_pn"), idx_pn);
  cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
  string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";

  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  chrns.push_back("chrX");
  chrns.push_back("chrY");
  map<int, seqan::CharString> rID_chrx;
  get_rID_chrx(bam_input, chrns, rID_chrx);
  string file1 = read_config(config_file, "file_alu_mate0") + pn;    
  string file1_prefix = read_config(config_file, "file_alu_mate1");
  string file2 = file1_prefix + "unsort/" + pn;    
  string file3 = file1_prefix + "sort_all/" + pn;    
  string file4 = file1_prefix + "location_all/" + pn;      
  string file5 = file1_prefix + "location_noRep/" + pn;      
  boost::timer clocki;    
  clocki.restart();

  if (opt == 1) {    // write discordant reads 
    ofstream fout1( file1.c_str() );
    fout1 << "flag qname this_chr this_pos bad_chr bad_pos num_hits\n"; // bad_*: potential, need check
    alu_mate_flag(bam_input, rID_chrx, fout1);
    fout1.close();
  } else if (opt == 2) {  // filter reads with alu_mate
    string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
    map<int, AluRefPos *> rID_alurefpos;
    Map fileMap;
    for (map<int, seqan::CharString>::iterator rc = rID_chrx.begin(); rc != rID_chrx.end(); rc++) {
      rID_alurefpos[rc->first] = new AluRefPos(file_alupos_prefix + toCString(rc->second), true);    
      fileMap[rc->first] = new ofstream( (file2 + "." + toCString(rc->second)).c_str() );
      assert(fileMap[rc->first]);
    }    
    check_type(file1, fileMap, rID_alurefpos);       
    for (Map::iterator fm = fileMap.begin(); fm != fileMap.end(); fm++) {
      delete rID_alurefpos[fm->first];
      delete fm->second;
    }
    rID_alurefpos.clear();
    cout << "./alu_insert.sh " << file2 << " " << file3  << endl;   // print cmds for sorting 
  } else if (opt == 3) { // scan for potential regions 
    scan_location(file3, file4, chrns);
    //filter_location_rep(file4, file5);
  }  
  
  cerr << "time used " << clocki.elapsed() << endl;
  return 0;  
}
