#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "utils.h"
#include "insert_utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

#define ALU_MIN_LEN 50  // min overlap in alu region
#define DEFAULT_READ_LEN 100 // if unknown, use default
#define CLIP_BP 10
#define FLAG_LONG -1
#define FLAG_RC -2
#define DISCORDANT_LEN 1000 
#define INS_COVERAGE_MAX 6    // approximate, not exact, remove some high coverage regions 
#define SCAN_WIN_LEN 400  // 0.99 quantile = 400
#define LEFT_PLUS_RIGHT 5 // minimum sum of left and right reads


inline string get_name(string path, string fn, string suffix){ return path + fn + suffix;}

inline void region_to_ref_pos(int region_begin, int region_end, int &ref_begin, int &ref_end){
  ref_begin = region_begin - DEFAULT_READ_LEN;
  ref_end = region_end + 2 * DEFAULT_READ_LEN;
}

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
	  cerr << "left " << mc->first << " " << rID_pre << " " << (mc->second).first << " -1 -1 " << (mc->second).second << " " << length(record.seq) << endl;
	same_chr_left.clear();
      }
      cerr << "done with " << rID_chrn[rID_pre] << " " << rID_pre << endl;
      rID_pre = record.rID;	
    }       
    if (ia++ % 1000000 == 0) cerr << ia << " " << same_chr_left.size() << endl;
    if (record.rID != record.rNextId) { // NB: redundant! one pair read twice, thus are written twice as well
      if (rID_chrn.find(record.rNextId) != rID_chrn.end())  // both chrn exist
	*(fileMap[record.rNextId]) << "dif_chr " << record.qName << " " << record.rNextId << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << numOfBestHits(record) << " " << length(record.seq) << endl;      
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
	  *(fileMap[record.rID]) << "mul_map " << record.qName << " " << record.rID <<  " " << record.beginPos << " " << record.rID << " " << record.pNext << " " << hits << " " << length(record.seq) << endl;
      } else {                   // left read not multi_mapped
	if ( (sc->second).second == FLAG_LONG) {  // NB, redundant! written twice when right reads are read
	  *(fileMap[record.rID]) << "ilong_RC " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " 1" << hits << " " << length(record.seq) << endl;
	  *(fileMap[record.rID]) << "ilong_RC " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rID << " " << record.pNext << " " << hits << " " << length(record.seq) << endl;
	} else if (hits == 1) {  // at least one pair is ok
	  *(fileMap[record.rID]) << "mul_map " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << (sc->second).second << " " << length(record.seq)  << endl;
	}
	same_chr_left.erase(sc); // rm left read
      }
     
    }    
  } 
  seqan::close(inStream);     
  return 0;
}

void read_line(ifstream &fin, list< READ_INFO *> &lr_reads, bool &readable){
  string line, type_flag, qname, alu_type;
  int this_chr, this_pos, bad_chr, bad_pos, hits, len_read;  
  stringstream ss;  
  readable = false;
  if (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> hits >> len_read >> alu_type;      
    // ignore mul_map, ilong_RC for now.
    if (type_flag == "dif_chr" or abs(this_pos - bad_pos) >= 2 * DISCORDANT_LEN) 
      lr_reads.push_back(new READ_INFO(this_pos, len_read, alu_type));
    readable = true;
  } 
}

void read_line(ifstream &fin, list< READ_INFO *> &lr_reads, bool &readable, list<int> &pos_to_scan){
  string line, type_flag, qname, alu_type;
  int this_chr, this_pos, bad_chr, bad_pos, hits, len_read;  
  stringstream ss;  
  readable = false;
  if (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos >> hits >> len_read >> alu_type;      
    // ignore mul_map, ilong_RC for now.
    if (type_flag == "dif_chr" or abs(this_pos - bad_pos) >= 2 * DISCORDANT_LEN) 
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

void alumate_counts_filter(string file1, string file2, vector<string> &chrns){
  bool readable;
  ofstream fout(file2.c_str());
  int this_pos, next_pos, lr_num, rr_num;
  map <string, int> alu_type_count;
  string line;
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    list< READ_INFO *> lr_reads;
    list <int> pos_to_scan;    
    list <int>::iterator pts;
    ifstream fin( (file1 + "." +  *ci).c_str());
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
      //cerr << "check point0 " << pos_to_scan.size() << " " << this_pos << endl;
      pos_to_scan.pop_front();
      lr_num = 0;
      rr_num = 0;      
      check_this_pos(lr_reads, this_pos, lr_num, rr_num, alu_type_count);
      if (lr_num + rr_num >= LEFT_PLUS_RIGHT) 
	fout << *ci << " " << this_pos  << " " << lr_num  << " " << rr_num << " " <<  major_type(alu_type_count) <<  endl;      
      ///if (pos_to_scan.size() < 5 or ( (!right_reads.empty()) and right_reads.back().first < next_pos + SCAN_WIN_LEN )) {
      while (pos_to_scan.size() < 5 or (lr_reads.back())-> beginPos < next_pos + SCAN_WIN_LEN ) {
	read_line(fin, lr_reads, readable, pos_to_scan); 
	if (!readable) break;
      }	
    }
    //cerr << "done " << *ci << endl;
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
    left_pos = (*ri).first;    // when to start ????
    //cout << ni << " " << rows.size() << " " << bei << " " << be_size <<  endl;
    while ( ( left_pos >= (*ei)) and bei < be_size) { bi++; ei++; ti++; bei++; }
    //cout << ni << " " << rows.size() << " " << (*ri).second << " " << *bi << " " << *ei << " " << *ti << " " <<  bei << " " << be_size <<  endl;
    if ( (len_overlap = min(*ei, left_pos + DEFAULT_READ_LEN) - max(*bi, left_pos)) > ALU_MIN_LEN ) 
      fout << (*ri).second << " " << *ti << endl;    
  }
  fout.close();
}

void reorder_column(string fn, MapFO &fileMap, bool has_header){
  string line, type_flag, qname, sub_type;
  int this_chr, this_pos, bad_chr, bad_pos, hits, len_read;  
  stringstream ss;  
  ifstream fin( fn.c_str() );
  if (has_header) getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> bad_chr >> bad_pos >> this_chr >> this_pos >> hits >> len_read >> sub_type;      
    if (type_flag == "left") continue;
    *(fileMap[this_chr]) << type_flag << " " << qname <<" " << this_chr <<" " << this_pos <<" " << bad_chr <<" " << bad_pos <<" " << hits << " " << len_read << " " << sub_type << endl;
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

bool consider_this_loci(string &line, int &ref_begin, int &ref_end, int &idx_pn_first, int idx_pn_this){
  int n_pn, n_reads, idx_pn;
  float coverage;
  string alu_type;
  stringstream ss;
  ss.str(line);
  ss >> ref_begin >> ref_end >> n_pn >> n_reads >> alu_type >> idx_pn_first;
  coverage = n_reads * DEFAULT_READ_LEN / (float) n_pn / (float)(ref_end - ref_begin + SCAN_WIN_LEN * 2);
  if ( coverage >= INS_COVERAGE_MAX ) return false;
  if ( idx_pn_first ==  idx_pn_this) return true;
  for (int i = 0; i < n_pn - 1; i++) {
    ss >> idx_pn;
    if (idx_pn_this == idx_pn) return true;
  }
  return false;
}

bool consider_this_loci(string &line, int &ref_begin) {
  int n_pn, n_reads, ref_end;
  float coverage;
  stringstream ss;
  ss.str(line);
  ss >> ref_begin >> ref_end >> n_pn >> n_reads;
  coverage = n_reads * DEFAULT_READ_LEN / (float) n_pn / (float)(ref_end - ref_begin + SCAN_WIN_LEN * 2);
  if ( coverage < INS_COVERAGE_MAX ) return true;
  return false;
}

void write_fastq_fr(string fout_reads, string fin_read_bam, string fin_read_map, int idx_pn, map <string, set< pair<int, int> > > & pos_of_clippedReads, bool &rm_old_fasta){
  // NB: clear first, otherwise error in append/write files
  if (rm_old_fasta) {
    system( ("rm " + fout_reads + "*fastq" ).c_str() );  
    rm_old_fasta = false;
  }
  int pre_pos = 0;
  int region_begin, region_end, idx_pn_first;
  string chrn, line; 
  stringstream ss;
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, fin_read_bam.c_str(), "r"));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );
  ofstream fastq_f, fastq_r;
  string fastq_fn;
  ifstream fin(fin_read_map.c_str());
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    getline(fin, line);
    ss.clear(); ss.str(line);
    ss >> chrn >> region_begin >> region_end >> idx_pn_first;
    if ( (!pre_pos) or (region_begin != pre_pos ) ) { 
      if (pre_pos) {
	fastq_f.close();
	fastq_r.close();
      }
      fastq_fn = fout_reads + chrn + "_"+ int_to_string(region_begin) + "_"+ int_to_string(region_end);
      pos_of_clippedReads[chrn].insert( make_pair(region_begin, region_end) );
      if (idx_pn_first == idx_pn) {
	fastq_f.open( ( fastq_fn + "_f.fastq").c_str() );
	fastq_r.open( ( fastq_fn + "_r.fastq").c_str() );
      } else {		
	// sometimes works as ios::write, if no clipped reads in previous PN.
	// fixme: easy to change, use map to keep track of reads that has been written into fasta
	fastq_f.open( ( fastq_fn + "_f.fastq").c_str(), ios::app); 
	fastq_r.open( ( fastq_fn + "_r.fastq").c_str(), ios::app);
      }
    }
    pre_pos = region_begin;
    if ( hasFlagRC(record) ) { 
      reverseComplement(record.seq);  // only necessary if i want to cut end reads 
      reverse(record.qual);
      writeRecord(fastq_r, record.qName, record.seq, record.qual, seqan::Fastq());
    } else {
      writeRecord(fastq_f, record.qName, record.seq, record.qual, seqan::Fastq());
    }
  }
  fastq_f.close();
  fastq_r.close();
  fin.close();
  seqan::close(inStream);
}

// this version only counts fastq files 
void get_filename_fastq(string fin_read_map, map <string, set< pair<int, int> > > & pos_of_clippedReads){
  int pre_pos = 0;
  string chrn, line;
  int region_begin, region_end, idx_pn_first;
  stringstream ss;
  ifstream fin(fin_read_map.c_str());
  while (getline(fin, line) ) {
    ss.clear(); ss.str(line);
    ss >> chrn >> region_begin >> region_end >> idx_pn_first;
    if ( (!pre_pos) or (region_begin != pre_pos ) ) 
      pos_of_clippedReads[chrn].insert( make_pair(region_begin, region_end) );
    pre_pos = region_begin;
  }
  fin.close();
}

void get_alu_type(string fn, map<int, string> &pos_aluType){
  stringstream ss;
  int region_begin, region_end;
  string line, tmp1, tmp2, alu_type;   
  pos_aluType.clear();
  ifstream fin(fn.c_str());
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> region_begin >> region_end >> tmp1 >> tmp2 >> alu_type;
    pos_aluType[region_begin] = alu_type;
  }
  fin.close();
}

void write_ref(string file_ref, string fa_name, TSeq &nnn, TSeq &fa_seq, TSeq alu_seq, bool alu_is_left) {
  TSeq fa_seq2;
  ofstream fout(file_ref.c_str());
  if (alu_is_left) {
    fa_seq2 = alu_seq;
    fa_seq2 += nnn;
    fa_seq2 += fa_seq;
    writeRecord(fout, "pl_" + fa_name, fa_seq2, seqan::Fasta());
    seqan::reverseComplement(alu_seq);
    fa_seq2 = alu_seq;
    fa_seq2 += nnn;
    fa_seq2 += fa_seq;
    writeRecord(fout, "nl_" + fa_name, fa_seq2, seqan::Fasta());
  } else {
    // Positive chain, alu Right
    fa_seq2 = fa_seq;
    fa_seq2 += nnn;
    fa_seq2 += alu_seq;
    writeRecord(fout, "pr_" + fa_name, fa_seq2, seqan::Fasta());
    seqan::reverseComplement(alu_seq);
    fa_seq2 = fa_seq;
    fa_seq2 += nnn;
    fa_seq2 += alu_seq;
    writeRecord(fout, "nr_" + fa_name, fa_seq2, seqan::Fasta());   
  }
  fout.close();
}

// only support 4 types for now 
string parse_alu_type(string alu_name){
  assert ( !alu_name.empty() );
  if ( alu_name.substr(0,4) == "AluY") return "AluY";
  if ( alu_name.substr(0,4) == "AluS") return "AluSx";
  if ( alu_name.substr(0,4) == "AluJ") {
    if ( alu_name.substr(0,5) == "AluJo") return "AluJo";
    if ( alu_name.substr(0,5) == "AluJb") return "AluJb";
    return "AluJo";
  }
  return "AluY";
}

void call_splazers(string &cmd, string &bin_splazers, string param, string fn1, string fn2, string fn3) {
  cmd = bin_splazers + " " + param + " " + fn1 + " " + fn2 + " -o " + fn3;
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  int opt, idx_pn;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];

  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  string path0 = read_config(config_file, "file_alu_mate0") ;    
  string path1 = read_config(config_file, "file_alu_mate1") ;    
  boost::timer clocki;    
  clocki.restart();

  map<int, string> ID_pn;
  get_pn(read_config(config_file, "file_pn"), ID_pn);
  map<int, seqan::CharString> rID_chrn;
  
  if (opt == 1) {    
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
    
    string file1, file2, file3;
    string header = "flag qname bad_chr bad_pos this_chr this_pos bad_num_hits len_read"; // bad_*: potential, need check
    /*
    // step 1, about 2h
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      file1 = file1_prefix + "." + toCString(rc->second);
      fileMap[rc->first] = new ofstream( file1.c_str() );
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << header << endl;
    }
    alu_mate_flag(bam_input, rID_chrn, fileMap);

   // step 2, if mate mapped to alu. 5 min. 
    string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      //delete fileMap[rc->first];
      string chrn = toCString(rc->second);
      file1 = file1_prefix + "." + chrn;
      file2 = file2_prefix + "." + chrn;
      AluRefPos *alurefpos = new AluRefPos(file_alupos_prefix + chrn);          
      keep_alu_mate(file1, file2, alurefpos, header);       
      delete alurefpos;      
    } 
    cerr << "searching alu mate, done\n";    
    
  // rewrite for counting, switch 2 column and sort by pos 
    header = "flag qname this_chr this_pos bad_chr bad_pos bad_num_hits len_read sub_type"; // bad_*: potential, need check 
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      file3 = file3_prefix + "." + toCString(rc->second);
      fileMap[rc->first] = new ofstream( file3.c_str() );
      assert(fileMap[rc->first]);
      *(fileMap[rc->first]) << header << endl ;      
    }    
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) 
      reorder_column(file2_prefix + "." + toCString(rc->second), fileMap, true);
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) 
      delete fileMap[rc->first];    

   // sort by col
    cerr << "sorting\n";
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
     file3 = file3_prefix + "." + toCString(rc->second);
     set< RowInfo, compare_row > rows;
     read_sort_by_col(file3, 4, true, rows);      
     ofstream fout(file3.c_str());
     fout << header << endl ;      
     for (set< RowInfo, compare_row >::iterator ri = rows.begin(); ri!=rows.end(); ri++) 
       fout << (*ri).second << endl;
     fout.close();
   }
    */

    // first pass for potential locations 
    alumate_counts_filter(file3_prefix, file3_prefix, chrns);
  } else if (opt == 2) { // combine positions from multiple individuals
    string opt2 = "tmp3";
    if (argc > 3) opt2 = argv[3];  // tmp3 or tmp5
    // 1. scan for potential regions 2. filter out rep regions 
    chrns.push_back("chrX");
    chrns.push_back("chrY");
    bool repmaskPos_read = false;
    RepMaskPos *repmaskPos;
    vector<string> fns;
    vector<int> idx_pns;
    ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
    while (fin >> idx_pn) {
      idx_pns.push_back(idx_pn);
      string pn = ID_pn[idx_pn];
      string file3 = get_name(path0, pn, ".tmp3" );
      string file4 = get_name(path0, pn, ".tmp4" );
      string file5 = get_name(path0, pn, ".tmp5" );      
      fns.push_back(file5);    
      if ( opt2 == "tmp3") {
	join_location(file3, file4, 50, 4);
	if (!repmaskPos_read) {
	  repmaskPos = new RepMaskPos(read_config(config_file, "file_repeatMask"), 100); // combine regions(dist < 100bp)
	  repmaskPos_read = true;
	}
	filter_location_rep(file4, file5, repmaskPos);
      }
    }
    fin.close();
    write_all_location(fns, idx_pns, chrns, path1+"tmp.insert_pos."); 
    join_location2(chrns, path1+"tmp.insert_pos.", path1+"insert_pos.", 40, 3); //at least 3 pn has insertion   

  } else if (opt == 3) { // run each pn seperately, write reads in insert regions 
    seqan::lexicalCast2(idx_pn, argv[3]);
    string fin_ins_pos = path1 + "insert_pos.";
    string fout_ins_read = path1 + "split_mapping_clip/";   
    if ( access( fout_ins_read.c_str(), 0 ) != 0 )
      system( ("mkdir " + fout_ins_read).c_str() );    
    string pn = get_pn(read_config(config_file, "file_pn"), idx_pn);
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";  
    seqan::BamIndex<seqan::Bai> baiIndex;
    assert (!read(baiIndex, bai_input.c_str()));
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);
    seqan::BamHeader header;
    seqan::BamAlignmentRecord record;
    seqan::Stream<seqan::Bgzf> inStream;
    assert(open(inStream, bam_input.c_str(), "r"));
    assert(!readRecord(header, context, inStream, seqan::Bam()));    
    get_rID_chrn(bam_input, chrns, rID_chrn);
    string fout_ins_read_bam = fout_ins_read + pn + ".bam" ;
    string fout_ins_read_map = fout_ins_read + pn + ".map" ;
    seqan::BamStream bamStreamOut(fout_ins_read_bam.c_str(), seqan::BamStream::WRITE);
    bamStreamOut.header = header;
    ofstream fout(fout_ins_read_map.c_str());
    int region_begin, region_end, ref_begin, ref_end, idx_pn_first;
    string line;
    for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
      string chrn = toCString(rc->second);
      ifstream fin( (fin_ins_pos + chrn).c_str());
      while (getline(fin, line)) {
	if ( consider_this_loci(line, region_begin, region_end, idx_pn_first, idx_pn) ) {
	  region_to_ref_pos(region_begin, region_end, ref_begin, ref_end);
	  vector <seqan::CharString> qnames;
	  readbam_clip(inStream, baiIndex, context, rc->first, ref_begin, ref_end, bamStreamOut, qnames, CLIP_BP);
	  for ( vector <seqan::CharString>::iterator qi = qnames.begin(); qi != qnames.end(); qi++) 
	    fout << chrn << " " << region_begin << " " << region_end << " " << idx_pn_first << " " << *qi << " " << ref_begin << " " << ref_end << endl;      
	}
      }
      fin.close();
    }
    fout.close();
    seqan::close(inStream);     
    seqan::close(bamStreamOut);

  }  else if (opt == 4) { // combine reads from all pns, write reference.fa
    string fin_read = path1 + "split_mapping_clip/";
    string fout_fastq = path1 + "split_mapping_fastq/"; // fastq and fa files
    if ( access( fout_fastq.c_str(), 0 ) != 0 ) system( ("mkdir " + fout_fastq).c_str() );
    string fout_splazer = path1 + "split_mapping_splazer/"; 
    if ( access( fout_splazer.c_str(), 0 ) != 0 ) system( ("mkdir " + fout_splazer).c_str() );
    string fout_fa = path1 + "split_mapping_fa/"; // fastq and fa files
    if ( access( fout_fa.c_str(), 0 ) != 0 ) system( ("mkdir " + fout_fa).c_str() );

    // step 1, write reads info for each position 
    bool rm_old_fasta = true;
    map<string, set< pair<int, int> > > pos_of_clippedReads;    
    ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
    while (fin >> idx_pn) {
      string pn = ID_pn[idx_pn];
      string fin_read_bam = fin_read + pn + ".bam" ;
      string fin_read_map = fin_read + pn + ".map" ;
      /** fastq files of broken reads, distinguish forward and RC reads 
	  useful if the last N bp of a read should be removed while calling splazers 
      */
      write_fastq_fr(fout_fastq, fin_read_bam, fin_read_map, idx_pn, pos_of_clippedReads, rm_old_fasta);
      /* new version, write either fastq or fasta files for broken reads */      
      //write_fastq_fr(fout_fa, fin_read_bam, fin_read_map, idx_pn, pos_of_clippedReads, rm_old_fasta, "fa");
      //get_filename_fastq(fin_read_map, pos_of_clippedReads); // only to update pos_of_clippedReads
      cout << "done with " << pn << endl;
    }
    fin.close();

    // step 2, write ref_genome.fa 
    string file_alu_cons = read_config(config_file, "file_alu_cons");
    string file_fa_prefix = read_config(config_file, "file_fa_prefix");
    seqan::FaiIndex faiIndex;
    unsigned fa_idx = 0;
    TSeq nnn;
    resize(nnn, 70, 'N');
    map<string, set< pair<int, int> > >::iterator pi;
    map<int, string> pos_aluType;    
    int region_begin, ref_begin, ref_end;
    for ( pi = pos_of_clippedReads.begin(); pi != pos_of_clippedReads.end(); pi++ ) {
      string chrn = pi->first;
      get_alu_type(path1+"insert_pos." + chrn, pos_aluType);
      assert (!read(faiIndex, (file_fa_prefix + chrn + ".fa").c_str()) );
      assert (getIdByName(faiIndex, chrn, fa_idx));
      for (set< pair<int, int> >::iterator ri = (pi->second).begin(); ri != (pi->second).end(); ri++) {
	region_begin =  (*ri).first;
	region_to_ref_pos(region_begin, (*ri).second, ref_begin, ref_end);
	//string file_ref = fout_fastq + chrn + "_" + int_to_string( region_begin ) + "_" + int_to_string ((*ri).second) + ".fa";
	string file_ref = fout_fastq + chrn + "_" + int_to_string( region_begin ) + "_" + int_to_string ((*ri).second);
	bool upper = true;
	TSeq fa_seq = (TSeq)fasta_seq(faiIndex, fa_idx, ref_begin, ref_end, upper);
	TSeq alu_seq;
	string alu_type = parse_alu_type(pos_aluType[region_begin]);
	read_fasta_by_name(alu_seq, file_alu_cons, alu_type); 
	stringstream ss;
	ss << ref_begin << "_" << ref_end << "_" << alu_type << "_" << length(alu_seq);
	write_ref(file_ref + "_l.fa", ss.str(), nnn, fa_seq, alu_seq, true);
	write_ref(file_ref + "_r.fa", ss.str(), nnn, fa_seq, alu_seq, false);
      }   
    }
    
    // step3: write command files to call splazers
    string bin_splazers = read_config(config_file, "bin_splazers");
    //  -id, allow indels, -m 1, only output best hits 
    string param1 = "-id -i 80 -m 1 -tr 95"; // fixme: some read length is 120 bp 
    string param2l = "-sm 20 -ep 3 -es 1 -minG 70 -of 4";  // alu is left
    string param2r = "-sm 20 -ep 1 -es 3 -minG 70 -of 4"; 
    string cmd;    
    ofstream fout_qsub("alu_insert_chr.txt"); 
    for ( pi = pos_of_clippedReads.begin(); pi != pos_of_clippedReads.end(); pi++ ) {
      string chrn = pi->first;
      fout_qsub << "sh /nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/rp_insert/alu_insert_" << chrn << ".txt\n";
      ofstream fout_cmd( ("alu_insert_" + chrn + ".txt").c_str() );
      for (set< pair<int, int> >::iterator ri = (pi->second).begin(); ri != (pi->second).end(); ri++) {
	int region_begin = (*ri).first;
	string file_prefix = fout_fastq + chrn + "_" + int_to_string( region_begin ) + "_" + int_to_string ((*ri).second);
	string file_prefix2 = fout_splazer + chrn + "_" + int_to_string( region_begin ) + "_" + int_to_string ((*ri).second);
	call_splazers(cmd, bin_splazers, param1 + " -f " + param2l, file_prefix + "_l.fa", file_prefix + "_f.fastq", file_prefix2 + "_fl.sam");
	fout_cmd << cmd << endl;
	call_splazers(cmd, bin_splazers, param1 + " -f " + param2r, file_prefix + "_r.fa", file_prefix + "_f.fastq", file_prefix2 + "_fr.sam");
	fout_cmd << cmd << endl;
	call_splazers(cmd, bin_splazers, param1 + " -r " + param2l, file_prefix + "_l.fa", file_prefix + "_r.fastq", file_prefix2 + "_rl.sam");
	fout_cmd << cmd << endl;
	call_splazers(cmd, bin_splazers, param1 + " -r " + param2r, file_prefix + "_r.fa", file_prefix + "_r.fastq", file_prefix2 + "_rr.sam");
	fout_cmd << cmd << endl;
      }   
      //cout << "done with " << chrn << endl;
      fout_cmd.close();
    }
    fout_qsub.close();
    cout << "NB: regions without clipped reads are skipped !\n";
    cout << "run alu_insert.*.txt to continue... \n";
  } 
  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
