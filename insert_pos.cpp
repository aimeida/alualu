#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

#define MAJOR_SPLIT_POS_FREQ 0.6
#define FLANK_REGION_LEN 80
typedef map <int, int> SplitPosInfo;

inline string get_name_pos1(string prefix, string chrn){ return prefix + chrn + ".highCov"; }

bool findSplit(int & beginPos,int &jump_len, seqan::BamAlignmentRecord &record) {
  char operation;
  int count;
  if ( length(record.cigar) < 3) return false;
  beginPos = record.beginPos;
  jump_len = 0;
  for (size_t i = 0; i < length(record.cigar); i++) {
    operation = record.cigar[i].operation;
    count = record.cigar[i].count;
    if ( operation == 'M' or operation == 'D') //  ??? or operation == 'S') 
      beginPos += count;
    if ( operation == 'N') {
      jump_len = count;
      break;
    }
  }
  return jump_len;
}

// move split pos further to ref_begin
void moveSplit_toRefBegin(int  & split_begin, TSeq const & ref, int nnn_begin, int jump_len){  
  while (split_begin < nnn_begin and split_begin + jump_len < (int)length(ref) and ref[split_begin] == ref[ split_begin + jump_len] ) 
    split_begin++;
} 

bool coverage_idx_pass(string &line, int &ref_begin, int &ref_end, int idx_pn_this, float freq_min, float freq_max, int pn_cnt, bool & maf_pass, stringstream & ss_highCov){
  int n_reads, idx_pn;
  float coverage;
  stringstream ss;
  ss.str(line);
  int n_pn = 0, flag = 0;
  ss >> n_reads >> ref_begin >> ref_end;
  while ( ss >> idx_pn) {
    n_pn ++;
    if (idx_pn_this == idx_pn) flag = 1;
  }
  coverage = n_reads * DEFAULT_READ_LEN / (float) n_pn / (float)(ref_end - ref_begin + SCAN_WIN_LEN * 2);
  if ( idx_pn_this == 0 and coverage >= INS_COVERAGE_MAX ) 
    ss_highCov << ref_begin << " " << ref_end << endl;

  float n_freq = (float)n_pn / pn_cnt; 
  maf_pass = (n_freq <= freq_max and n_freq >= freq_min);
  return maf_pass and coverage < INS_COVERAGE_MAX and flag ;
}

// write all reads that is useful for split mapping or inferring insert sequence !
bool readbam_loci(string chrn, seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, size_t region_begin, int region_end, ofstream &frecord) {
  
  int ref_begin, ref_end; // range for broken reads
  refPos_by_region(region_begin, region_end, ref_begin, ref_end);  
  bool hasAlignments = false;  
  if (!jumpToRegion(inStream, hasAlignments, context, rID, ref_begin, ref_end, baiIndex)) return false;
  if (!hasAlignments) return false;
  seqan::BamAlignmentRecord record;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( record.rID != rID || record.beginPos >= ref_end) break;
    int read_end = record.beginPos + getAlignmentLengthInRef(record);
    if ( read_end <= ref_begin) continue; 
    if ( ! QC_insert_read(record) ) continue;
    if ( has_soft_last(record, CLIP_BP) or has_soft_first(record, CLIP_BP) ) 
      frecord << chrn << " " << region_begin << " " << region_end << " " << record.beginPos << " " << read_end 
	      << " " << get_cigar(record) << " " << hasFlagRC(record) << " " << record.seq << endl;
  } 
  return true;
}

void reads_insert_loci(int idx_pn, string chrn, string bam_input, string bai_input, string fin_pos, string fout_reads_fa, float freq_min, float freq_max, ofstream &fout_pos, int pn_cnt){

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
  map<int, seqan::CharString> rID_chrn;
  vector<string> chrns;
  chrns.push_back(chrn);
  get_rID_chrn(bam_input, chrns, rID_chrn);

  ofstream frecord(fout_reads_fa.c_str());
  frecord << "chrn region_begin region_end beginPos endPos cigar hasFlagRC seq\n" ;
  int region_begin, region_end;
  string line;
  for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
    string chrn = toCString(rc->second); 
    ifstream fin( (fin_pos + chrn).c_str());
    stringstream ss_highCov;
    if (!fin) continue;  // ignore random chr
    int cnt = 0;
    bool maf_pass;
    bool loci_pass;
    while (getline(fin, line)) {
      loci_pass = coverage_idx_pass(line, region_begin, region_end, idx_pn, freq_min, freq_max, pn_cnt, maf_pass, ss_highCov);
      if ( maf_pass and !idx_pn) fout_pos << line << endl;
      if ( loci_pass) readbam_loci(chrn, inStream, baiIndex, context, rc->first, region_begin, region_end, frecord);
      cnt++;	  
    }
    fin.close();
    cout << cnt << " loci at " << chrn << " are read for " << idx_pn << endl;

    if (!idx_pn) {
      ofstream fPos;
      fPos.open( get_name_pos1(fin_pos, chrn).c_str());
      fPos << ss_highCov.str();
      fPos.close();
      ss_highCov.clear();
    }
  }
  seqan::close(inStream);     
}

void read_file_pn_used(string fn, set<string> & pns_used) {
  ifstream fin(fn.c_str()) ;
  assert(fin);
  string pn;
  while (fin >> pn) pns_used.insert(pn);
  fin.close();
}

void read_file_pn_used(string fn, set<int> & ids_used, map<string, int> & pn_ID) {
  ifstream fin(fn.c_str()) ;
  assert(fin);
  string pn;
  while (fin >> pn) ids_used.insert(pn_ID[pn]);
  fin.close();
}

void parse_cigar(string cigar, list <char> & opts, list <int> & cnts){
  string cnt;
  int cnt_int;
  opts.clear();
  cnts.clear();
  for (size_t i = 0; i < cigar.size(); i++) {
    if ( !isdigit(cigar[i]) ) {
      opts.push_back(cigar[i]);
      if (!cnt.empty()) {
	seqan::lexicalCast2(cnt_int, cnt);
	cnts.push_back(cnt_int);
      }
      cnt = "";
    } else {
      cnt += cigar[i];
    }
  }
  seqan::lexicalCast2(cnt_int, cnt);
  cnts.push_back(cnt_int);  
}

float major_key_freq(map <int, int> &m, int & key)  {
  if ( m.empty() ) return 0;
  int max_count = 0, sum_count = 0;
  key = m.begin()->first;
  for (map <int, int>::iterator it = m.begin(); it != m.end(); it++) { 
    sum_count +=  it->second;
    if ( it->second > max_count ) {
      max_count = it->second;
      key = it->first;
    }
  }
  //cout << "mk_freq " << max_count/(float) sum_count << endl;
  return max_count/(float) sum_count;
}

float left_major_pos(map<int, int> &m, int pos_dif, int &pos) {
  if ( m.empty() ) return 0;
  map <int, int>::iterator it = m.begin();
  pos = it->first;
  int max_count = m[pos];
  int sum_count = max_count;
  for (; it != m.end(); it++) { 
    sum_count +=  it->second;
    if ( it->first <= pos + pos_dif) max_count += it->second;    
  }
  return max_count/(float) sum_count;
}

float right_major_pos(map<int, int> &m, int pos_dif, int &pos) {
  if ( m.empty() ) return 0;
  int max_count, sum_count;
  map <int, int>::iterator it;
  for (it = m.begin(); it != m.end(); it++) { pos = it->first; max_count = it->second;}
  sum_count = max_count;
  size_t i;
  for (i = 0, it = m.begin(); i < m.size() - 1; i++, it++) {
    sum_count +=  it->second;
    if ( pos - it->first <= pos_dif ) max_count += it->second;
  }
  return max_count/(float) sum_count;
}

bool find_insert_pos(string file_clip, int region_begin, int region_end, int flanking_len, float major_freq, int & clipLeft, int & clipRight, set <string> &pn_set){
  pn_set.clear();
  stringstream ss;
  string line, pn, cigar, seq;
  int beginPos, endPos, hasFlagRC; 
  int ni = 0;
  map <int, int> ni_pos;
  map <int, string> ni_pn;
  ifstream fin;
  fin.open( file_clip.c_str());
  SplitPosInfo splitPosLeft, splitPosRight;
  while ( getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> pn >> beginPos >> endPos >> cigar;
    list <char> cigar_opts;
    list <int> cigar_cnts;
    parse_cigar(cigar, cigar_opts, cigar_cnts);
    if ( *cigar_opts.begin() == 'S' and *cigar_cnts.begin() >= CLIP_BP) {
      if ( beginPos >= region_begin - flanking_len and beginPos < region_end + flanking_len) {
	addKey(splitPosLeft, beginPos);
	// if ( region_begin == 84952) cout << "debug begin " << beginPos << " " << pn << endl;
	ni_pos[ni] = beginPos;
	ni_pn[ni] = pn;
      } // to add: else {add info for next region}
    } else if ( cigar_opts.back() == 'S' and cigar_cnts.back() >= CLIP_BP) {
      if ( endPos >= region_begin - flanking_len and endPos < region_end + flanking_len) {
	addKey(splitPosRight, endPos);
	// if ( region_begin == 84952) cout << "debug end " << endPos <<  " " << pn << endl;	  
	ni_pos[ni] = endPos;
	ni_pn[ni] = pn;
      }
    }        
    ni++;
  }
  fin.close();
  
  //if ( region_begin == ) cout << "debug " << splitPosLeft.size() << " " << splitPosRight.size() << endl;  
  int pos1, pos2;
  clipLeft = 0; clipRight = 0;
  if ( major_key_freq(splitPosLeft, pos1) >= major_freq or left_major_pos(splitPosLeft, CLIP_BP, pos1) >= major_freq )
    clipLeft = pos1;  
  if ( major_key_freq(splitPosRight, pos2) >= major_freq or right_major_pos(splitPosRight, CLIP_BP, pos2) >= major_freq )
    clipRight = pos2;
  
  if (!clipLeft) clipLeft = clipRight;
  if (!clipRight) clipRight = clipLeft;
  if ((!clipLeft) or clipRight - clipLeft > CLIP_BP) return false;
  
  map <int, int>::iterator si = ni_pos.begin();
  map <int, string>::iterator pi = ni_pn.begin();  
  for ( size_t i = 0; i < ni_pos.size(); i++ ) {
    if ( abs( si->second - clipLeft ) <= CLIP_BP or  abs( si->second - clipRight ) <= CLIP_BP )
      pn_set.insert(pi->second);
    si++;
    pi++;
  }
  // ignore private mutations 
  return pn_set.size() > 1;  
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);

  boost::timer clocki;    
  clocki.restart();  

  int opt;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];

  string path1 = read_config(config_file, "file_alu_insert1") ;    
  string fout_path = read_config(config_file, "file_clip_reads");
  check_folder_exists(fout_path);
  vector<string> chrns;
  string tmp_path;
  for (int i = 1; i < 23; i++) {
    string chrn = "chr" + int_to_string(i);
    chrns.push_back(chrn);
    tmp_path = fout_path + chrn + "_pos/";
    check_folder_exists( tmp_path);
    tmp_path = fout_path + chrn + "/";
    check_folder_exists( tmp_path );      
  }

#ifdef DEBUG_MODE   // only look at chr1 for now 
    chrns.clear();
    chrns.push_back("chr1");
#endif       

  map<int, string> ID_pn;
  get_pn(read_config(config_file, "file_pn"), ID_pn);      
  string file_pn_used_prefix = path1 + "pn.insert_pos.";
  string fin_pos = path1 + "insert_pos.";
  float freq_min = 0.02, freq_max = 1;  // small groups 
  stringstream ss;
  ss << "_" << setprecision(3) << freq_min << "_" << setprecision(3) << freq_max;
  string fout_suffix = ss.str();
  ss.clear();
  
  if (opt == 1) { // less than 2 to 60 mins for each pn.
    int idx_pn;
    seqan::lexicalCast2(idx_pn, argv[3]);
    assert (argc == 4);
    string pn = ID_pn[idx_pn];    
    cout << "reading " << pn << endl;
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";      
    ofstream fout_pos;
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      set<string> pns_used;
      read_file_pn_used(file_pn_used_prefix + *ci, pns_used); // some pn are not used 
      // if this pn is ignored, due to too many reads
      if ( idx_pn and pns_used.find(pn) == pns_used.end() ) continue; 
      if (!idx_pn)  fout_pos.open( (fin_pos + *ci + fout_suffix).c_str());
      string fout_reads_fa = fout_path + *ci + "/" + pn + fout_suffix;
      reads_insert_loci(idx_pn, *ci, bam_input, bai_input, fin_pos, fout_reads_fa, freq_min, freq_max, fout_pos, ID_pn.size() );
      if (!idx_pn) fout_pos.close();
    }
    ////writeRecord(fout, record.qName, record.seq, record.qual, seqan::Fastq());
    ////writeRecord(fout, record.qName, record.seq, seqan::Fasta());
  } else if ( opt == 2 ) { 
    fstream fout;
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      set<string> pns_used;
      set<string> beginPos_fn;
      read_file_pn_used(file_pn_used_prefix + *ci, pns_used);
      for (set<string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) {
	string file_input_clip, file_output_clip;
	file_input_clip = fout_path + *ci + "/" + *pi + fout_suffix;
	ifstream fin(file_input_clip.c_str());
	assert(fin);
	stringstream ss;
	string line, tmp1, cigar, seq, region_begin, region_end, beginPos, endPos, hasFlagRC;
	getline(fin, line);
	string region_begin_pre="";
	while ( getline(fin, line) ) {
	  ss.clear(); ss.str(line); 
	  ss >> tmp1 >> region_begin >> region_end >> beginPos >> endPos >> cigar >> hasFlagRC >> seq;
	  file_output_clip = fout_path + *ci + "_pos/" + region_begin + "_" + region_end;
	  if (region_begin != region_begin_pre) {
	    region_begin_pre = region_begin;
	    fout.close();
	    if ( beginPos_fn.find(region_begin) == beginPos_fn.end() ) {
	      fout.open(file_output_clip.c_str(), fstream::out);
	      beginPos_fn.insert(region_begin);
	    } else {
	      fout.open(file_output_clip.c_str(), fstream::app|fstream::out);
	    }
	  }
	  fout << *pi << " " << beginPos << " " << endPos << " " << cigar << " " << hasFlagRC << " " << seq << endl;
	}	  
	fin.close();
	//cout << "read " << *pi << endl;       	
      }
    }    
  } else if ( opt == 3 ) {  // find exact split position

    cerr << "NB: private insertion is ignored here !\n";
    string cnt, line;
    int region_begin, region_end;
    stringstream ss;
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      ifstream fin( (fin_pos + *ci + fout_suffix).c_str());
      assert (fin);
      ofstream fout( (path1 + "clip/" + *ci + fout_suffix).c_str() ); 
      while ( getline(fin, line) ) {
	ss.clear(); ss.str(line); 
	ss >> cnt >> region_begin >> region_end;
	string file_clip = fout_path + *ci + "_pos/" + int_to_string(region_begin) + "_" + int_to_string(region_end);
	int clipLeft, clipRight;
	//int clipLeft_next, clipRight_next;
	set <string> pn_set;
	if (find_insert_pos(file_clip, region_begin, region_end, FLANK_REGION_LEN, MAJOR_SPLIT_POS_FREQ, clipLeft, clipRight, pn_set)) {
	  fout << region_begin << " " << clipLeft << " " << clipRight;
	  for (set<string>::iterator pni = pn_set.begin(); pni != pn_set.end(); pni++)
	    fout << " " << *pni ;
	  fout << endl;
	}
	//if ( region_begin == 226767) exit(0);	  
      }
      fin.close();
      fout.close();
    }
  }
  
  return 0;
}
