#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include <seqan/basic.h>
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

// move split pos further to refBegin
void moveSplit_toRefBegin(int  & split_begin, TSeq const & ref, int nnn_begin, int jump_len){  
  while (split_begin < nnn_begin and split_begin + jump_len < (int)length(ref) and ref[split_begin] == ref[ split_begin + jump_len] ) 
    split_begin++;
} 

bool coverage_idx_pass(string &line, int &refBegin, int &refEnd, int idx_pn_this, float freq_min, float freq_max, int pn_cnt, bool & maf_pass, stringstream & ss_highCov){
  int n_reads, idx_pn;
  float coverage;
  stringstream ss;
  ss.str(line);
  int n_pn = 0, flag = 0;
  ss >> n_reads >> refBegin >> refEnd;
  while ( ss >> idx_pn) {
    n_pn ++;
    if (idx_pn_this == idx_pn) flag = 1;
  }
  coverage = n_reads * DEFAULT_READ_LEN / (float) n_pn / (float)(refEnd - refBegin + SCAN_WIN_LEN * 2);
  if ( idx_pn_this == 0 and coverage >= INS_COVERAGE_MAX ) 
    ss_highCov << refBegin << " " << refEnd << endl;

  float n_freq = (float)n_pn / pn_cnt; 
  maf_pass = (n_freq <= freq_max and n_freq >= freq_min);
  return maf_pass and coverage < INS_COVERAGE_MAX and flag ;
}

// write all reads that is useful for split mapping or inferring insert sequence !
bool readbam_loci(string chrn, seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, size_t region_begin, int region_end, ofstream &frecord) {
  
  int refBegin, refEnd; // range for broken reads
  refPos_by_region(region_begin, region_end, refBegin, refEnd);  
  bool hasAlignments = false;  
  if (!jumpToRegion(inStream, hasAlignments, context, rID, refBegin, refEnd, baiIndex)) return false;
  if (!hasAlignments) return false;
  seqan::BamAlignmentRecord record;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( record.rID != rID || record.beginPos >= refEnd) break;
    int read_end = record.beginPos + getAlignmentLengthInRef(record);
    if ( read_end <= refBegin) continue; 
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

bool pn_has_clipPos( int & clipLeft, int & clipRight, vector <pair<char, int> > & splitori_pos, float major_freq) {
  map <int, int> clipLefts, clipRights;
  for ( vector <pair<char, int> >::iterator si = splitori_pos.begin(); si != splitori_pos.end(); si++) {
    if ((*si).first == 'L') addKey(clipLefts, (*si).second);
    else if ((*si).first == 'R') addKey(clipRights, (*si).second);
  }
  clipLeft = 0; clipRight = 0;
  if (major_key_freq(clipLefts, clipLeft) >= major_freq or
      major_key_freq(clipRights, clipRight) >= major_freq ) {
    if (!clipLeft) clipLeft = clipRight;
    if (!clipRight) clipRight = clipLeft;
    return clipLeft > 0 and clipRight > clipLeft and clipRight - clipLeft < INSERT_POS_GAP;
  }
  return false;
}

bool adjust_cigar(seqan::CharString & ref_fa, int refBegin, string file_clip,float major_freq, int & clipLeft, int & clipRight, map <string, int> &pn_splitCnt, int estimated_clipPos, int max_dist_to_splitPos){
  stringstream ss;
  string line, pn, cigar, read_seq;
  int beginPos, endPos, hasRCFlag, _clipLeft, _clipRight, align_len; 
  map <string, vector <pair<char, int> > > pn_splitori_pos;
  ifstream fin( file_clip.c_str());
  while ( getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> pn >> beginPos >> endPos >> cigar >> hasRCFlag >> read_seq;
    list <char> cigar_opts;
    list <int> cigar_cnts;
    parse_cigar(cigar, cigar_opts, cigar_cnts);
    if ( abs(beginPos - estimated_clipPos) >= max_dist_to_splitPos and
	 abs(endPos - estimated_clipPos) >= max_dist_to_splitPos )
      continue;
    
    if ( *cigar_opts.begin() == 'S' and *cigar_cnts.begin() >= CLIP_BP) { // splitEnd, S*M*
      align_len = read_seq.size() - *cigar_cnts.begin();
      _clipRight = beginPos;      
      int i = 1;  // i = 0 is the last match by BWA
      while ( *cigar_cnts.begin() - i >= 0 and ref_fa[_clipRight - refBegin - i] == read_seq[*cigar_cnts.begin() - i] ) i++;
      if ( *cigar_cnts.begin() - (i-1) < CLIP_BP) 
	continue;
      _clipRight -= (i-1);
      pn_splitori_pos[pn].push_back( make_pair('R', _clipRight) );

      //cout << "E " << _clipRight << " " << i-1 << " " << cigar << " " << estimated_clipPos << " " << beginPos << " " << endPos << " " << pn << endl;
      /*
      if ( i > 1) { 	
	cout << infix(ref_fa, _clipRight - refBegin, _clipRight - refBegin + align_len) << endl;
	cout << read_seq.substr(*cigar_cnts.begin(), align_len) << endl ;
      }    */
    }

    if ( cigar_opts.back() == 'S' and cigar_cnts.back() >= CLIP_BP) {    // M*S*
      align_len = read_seq.size() - cigar_cnts.back();
      _clipLeft = endPos;
      int i = 0;
      while (align_len + i < read_seq.size() and ref_fa[_clipLeft - refBegin + i] == read_seq[align_len + i]) i++;
      if ( cigar_cnts.back() - i < CLIP_BP)
	continue;
      _clipLeft += i;
      pn_splitori_pos[pn].push_back( make_pair('L', _clipLeft) );

      //cout << "B " << _clipLeft << " " << i << " " << cigar << " " << estimated_clipPos << " " << beginPos << " " << endPos << " " << pn << endl;
      /*
      if ( i > 0) {      // move to right until can't move
	cout << infix(ref_fa, _clipLeft - refBegin - align_len, _clipLeft - refBegin) << endl;
	cout << read_seq.substr(0, align_len) << endl ;
	} */
    }        
  }
  fin.close();  
  
  map <int, int> clipLefts, clipRights;
  map <string, vector <pair<char, int> > >::iterator pi; 
  for (pi = pn_splitori_pos.begin(); pi != pn_splitori_pos.end(); pi++) {
    if (pn_has_clipPos(_clipLeft, _clipRight, pi->second, major_freq)) {
      addKey(clipLefts, _clipLeft);
      addKey(clipRights, _clipRight);
      cout << _clipLeft << " " << _clipRight << " " << pi->first << endl;
    }    
  }
  
  if ( major_key_freq(clipLefts, _clipLeft) < major_freq or
       major_key_freq(clipRights, _clipRight) < major_freq ) 
    return false;

  clipLefts.clear(); 
  clipRights.clear();
  pn_splitCnt.clear(); 
  clipLeft = _clipLeft;
  clipRight = _clipRight;
  for (pi = pn_splitori_pos.begin(); pi != pn_splitori_pos.end(); pi++) {
    vector <pair<char, int> >::iterator pii;
    for ( pii = (pi->second).begin(); pii != (pi->second).end(); pii++) {
      if ( ( (*pii).first == 'L' and abs((*pii).second - clipLeft) <= 3 ) or 
	   ( (*pii).first == 'R' and abs((*pii).second - clipRight) <= 3 ) )
	addKey(pn_splitCnt, pi->first);      
    }
  }  

  return pn_splitCnt.size() > 1; // ignore private mutations 
}

bool find_insert_pos(seqan::FaiIndex &faiIndex, unsigned fa_idx, string file_clip, int region_begin, int region_end, int flanking_len, float major_freq, int & clipLeft, int & clipRight, map <string,int> &pn_splitCnt) {
  pn_splitCnt.clear();
  stringstream ss;
  string line, pn, cigar, seq;
  int beginPos, endPos;
  ifstream fin;
  fin.open( file_clip.c_str());
  SplitPosInfo splitPosLeft, splitPosRight;
  while ( getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> pn >> beginPos >> endPos >> cigar;
    list <char> cigar_opts;
    list <int> cigar_cnts;
    parse_cigar(cigar, cigar_opts, cigar_cnts);
    /** one read can have both splitPosLeft and splitPosRight, due to error in cigar string. 
	eg. 556287 556309 S35M22S63 */    
    if ( *cigar_opts.begin() == 'S' and *cigar_cnts.begin() >= CLIP_BP) {
      if ( beginPos >= region_begin - flanking_len and beginPos < region_end + flanking_len) 
	addKey(splitPosRight, beginPos);
    }     
    if ( cigar_opts.back() == 'S' and cigar_cnts.back() >= CLIP_BP) {
      if ( endPos >= region_begin - flanking_len and endPos < region_end + flanking_len) 
	addKey(splitPosLeft, endPos);
        //cout << "debug begin " << beginPos << " " << cigar << " " << pn << endl;
    }        
  }
  fin.close();  
  int pos1, pos2;
  clipLeft = 0; clipRight = 0;
  if ( major_key_freq(splitPosLeft, pos1) >= major_freq or left_major_pos(splitPosLeft, CLIP_BP, pos1) >= major_freq )
    clipLeft = pos1;  
  if ( major_key_freq(splitPosRight, pos2) >= major_freq or right_major_pos(splitPosRight, CLIP_BP, pos2) >= major_freq )
    clipRight = pos2;
  //cout << "debug " << splitPosLeft.size() << " " << splitPosRight.size() << endl;  
  //cout << "debug2 " << clipLeft << " " << clipRight << endl;
  //print_map(splitPosRight);
  if (!clipLeft and !clipRight) return false;

  int estimated_clipPos = clipLeft ? clipLeft : clipRight;
  int ref_flanking = 200;
  int refBegin = estimated_clipPos - ref_flanking;
  int refEnd = estimated_clipPos + ref_flanking;
  seqan::CharString ref_fa = fasta_seq(faiIndex, fa_idx, refBegin, refEnd, true);
  return adjust_cigar(ref_fa, refBegin, file_clip, major_freq, clipLeft, clipRight, pn_splitCnt, estimated_clipPos, INSERT_POS_GAP);
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
  } else if ( opt == 2 ) { // for each insertion loci, combine reads from different pns 

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
    
    seqan::FaiIndex faiIndex;
    unsigned fa_idx = 0;    
    string file_fa_prefix = read_config(config_file, "file_fa_prefix");    
    string cnt, line;
    int region_begin, region_end;
    stringstream ss;
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      assert (!read(faiIndex, (file_fa_prefix + *ci + ".fa").c_str()) );
      assert (getIdByName(faiIndex, *ci, fa_idx));
      ifstream fin( (fin_pos + *ci + fout_suffix).c_str());
      assert (fin);
      ofstream fout( (path1 + "clip/" + *ci + fout_suffix).c_str() ); 
      while ( getline(fin, line) ) {
	ss.clear(); ss.str(line); 
	ss >> cnt >> region_begin >> region_end;
	string file_clip = fout_path + *ci + "_pos/" + int_to_string(region_begin) + "_" + int_to_string(region_end);
	int clipLeft, clipRight;
	map <string, int> pn_splitCnt;
	if (find_insert_pos(faiIndex, fa_idx, file_clip, region_begin, region_end, FLANK_REGION_LEN, MAJOR_SPLIT_POS_FREQ, clipLeft, clipRight, pn_splitCnt)) {
	  fout << region_begin << " " << region_end << " " << clipLeft << " " << clipRight;
	  for (map<string, int>::iterator pni = pn_splitCnt.begin(); pni != pn_splitCnt.end(); pni++)
	    fout << " " << pni->first << "," << pni->second;
	  fout << endl;
	}
      }
      fin.close();
      fout.close();
    }
  } else if ( opt == 0 ) {  // find exact split position
    string file_clip;
    bool insert_pos_found;
    seqan::FaiIndex faiIndex;
    unsigned fa_idx = 0;    
    string file_fa_prefix = read_config(config_file, "file_fa_prefix");    
    assert (!read(faiIndex, (file_fa_prefix + "chr1.fa").c_str()) );
    assert (getIdByName(faiIndex, "chr1", fa_idx));
    file_clip = fout_path + "chr1_pos/1574161_1574329";
    int region_begin = 1574161;
    int region_end = 1574329;
    int clipLeft, clipRight;
    map <string, int> pn_splitCnt;
    insert_pos_found = find_insert_pos(faiIndex, fa_idx, file_clip, region_begin, region_end, FLANK_REGION_LEN, MAJOR_SPLIT_POS_FREQ, clipLeft, clipRight, pn_splitCnt);
  }
  return 0;
}
