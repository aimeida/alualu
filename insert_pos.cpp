#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

typedef seqan::Dna5String TSeq;
#define NNN_LEN 70
#define MAJOR_SPLIT 0.6
#define MAJOR_SPLIT_ONESIDE 0.8
#define SPLIT_OFFSET 3
#define MIN_SPLIT_READ 3

class SPLIT_INFO {
public:
  int split_to_refBegin; // dist to ref_begin
  bool alu_is_left;
  SPLIT_INFO( int p, bool ail) : split_to_refBegin(p), alu_is_left(ail) {}
};

template <class V>
struct compare_row {
  bool operator()(const V & a, const V & b) const {
    return (a.first == b.first) ? (a.second < b.second) :  (a.first < b.first);
  }
};

void parse_fa_name(string const &fa_name, int &ref_begin, int &ref_end, int &alu_len){ // eg.nl_95410890_95411310_AluY_299
  stringstream ss(fa_name);
  string token;
  getline(ss, token, '_');
  getline(ss, token, '_');
  seqan::lexicalCast2(ref_begin, token);
  getline(ss, token, '_');
  seqan::lexicalCast2(ref_end, token);
  getline(ss, token, '_');
  getline(ss, token, '_');
  seqan::lexicalCast2(alu_len, token);
}

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
    //fixme:: if ( operation == 'S') cerr << "soft ?? " << record.qName << endl;    
    if ( operation == 'N') {
      jump_len = count;
      break;
    }
  }
  return jump_len;
}

inline string get_sam_name(string &file_sam_prefix) {
  return file_sam_prefix + ".sam";
} 

inline string get_sam_name(string &file_sam_prefix, string suffix){
  return file_sam_prefix + "_" +  suffix + ".sam";
} 

void combine_sam(string file_sam_prefix){
  list <string> sam_suffix;
  list <string>::iterator si;
  sam_suffix.push_back("fl");
  sam_suffix.push_back("fr");
  sam_suffix.push_back("rl");
  sam_suffix.push_back("rr");
  si = sam_suffix.begin();
  while (si != sam_suffix.end()){
    int ie = is_nonempty_file( get_sam_name(file_sam_prefix, *si).c_str());
    if ( ie == 0 ) system ( ("rm "+ get_sam_name(file_sam_prefix, *si) ).c_str() );    // exist, but empty
    if ( ie > 0 ) si++;
    else si = sam_suffix.erase(si);  // not exist
  }
  if (!sam_suffix.empty()) { 
    si = sam_suffix.begin();    
    system( ("mv " + get_sam_name(file_sam_prefix, *si++) + " " +  get_sam_name(file_sam_prefix)).c_str() );
    while ( si != sam_suffix.end() ) {
      system( ("cat " + get_sam_name(file_sam_prefix, *si) + " >> " +  get_sam_name(file_sam_prefix)).c_str() );
      system( ("rm " + get_sam_name(file_sam_prefix, *si) ).c_str() );
      si++;
    }
  } 
}

// move split pos further to ref_begin
void moveSplit_toRefBegin(int  & split_begin, TSeq const & ref, int nnn_begin, int jump_len){  
  while (split_begin < nnn_begin and split_begin + jump_len < (int)length(ref) and ref[split_begin] == ref[ split_begin + jump_len] ) 
    split_begin++;
} 

void read_sam(string fin_sam, string fin_ref_prefix, vector< SPLIT_INFO *> & vec_splitPos, int & ref_begin) {
  string line, qName, col2, fa_name;
  stringstream ss;
  vector <string> fa_names;
  ifstream fin(fin_sam.c_str());
  assert(fin);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str(line);
    ss >> qName >> col2 >> fa_name;
    fa_names.push_back(fa_name);
  }  
  fin.close();
  
  int split_begin, jump_len, ref_end, alu_len;  // relative to ref (hg18 + NNN + Alu) 
  parse_fa_name(fa_name, ref_begin, ref_end, alu_len);
  int nnn_begin = 0, nnn_end = 0;
  int ref_len_aluN = alu_len + NNN_LEN;
  vector <string>::iterator fi = fa_names.begin();
  map<string, TSeq> name_to_ref;
  seqan::BamStream bamStream(fin_sam.c_str());  
  seqan::BamAlignmentRecord record;
  while (!atEnd(bamStream)){
    readRecord(record, bamStream);
    fa_name = *fi++;
    bool alu_is_left = (fa_name[1] == 'l');
    nnn_begin = alu_is_left ? alu_len : (ref_end - ref_begin);
    nnn_end = nnn_begin + NNN_LEN;
    if ( !findSplit(split_begin, jump_len, record)) continue;  // split_end = split_begin + jump_len
    if ( split_begin > nnn_begin || split_begin + jump_len < nnn_end) continue;

    if ( name_to_ref.find(fa_name) == name_to_ref.end() ) {
      string fin_ref = fin_ref_prefix + "_" + fa_name[1] + ".fa";
      read_fasta_by_name(name_to_ref[fa_name], fin_ref, fa_name);
    }
    moveSplit_toRefBegin(split_begin, name_to_ref[fa_name], nnn_begin, jump_len);        
    if ( alu_is_left )  vec_splitPos.push_back( new SPLIT_INFO(split_begin + jump_len - ref_len_aluN, true) ) ;
    else  vec_splitPos.push_back( new SPLIT_INFO(split_begin, false) );     
  }
  seqan::close(bamStream);
}

float majorKey(map <int, int> &m, int & key)  
{
  int max_count = 0, sum_count = 0;
  if ( m.empty() ) return 0;
  key = m.begin()->first;
  for (map <int, int>::iterator it = m.begin(); it != m.end(); it++) { 
    sum_count +=  it->second;
    if ( it->second > max_count ) {
      max_count = it->second;
      key = it->first;
    }
  }
  return max_count/(float) sum_count;
};

void major_pos( map <int, int > & map_pos, float th, int & main_pos ) {  
  main_pos = 0;
  int mp;
  float main_th = majorKey( map_pos, mp);
  if ( main_th >= th ) main_pos = mp;
  else if ( main_th > 0 ) {
    int max_count = 0, sum_count = 0;
    for (map <int, int >::iterator it = map_pos.begin(); it != map_pos.end(); it++) {
      int k = it->first;
      int v = it->second;
      sum_count +=  it->second;
      if ( abs(k - mp) <= SPLIT_OFFSET ) max_count += v;
    }
    if ( max_count/(float) sum_count >= th )  main_pos = mp;
  }  
}

bool major_pos_exist( vector< SPLIT_INFO *> & vec_splitPos, int & split_pa, int & split_pb){
  map < int, int > pos_left, pos_right;
  vector < SPLIT_INFO *>::iterator si;
  for ( si = vec_splitPos.begin(); si != vec_splitPos.end(); si++ )
    if ( (*si)->alu_is_left ) addKey(pos_left, (*si)->split_to_refBegin ) ;
    else  addKey (pos_right,  (*si)->split_to_refBegin ) ;
  int pa, pb;
  major_pos(pos_left, MAJOR_SPLIT_ONESIDE, pa);
  major_pos(pos_right, MAJOR_SPLIT_ONESIDE, pb);
  if ( pa || pb) {
    split_pb = max(pa, pb);
    split_pa = min(pa, pb) ? min(pa, pb) :split_pb;
    return true;
  } 
  for (map < int, int >::iterator pi = pos_right.begin(); pi != pos_right.end(); pi++) 
    addKey(pos_left, pi->first, pi->second);
  major_pos(pos_left, MAJOR_SPLIT, pa);
  if (pa) {
    split_pa = pa;
    split_pb = pa;
    return true;
  }
  return false;
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  int opt;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];

  boost::timer clocki;    
  clocki.restart();
  
  string path1 = read_config(config_file, "file_alu_insert1") ;    
  string path_sam = path1 + "split_mapping_splazer/";  
  string path_map = path1 + "split_mapping_clip/";
  string path_ref = path1 + "split_mapping_fastq/";

  map<int, string> ID_pn;
  get_pn(read_config(config_file, "file_pn"), ID_pn);  
  map <string, set< pair<int, int> > > pos_of_clippedReads;    
  map<string, set< pair<int, int> > >::iterator pi;
  ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
  int idx_pn;
  while (fin >> idx_pn) 
    get_filename_fastq(path_map + ID_pn[idx_pn] + ".map", pos_of_clippedReads); // only to update pos_of_clippedReads
  fin.close();    
  
  if (opt == 1) { // try diff functions 
    
    

  } else if (opt == 2) { // combine sam    
    for ( pi = pos_of_clippedReads.begin(); pi != pos_of_clippedReads.end(); pi++ ) {
      for (set< pair<int, int> >::iterator ri = (pi->second).begin(); ri != (pi->second).end(); ri++) {	
	string file_sam_prefix = path_sam + pi->first + "_" + int_to_string( (*ri).first ) + "_" + int_to_string ((*ri).second);
	combine_sam(file_sam_prefix); 
      }
    }    
    cout << "combine *sam files for each loci, done\n";

  } else if (opt == 3) { // write pos file for split position 
    string fin_sam, fin_ref_prefix;
    set < pair<int, int>, compare_row < pair<int, int> > > set_insertPos;
    vector< SPLIT_INFO *> vec_splitPos;
    for ( pi = pos_of_clippedReads.begin(); pi != pos_of_clippedReads.end(); pi++ ) {
      string chrn = pi->first;
      int nf_all = 0, nf_pos = 0;
      set_insertPos.clear();
      for (set< pair<int, int> >::iterator ri = (pi->second).begin(); ri != (pi->second).end(); ri++) {	
	fin_sam = path_sam + chrn + "_" + int_to_string( (*ri).first ) + "_" + int_to_string ((*ri).second) + ".sam";
	if (is_nonempty_file(fin_sam) <= 0) continue;
	nf_all ++ ;
	clear_vec_of_ptr(vec_splitPos);
	fin_ref_prefix = path_ref + chrn + "_" + int_to_string( (*ri).first ) + "_" + int_to_string ((*ri).second);
	int split_pa, split_pb, ref_begin;
	read_sam(fin_sam, fin_ref_prefix, vec_splitPos, ref_begin);	
	if (vec_splitPos.size() >= MIN_SPLIT_READ and major_pos_exist(vec_splitPos, split_pa, split_pb) ) {
	  nf_pos ++;
	  cout << ref_begin << " " << ref_begin + split_pa << " " << ref_begin + split_pb << endl;
	  set_insertPos.insert( make_pair(ref_begin + split_pa, ref_begin + split_pb) );	  
	}	
      }
      
      set < pair<int, int>, compare_row < pair<int, int> > >::iterator si;
      ofstream fout( (path_sam + chrn + ".pos").c_str() );
      for (si = set_insertPos.begin(); si != set_insertPos.end(); si++)
	fout << (*si).first << " " << (*si).second << endl;
      fout.close();
      
      cout << nf_pos << " out of " <<  nf_all << " regions have unambiguous split position.\n";
      cout << "written into: " << path_sam << chrn << ".pos" << endl;
    }
  }   
}
