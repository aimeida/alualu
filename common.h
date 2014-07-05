// for header files and debug print 
#ifndef COMMON_H
#define COMMON_H

#include <map>
#include <set>
#include <list>
#include <queue>
#include <vector>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <math.h> 
#include <dirent.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include <boost/timer.hpp>

using namespace std;

// macro for alu_delete
#define LOG10_RATIO_UB 3 // max log(OR) for a single read, estimated from quantile(0.95)/quantile(0.05), 1e3 to 1e4
#define LOG10_GENO_PROB 5 // min genotype prob by reads 
#define ALU_FLANK 600
#define BOUNDARY_OFFSET 10 
#define CLIP_BP 10  // min length to be called soft clip 
#define COVERAGE_HIGH 90
// macro for alu_insert
#define DISCORDANT_LEN 2000 
#define INS_COVERAGE_MAX 6 // a random number, remove some high coverage regions 
#define SCAN_WIN_LEN 400   // 0.99 quantile. 450 for 0.999 quantile. only approximate,different reading group differs
#define LEFT_PLUS_RIGHT 3  // minimum sum of left and right reads
#define MAX_LEN_REGION 100 // max length of insert region
#define MAX_POS_DIF 50  // combine scanned insertion loci 
#define ALIGN_END_CUT 10 // cut the last 10 bp while realign to reference
#define CLIP_BP_LEFT 5
#define CLIP_BP_RIGHT 30
#define CLIP_Q 30 

inline int ceil_by_resolution(int x, int r) {
  return r * ceil ( (float) x / (float) r );
}

typedef pair<int, string > IntString;
inline bool compare_IntString(const IntString & a, const IntString & b) {
  return (a.first == b.first) ? (a.second < b.second) :  (a.first < b.first);
}

inline void check_folder_exists(string path) {
  if ( access( path.c_str(), 0 ) != 0 ) system( ("mkdir " + path).c_str() );    
}

inline void move_files(string path_move, string fns) {
  check_folder_exists(path_move);
  system(("mv " + fns + " " + path_move).c_str());    
}

inline int get_col_idx(string fn, string col_name) {
  int coln = 0;
  ifstream fin( fn.c_str());
  assert(fin); 
  string header, tmp1;
  getline(fin, header);
  stringstream ss;
  ss.clear(); ss.str( header );
  while ( ss >> tmp1) { 
    coln++;
    if (tmp1 == col_name) {
      fin.close();
      return coln;
    }
  }
  fin.close();
  return 0;
}

template < typename K >
bool compare_first(const K & a, const K & b) {
  return (a.first == b.first) ? (a.second < b.second) :  (a.first < b.first);
}

template<typename A, typename B>
  std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
  return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
  std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
  std::multimap<B,A> dst;
  std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
		 flip_pair<A,B>);
  return dst;
}

template <class K, class V> 
void addKey(map <K,V> &m, K key, int cnt=1)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) == m.end()) {
    m[key] = cnt;
  } else (it->second) += cnt;
};

template <class K, class V> 
bool key_exists(map <K,V> &m, K key)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) != m.end()) return true;
  return false;
};

template <class K, class V> 
void get_mapVal(map <K,V> &m, K key, V & default_val)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) != m.end()) 
    default_val = it->second;
};

template <class K, class V>
void debugprint_map(map <K,V> &m, size_t n_max = 0) {
  cout << "size of map " << m.size() << endl;
  typename map <K,V>::iterator it;
  if (n_max == 0) n_max = m.size();
  size_t ni = 0;
  for (it = m.begin(); it != m.end(),ni < n_max; it++, ni++) cerr << it->first << ":" << it->second << " ";
  cerr << endl;  
};

template <class K>
void debugprint_vec(vector <K> &m) {
  typename vector <K>::iterator it;
  for (it = m.begin(); it != m.end(); it++) cerr << *it << " ";
  cerr << endl;
};

template <typename VT>
void sort_file_by_col(string fn, int coln, bool has_header){
  assert( coln > 0 );
  ifstream fin( fn.c_str());
  assert(fin);
  string line, tmpv, header;
  stringstream ss;
  VT valn;
//  typedef std::list< pair<VT, std::string> > LT;
  std::list< pair<VT, std::string> > rows;
  if (has_header) getline(fin, header);
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    if (coln == 1) ss >> valn;
    else {
      for (int i = 0; i < coln-1; i++) ss >> tmpv;
      ss >> valn;
    }
    rows.push_back( make_pair(valn, line) );
  }
  fin.close();
  rows.sort(compare_first < pair<VT, string> >);
  ofstream fout( fn.c_str());
  if (has_header) fout << header << endl;  
  // tell it it's typename !!!!
  for (typename std::list< pair<VT, std::string> >::iterator ri = rows.begin(); ri != rows.end(); ri++)
    fout << (*ri).second << endl;
  fout.close();
};

class ConfigFileHandler{
 public:
  map <string, string> configs;
  ConfigFileHandler(string config_file);
  string get_conf(string key);
};

string int_to_string(int i);
string get_pn(string pn_file, int idx_pn);
void get_pn(string pn_file, map<int, string> &ID_pn);
void read_file_pn_used(string fn, std::set <string> & pns_used);
void read_file_pn_used(string fn, vector <string> & pns_used);
int check_file_size(string fn);

#endif /*COMMON_H*/
