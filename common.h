// for header files and debug print 
#ifndef COMMON_H
#define COMMON_H

#include <map>
#include <set>
#include <queue>
#include <vector>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <math.h> 
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include <boost/timer.hpp>

using namespace std;

// macro for alu_delete
#define DEFAULT_READ_LEN 100 // if unknown, use default
#define LOG10_RATIO_UB 3 // max log(OR) for a single read, estimated from quantile(0.95)/quantile(0.05), 1e3 to 1e4
#define LOG10_GENO_PROB 5 // min genotype prob by reads 
#define ALU_FLANK 600
#define BOUNDARY_OFFSET 10 
#define CLIP_BP 10  // min length to be called soft clip 
#define COVERAGE_HIGH 90
// macro for alu_insert
#define ALU_MIN_LEN 50  // min overlap in alu region
#define INSERT_POS_GAP 40  // (seq, ref/alu_insertion, seq), max length of ref ( replaced by alu_insert)
#define DISCORDANT_LEN 2000 
#define INS_COVERAGE_MAX 6 // a random number, remove some high coverage regions 
#define SCAN_WIN_LEN 400   // 0.99 quantile. 450 for 0.999 quantile. only approximate,different reading group differs
#define LEFT_PLUS_RIGHT 4  // minimum sum of left and right reads
#define MAX_LEN_REGION 100 // max length of insert region
#define MAX_POS_DIF 50  // combine scanned insertion loci 
#define NUM_PN_JOIN 10  // max number of pns to sort and join
// macros for insert_pos
#define MAJOR_SPLIT_POS_FREQ 0.6
#define FLANK_REGION_LEN 80
#define ALIGN_END_CUT 10 // cut the last 10 bp while realign to reference
// macros for ins_del
#define REF_EXT_LEN 150 // extend exact insert pos to both sides. This data, max read len = 135
#define ALU_INSERT_FLANK 400  // 600 (alu_delete.cpp) - 200 (min alu length)
#define ALUCONS_SIMILARITY 0.8
#define MIN_READS_CNT 3 // min interesting reads, in order to consider this pos
#define MIN_MATCH_LEN 60 // in order to use this read

string int_to_string(int i);
string get_pn(string pn_file, int idx_pn);
void get_pn(string pn_file, map<int, string> &ID_pn);
int is_nonempty_file(string fn);

class ConfigFileHandler{
 public:
  map <string, string> configs;
  ConfigFileHandler(string config_file);
  string get_conf(string key);
};

template <class K, class V> 
void addKey(map <K,V> &m, K key, int cnt=1)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) == m.end()) {
    m[key] = cnt;
  } else (it->second) += cnt;
};

template <class K, class V>
void print_map(map <K,V> &m, size_t n_max = 0) {
  cout << "size of map " << m.size() << endl;
  typename map <K,V>::iterator it;
  if (n_max == 0) n_max = m.size();
  size_t ni = 0;
  for (it = m.begin(); it != m.end(),ni < n_max; it++, ni++) cerr << it->first << ":" << it->second << " ";
  cerr << endl;  
};

template <class K>
void print_vec(vector <K> &m) {
  typename vector <K>::iterator it;
  for (it = m.begin(); it != m.end(); it++) cerr << *it << " ";
  cerr << endl;
};


#endif /*COMMON_H*/
