#define DEBUG_MODE  // test only chr1 for now
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "common.h"
#include "utils.h"
#ifndef INSERT_UTILS_H
#define INSERT_UTILS_H

#define ALU_MIN_LEN 50  // min overlap in alu region
#define DEFAULT_READ_LEN 100 // if read length unknown, use this default
#define CLIP_BP 10 
#define DISCORDANT_LEN 2000 
#define INS_COVERAGE_MAX 6 // a random number, remove some high coverage regions 
#define SCAN_WIN_LEN 400   // 0.99 quantile. 450 for 0.999 quantile. However different reading group differs, only approximate
#define LEFT_PLUS_RIGHT 4  // minimum sum of left and right reads
#define MAX_LEN_REGION 100 // max length of insert region
#define MAX_POS_DIF 50  // combine scanned insertion loci 
#define NUM_PN_JOIN 10  // max number of pns to sort and join

// to use, but not in use yet
#define READ_CUT_BP 5      // cut 5 bp in the end of the reads
#define MAX_READS_CNT_CLIP 50 // max. number of reads to infer insert position, don't need more than that
#define MAX_READS_CNT_CONS 50 // max. number of reads to construct consensus region 
// following #define not used any more 
#define FLAG_LONG -1
#define FLAG_RC -2

typedef map<int, ofstream* > MapFO;
typedef map<string, ofstream* > MapSFO;
typedef pair<int, string > RowInfo;
typedef seqan::Dna5String TSeq;

class READ_INFO {
public:
  int beginPos, endPos;
  string alu_type;
  READ_INFO(int p, int lr, string alu_type) : beginPos(p), endPos(p+lr-2), alu_type(alu_type) {}
};

struct compare_row {
  bool operator()(const RowInfo& a, const RowInfo& b) const {
    return (a.first == b.first) ? (a.second < b.second) :  (a.first < b.first);
  }
};

inline bool compare_list(const RowInfo& a, const RowInfo& b) {
  return (a.first == b.first) ? (a.second < b.second) :  (a.first < b.first);
}

inline void refPos_by_region(int region_begin, int region_end, int & ref_begin, int & ref_end){
  ref_begin = region_begin - 2 *  DEFAULT_READ_LEN;
  ref_end = region_end + 2 * DEFAULT_READ_LEN;
}

inline void searchPos_by_region(int region_begin, int region_end, int & search_begin, int & search_end){
  search_begin = region_begin - SCAN_WIN_LEN;
  search_end = region_end + SCAN_WIN_LEN;
}

inline void regionPos_by_ref(int &region_begin, int &region_end, int ref_begin, int ref_end){
  region_begin = ref_begin + DEFAULT_READ_LEN;
  region_end = ref_end - 2 * DEFAULT_READ_LEN;
}

bool alu_mate_flag( string bam_input, map<int, seqan::CharString> const &rID_chrn, MapFO &fileMap); 
bool alu_mate_flag_depreciate( string bam_input, map<int, seqan::CharString> const &rID_chrn, MapFO &filoeMap); 

// to develop or not in use any more 
string parse_alu_type(string alu_name);
void get_alu_type(string fn, map<int, string> &pos_aluType);

#endif /*INSERT_UTILS_H*/
