#define DEBUG_MODE  // test only chr1 for now
#define SEQAN_HAS_ZLIB 1
#include "utils.h"
#ifndef INSERT_UTILS_H
#define INSERT_UTILS_H

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

bool alu_mate_flag( BamFileHandler * bam_fh, MapFO &fileMap); 

string parse_alu_type(string alu_name);

#endif /*INSERT_UTILS_H*/
