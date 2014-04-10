#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "common.h"
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
};

bool readbam_clip(seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, size_t ref_begin, int ref_end, seqan::BamStream &bamStreamOut, vector<seqan::CharString> & qnames, size_t offset);
void write_from_bam(string path_fout, string fin_read_bam, string fin_read_map, map <string, set< pair<int, int> > > & chr_clipReads, string format);

void write_fastq_seqcons(string path_fout, string fin_read_bam, string fin_read_map, map <string, set< pair<int, int> > > & chr_clipReads, string format);
void write_fastq_seqcons_nopair(string path_fout, string fin_read_bam, string fin_read_map, map <string, set< pair<int, int> > > & chr_clipReads, string format);

#endif /*INSERT_UTILS_H*/
