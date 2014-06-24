// read distribution of insert length, etc
#define SEQAN_HAS_ZLIB 1
#include "utils.h"
#ifndef DELETE_UTILS_H
#define DELETE_UTILS_H

enum T_READ {useless_read, unknow_read, mid_read, clip_read}; // used for deletions

class GENO_PROB {
 public:
  float g0, g1, g2;
 GENO_PROB(float g0, float g1, float g2) : g0(g0), g1(g1), g2(g2) {}
};

inline int min_align_score(int align_len) { return round(0.75 * align_len - 2);}
bool get_align_pos(int aluBegin, int aluEnd, int beginPos, int endPos, int &ref_a, int &ref_b, int &read_a, int &read_b, seqan::BamAlignmentRecord &record);
bool split_global_align(seqan::CharString &fa_seq, seqan::BamAlignmentRecord &record, int read_a, int read_b, int &score, int &align_len);
T_READ classify_read(seqan::BamAlignmentRecord & record, int aluBegin, int aluEnd, FastaFileHandler *fasta_fh, bool only_tLen_info = false);
void filter_by_llh_noPrivate(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> &chrns, int col_00);
bool combine_pns_vcf_noPrivate(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> & chrns, int col_00);
string debug_print_tread(T_READ td);


#endif 
