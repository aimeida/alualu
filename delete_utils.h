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

bool get_align_pos(int aluBegin, int aluEnd, int beginPos, int endPos, int &ref_a, int &ref_b, int &read_a, int &read_b, seqan::BamAlignmentRecord &record);
bool split_global_align(seqan::CharString &fa_seq, seqan::BamAlignmentRecord &record, int read_a, int read_b, int &score, int &align_len);
T_READ classify_read(seqan::BamAlignmentRecord & record, int aluBegin, int aluEnd, FastaFileHandler *fasta_fh, bool debug = false);
void combine_pns_llh(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> &chrns, int col_00);
bool combine_pns_vcf(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> & chrns, int col_00);
void write_rm2(string f_input, string f_output, map < string, std::set<int> > & chrn_aluBegin, float chisq_th);
void filtered_vcf(string f_input, string f_output, int offset, map <string, set<int> > &chrn_aluBegin);
string tread_toStr(T_READ td);

#endif 
