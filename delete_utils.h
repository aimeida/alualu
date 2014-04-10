// read distribution of insert length, etc
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/basic.h>
#include "common.h"
#include "utils.h"
#ifndef DELETE_UTILS_H
#define DELETE_UTILS_H

#define DEFAULT_READ_LEN 100 // if unknown, use default
#define LOG_RATIO_UB 7 // max log(OR) for a single read, estimated from quantile(0.95)/quantile(0.05), 1e3 to 1e4
#define ALU_FLANK 600
#define BOUNDARY_OFFSET 10 
#define CLIP_BP 10  // min length to be called soft clip 

enum T_READ {useless_read, unknow_read, mid_read, clip_read};
typedef map<int, int> MapII;
typedef map<int, int>::iterator MapIIt;

inline int min_align_score(int align_len) { return round(0.75 * align_len - 2);}
bool get_align_pos(int aluBegin, int aluEnd, int beginPos, int endPos, int &ref_a, int &ref_b, int &read_a, int &read_b, seqan::BamAlignmentRecord &record);
bool split_global_align(seqan::CharString &fa_seq, seqan::BamAlignmentRecord &record, int read_a, int read_b, int &score, int &align_len);
string enum_to_string(T_READ x);
T_READ classify_read(seqan::BamAlignmentRecord & record, int align_len, int aluBegin, int aluEnd, seqan::FaiIndex &faiIndex, unsigned fa_idx, bool read_is_left);
bool check_delete_region(string const &bam_input, string const &bai_input, string const & fa_input, string chrn, int beginPos, int endPos );
bool normalize_prob(float *log_p);
 
class EmpiricalPdf
{
  map<int, float> prob_vec; 
  int min_len, max_len, bin_width;
  float min_prob;
 public:
  EmpiricalPdf(string pdf_file);
  float pdf_obs(int insertlen);
};
void genotype_prob(map <string, vector<int> >  &insertlen_rg, map <string, EmpiricalPdf *> &empiricalpdf_rg, int alu_len, float *log_p);


#endif 
