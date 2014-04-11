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
#define LOG10_RATIO_UB 5 // max log(OR) for a single read, estimated from quantile(0.95)/quantile(0.05), 1e3 to 1e4
#define ALU_FLANK 600
#define BOUNDARY_OFFSET 10 
#define CLIP_BP 10  // min length to be called soft clip 
#define COVERAGE_HIGH 90

enum T_READ {useless_read, unknow_read, mid_read, clip_read};

inline int min_align_score(int align_len) { return round(0.75 * align_len - 2);}
inline void genoProb_print(float *gp, ostream & os, int precision){ os << " " << setprecision(precision) << gp[0] << " " << gp[1] << " " << gp[2]; }
inline bool p00_is_dominant(float * log10_p, int min_log10p) { return  max( log10_p[2] - log10_p[0], log10_p[1] - log10_p[0]) <= min_log10p; }

inline string phred_log (float p) { return p ? (int_to_string (-(int)(log10 (p) * 10)) ) : "255"; }
inline string phred_scaled(float p0, float p1, float p2){ return phred_log(p0)  + "," + phred_log(p1) + "," + phred_log(p2); }

bool get_align_pos(int aluBegin, int aluEnd, int beginPos, int endPos, int &ref_a, int &ref_b, int &read_a, int &read_b, seqan::BamAlignmentRecord &record);
bool split_global_align(seqan::CharString &fa_seq, seqan::BamAlignmentRecord &record, int read_a, int read_b, int &score, int &align_len);
string enum_to_string(T_READ x);
T_READ classify_read(seqan::BamAlignmentRecord & record, int align_len, int aluBegin, int aluEnd, seqan::FaiIndex &faiIndex, unsigned fa_idx, bool read_is_left);
void log10P_to_P(float *log_gp, float *gp, int max_logp_dif);
bool check_delete_region(string const &bam_input, string const &bai_input, string const & fa_input, string chrn, int beginPos, int endPos );

class EmpiricalPdf
{
  map<int, float> prob_vec; 
  int min_len, max_len, bin_width;
  float min_prob;
 public:
  EmpiricalPdf(string pdf_file);
  float pdf_obs(int insertlen);
};


#endif 
