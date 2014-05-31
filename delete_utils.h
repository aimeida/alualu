// read distribution of insert length, etc
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/basic.h>
#include "common.h"
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
inline void genoProb_print(float *gp, ostream & os, int precision){ os << " " << setprecision(precision) << gp[0] << " " << gp[1] << " " << gp[2]; }
inline bool p00_is_dominant(float * log10_p, int min_log10p) { return  max( log10_p[2] - log10_p[0], log10_p[1] - log10_p[0]) <= min_log10p; }

inline string phred_log (float p) { return p ? (int_to_string (-(int)(log10 (p) * 10)) ) : "255"; }

string phred_scaled(float p0, float p1, float p2);
bool get_align_pos(int aluBegin, int aluEnd, int beginPos, int endPos, int &ref_a, int &ref_b, int &read_a, int &read_b, seqan::BamAlignmentRecord &record);
bool split_global_align(seqan::CharString &fa_seq, seqan::BamAlignmentRecord &record, int read_a, int read_b, int &score, int &align_len);
T_READ classify_read(seqan::BamAlignmentRecord & record, int align_len, string chrn, int aluBegin, int aluEnd, FastaFileHandler *fasta_fh, bool read_is_left);
void log10P_to_P(float *log_gp, float *gp, int max_logp_dif);

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
