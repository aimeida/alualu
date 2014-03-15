// read distribution of insert length, etc
#include "common.h"

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
bool normalize_prob(float *log_p);
