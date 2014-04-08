#include "diststat.h"

EmpiricalPdf::EmpiricalPdf(string pdf_file){ 
  int pos;
  float posp;
  bin_width = 0;
  ifstream fin( pdf_file.c_str());
  assert (fin);
  fin >> pos >> posp;
  prob_vec[pos] = posp;
  min_len = pos;
  min_prob = posp;
  while (fin >> pos >> posp) {
    prob_vec[pos] = posp;
    min_prob = min_prob > posp ? posp : min_prob;      
    if (!bin_width) bin_width = pos - min_len;
  }
  max_len = pos;
  fin.close();
  cerr << "read pdf dist " << pdf_file << endl;
  cerr << min_len << " " << max_len << " " << bin_width << " " << min_prob << endl;
}

float EmpiricalPdf::pdf_obs(int insertlen) {
  if (insertlen >= max_len || insertlen <= min_len) return min_prob;
  int nearby_pos = min_len + (insertlen-min_len)/bin_width*bin_width;
  map<int, float >::iterator it = prob_vec.find(nearby_pos);
  if ( it != prob_vec.end())  return it->second;
  return min_prob;
}   

void add_read_prob(){
  
}

void genotype_prob(map < string, vector<int> >  &insertlen_rg, map <string, EmpiricalPdf *> &empiricalpdf_rg, int alu_len, float *log_p){
  for (int i = 0; i < 3; i++) log_p[i] = 0;
  EmpiricalPdf *empiricalpdf;
  for (map < string, vector<int> >::iterator irg=insertlen_rg.begin(); irg!=insertlen_rg.end(); irg++) {
    empiricalpdf = empiricalpdf_rg[irg->first];
    for (vector<int>::iterator ir = irg->second.begin(); ir != irg->second.end(); ir++ ) {
      float p_y = empiricalpdf->pdf_obs(*ir);
      float p_z = empiricalpdf->pdf_obs(alu_len + *ir);      
      log_p[0] += log(p_y);
      log_p[1] += log(0.67 * p_y + 0.33 * p_z);
      log_p[2] += log(p_z); 
    }
  }
  normalize_prob(log_p);
}

bool normalize_prob(float *log_p){
  float r20 = exp(log_p[2] - log_p[0]); 
  float r10 = exp(log_p[1] - log_p[0]); 
  if ( max(r20, r10) < 1e-5 ) return false;
  float r_sum = r20 + r10 + 1;
  log_p[0] = 1./r_sum;
  log_p[1] = r10/r_sum;
  log_p[2] = r20/r_sum;
  return true;
}



// int main( int argc, char* argv[] )
// {
//   cerr << "test code here \n";
//   return 0;
// }
