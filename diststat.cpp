
#include "diststat.h"

EmpiricalPdf::EmpiricalPdf(string pdf_file){ 
  int pos;
  float posp;
  bin_width = 0;
  ifstream fin( pdf_file.c_str());
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
  // cerr << "read pdf dist " << min_len << " " << max_len << " " << bin_width << " " << min_prob << endl;
}

float EmpiricalPdf::pdf_obs(int insertlen) {
  if (insertlen >= max_len || insertlen <= min_len) return min_prob;
  int nearby_pos = min_len + (insertlen-min_len)/bin_width*bin_width;
  map<int, float >::iterator it = prob_vec.find(nearby_pos);
  if ( it != prob_vec.end())  return it->second;
  return min_prob;
}   

void genotype_prob(vector<int> &reads_insert_len, ReadsPosStore &reads_pos, int alu_len, float *log_p){
  ReadsPosStore::iterator rp;
  vector<int>::iterator ri; 
  int insertlen, left_read_end;
  for (rp = reads_pos.begin(); rp != reads_pos.end(); rp++){
    if (rp->second.size() !=4 ) continue; // both ends mapped to this region && unique mapping
    ri = rp->second.begin();
    left_read_end = *(++ri);
    insertlen = *(++ri) - left_read_end;
    cerr << "debug " << insertlen <<  endl;
    print_vec(rp->second);
    //    for (vector<int>::iterator pi = reads_insert_len.begin(); pi != reads_insert_len.end(); pi++){
    //      float p_y = empiricalpdf->pdf_obs(*pi);
//      float p_z = empiricalpdf->pdf_obs(alu_len + *pi);
//      log_p[0] += log(p_y);
//      log_p[1] += log(0.5 * p_y + 0.5 * p_z);
//      log_p[2] += log(p_z); 
//    }
//    post_prob(log_p, 3);
  }
}

void normalize_prob(float *log_p, int len){
  float *p = log_p;
  float p_sum = 0;
  int i = 0;
  for (p = log_p, i = 0; i < len; i++, p++)  cerr << "start " << *p << endl;

  for (p = log_p, i = 0; i < len; i++, p++) {
    *p = exp(*p);
    p_sum += *p;
  }
  
  cerr << "need ratio trick for calculation !! \n";
  for (p = log_p, i = 0; i < len; i++, p++) *p /= p_sum;
  for (p = log_p, i = 0; i < len; i++, p++)  cerr << "print " << *p << endl;
}
