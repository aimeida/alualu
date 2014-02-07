// read distribution of insert length, etc
#include "common.h"

class ProbByQuantile
{
  vector<int> prob_vec; 
 public:
  ProbByQuantile(string quantile_file){ 
    cout << "build dist from file:" << quantile_file << endl;
    prob_vec.push_back(1);
  }
  float prob_inputlen(int inputlen) {
    if (inputlen > 1000) return 1.;
    return 0.8;
  }  
};
