#ifndef UTILS_H
#define UTILS_H

string read_config(string config_file, string key);

class AluRefPos
{
  int flanking;
  queue<int> beginP, endP;
 public:
  AluRefPos(string file_alupos, int i);
  int updatePos(int &beginPos, int &endPos);
  ~AluRefPos(void);
};


void insertLen_of_nonUniq_mapping(vector<int> &starts_ends, vector<int> &reads_insert_len);

/*
// failed at compile !!
// template in header file only 
template <typename Type>
void printVec(vector<Type> &m){
  vector<Type>::iterator j;
  for (j=m.begin(); j!=m.end(); j++)  cerr << *j << "\t";
  cerr << endl;
};
*/

#endif /*UTILS_H*/
