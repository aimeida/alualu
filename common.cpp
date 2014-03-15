#include "common.h"

string int_to_string(int i) {
  stringstream ss;
  ss << i;
  return ss.str();
}

string get_pn(string pn_file, int idx_pn){
  ifstream fin( pn_file.c_str());
  string pn;
  int i = 0;
  while (fin >> pn) {
    if (i == idx_pn ) {
      fin.close();  
      return pn;
    }
    i++;
  }
  return "na";
}

void print_vec(vector<int> &m){  
  for (vector<int>::iterator j=m.begin(); j!=m.end(); j++)  cerr << *j << "\t";
  cerr << endl;
}

  //for (map<int, float>::iterator it=prob_vec.begin(); it!=prob_vec.end(), i < 10; ++it, i++)
  //cout << it->first << " => " << it->second << '\n';
