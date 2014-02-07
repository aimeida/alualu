// read parameters used in the program 
#include "common.h"

template <typename Type>
Type read_config(string config_file, string key){
  string _key;
  Type _value;
  ifstream fin( config_file.c_str());
  if (!fin) 
    try {
      throw 1;    
    } catch(int e) {
    cout << "config_file not exists!!!" << endl;
  }   
  while (fin >> _key >> _value) 
    if (_key == key) return _value;
  fin.close();
  
  return (Type) 0;
}
