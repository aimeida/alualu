#include "common.h"
#include "utils.h"

string read_config(string config_file, string key){
  string _key, _value;
  ifstream fin( config_file.c_str());
  if (!fin) 
    try {
      throw 1;    
    } catch(int e) {
      cerr << "#### ERROR #### file: "<< config_file << " not exists!!!" << endl;
    }   
  while (fin >> _key >> _value) {
    if (_key == key) 
      return _value;
  }
  fin.close();  
  try {
    throw 1;    
  } catch(int e) {
    cerr << "#### ERROR #### key: " << key <<" doesn't exist in file:" << config_file << endl;
    exit(0);
  }     
}

///// 1-based to 0-based.
////if (!seqan::lexicalCast2(beginPos, beginP) || beginPos <= 0) {
AluRefPos::AluRefPos(string file_alupos, int i) {
  flanking = i;
  ifstream fin( file_alupos.c_str());
  if (!fin) 
    try {
      throw 1;    
    } catch(int e) {
      cerr << "#### ERROR #### file: "<< file_alupos << " not exists!!!" << endl;
    }     
  string _tmp1, _tmp2, alu_type, chain;
  int bp, ep;
  while (fin >> _tmp1 >> bp >> ep >> alu_type >> _tmp2 >> chain){
    beginP.push(bp);
    endP.push(ep);
  }
  fin.close();
  cerr << "reading from " << file_alupos << " with " << beginP.size() << " loci\n" ;
}

int AluRefPos::updatePos(int &beginPos, int &endPos){
  beginPos = beginP.front() - flanking;
  endPos = endP.front() + flanking;
  beginP.pop();
  endP.pop();
  //cerr << "debug!! " << beginPos + flanking << endl;  
  return beginPos > 0 ? 1 : 0;
}  

AluRefPos::~AluRefPos(void){
}

void insertLen_of_nonUniq_mapping(vector<int> &starts_ends, vector<int> &reads_insert_len) {
  // starts_ends = [start_0, end_0, start_1, end_1, start_2, end_2 ....]
  vector <int> pre_ends;
  vector <int>::iterator vi, pi;
  int cur_begin, cur_end, insert_len;
  vi = starts_ends.begin();
  pre_ends.push_back(*(++vi));
  vi++;
  while ( vi != starts_ends.end() ) {
    cur_begin = *vi++;
    cur_end = *vi++; 
    for (pi = pre_ends.begin(); pi != pre_ends.end(); pi++) {
      insert_len = cur_begin - *pi;
      if (insert_len > 0) reads_insert_len.push_back(insert_len);
    }
    pre_ends.push_back(cur_end);
  }
  pre_ends.clear();    
}

