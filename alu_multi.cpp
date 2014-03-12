#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <sstream>
#include "common.h"
#include "utils.h"

void combine_multi(vector <string> &pns, string &file_pn_delete_prefix, ofstream &fout, string chrx){
  string _chrx, tmp1, tmp2;
  float p0, p1, p2, coverage;
  int aluBegin;
  map <int, int> del_count;
  map <int, int>::iterator dc;

  ifstream fin;
  for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
    fin.open( (file_pn_delete_prefix + *pi).c_str() );
    assert(fin);
    while (fin >> _chrx >> tmp1 >> aluBegin >> tmp2 >> p0 >> p1 >> p2 >> coverage) {
      if (chrx != _chrx) continue;
      if (p1 > p0 or p2 > p0) addKey(del_count, aluBegin, p1 > p2 ? 1 : 2); // only consider potential deletion region 
    }
    fin.close();
  }
  //print_map(del_count);
  float pn_chr = 2. * pns.size();
  float fp, f0, f1, f2;
  map <int, float> del_llh;
  map <int, float>::iterator dl;
  map <int, vector <float> > del_debug;
  for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
    fin.open( (file_pn_delete_prefix + *pi).c_str() );
    while (fin >> _chrx >> tmp1 >> aluBegin >> tmp2 >> p0 >> p1 >> p2 >> coverage) {
      if (_chrx != chrx) continue;
     fp = del_count[aluBegin] /  pn_chr;      
      /* assume HW */
      f0 = (1 - fp) * (1 - fp);
      f1 = 2 * fp * (1 - fp);
      f2 = fp * fp;
      float tmp_r = (f0 * p0 + f1 * p1 + f2 * p2) / p0;
      if ( (dl = del_llh.find(aluBegin)) == del_llh.end() ) {
	del_llh[aluBegin] = log( tmp_r);
      } else {
	dl->second +=  log( tmp_r);
      }
    }
    fin.close();
  }
  
  for (dl = del_llh.begin(); dl != del_llh.end(); dl++) {
    if (dl->second > 0) fout << chrx << " " << dl->first << " " << dl->second * 2. << " " << del_count[dl->first] << endl;
  }
}

void check_region(string &bam_file, string &chrx, int beginPos, int endPos){  
}

int main( int argc, char* argv[] ) {
  if (argc < 2) exit(1);
  string config_file = argv[1];

  string file_pn_delete_prefix = read_config(config_file, "file_pn_delete_prefix");
  string f_out = read_config(config_file, "test_out");
  vector <string> pns;
  string pn;
  ifstream fin(read_config(config_file, "test_pn").c_str());
  stringstream ss;
  while (fin >> pn)   pns.push_back(pn);
  fin.close();

  ofstream fout(f_out.c_str());
  for ( int ci = 1; ci < 23; ci ++ ){
  //  for ( int ci = 1; ci < 2; ci ++ ){
    ss.str("");
    ss << ci;
    string chrn = "chr" + ss.str();
    combine_multi(pns, file_pn_delete_prefix, fout, chrn);    
  }
  cerr << "output to " << f_out << endl;
  fout.close();   
}
