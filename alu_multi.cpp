#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "utils.h"

void combine_multi(vector <string> &pns, string &file_delete_multi_prefix, string &f_out, string chrx){
  string _chrx, tmp1, tmp2;
  double p0, p1, p2;
  int aluBegin;
  map <int, int> del_count;
  map <int, vector <double> > del_debug;
  map <int, double> del_llh;
  ifstream fin;
  for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
    fin.open( (file_delete_multi_prefix + *pi).c_str() );
    assert(fin);
    while (fin >> _chrx >> tmp1 >> aluBegin >> tmp2 >> p0 >> p1 >> p2) {
      if (chrx != _chrx) continue;
      if (p1 > p0 or p2 > p0) addKey(del_count, aluBegin, p1 > p2 ? 1 : 2);
    }
    fin.close();
  }
  //print_map(del_count);
  float pn_chr = 2. * pns.size();
  float fp, f0, f1, f2;
  map <int, int>::iterator dc;
  map <int, double>::iterator dl;
  for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
    fin.open( (file_delete_multi_prefix + *pi).c_str() );
    while (fin >> _chrx >> tmp1 >> aluBegin >> tmp2 >> p0 >> p1 >> p2) {
      if (chrx != _chrx) continue;
      if ( (dc = del_count.find(aluBegin)) == del_count.end() ) continue;      
      if ( (dl = del_llh.find(aluBegin)) == del_llh.end() ) {
	fp = dc->second / pn_chr;
	f0 = (1 - fp) * (1 - fp);
	f1 = 2 * fp * (1 - fp);
	f2 = fp * fp;
	del_llh[dc->first] = log( (f0 * p0 + f1 * p1 + f2 * p2) / p0);
      } else {
	dl->second +=  log( (f0 * p0 + f1 * p1 + f2 * p2) / p0 );
      }
      
      del_debug[dc->first].push_back( (f0 * p0 + f1 * p1 + f2 * p2) / p0 );
    }
    fin.close();
  }
  
  for (dl = del_llh.begin(); dl != del_llh.end(); dl++) {
    cerr << dl->first << " " << dl->second * 2. << " " << del_count[dl->first] << endl; 
    print_vec(del_debug[dl->first]);
  }
}

void check_region(string &bam_file, string &chrx, int beginPos, int endPos){  
}

int main( int argc, char* argv[] ) {
  if (argc < 2) exit(1);
  string config_file = argv[1];
  int n_pn;
  seqan::lexicalCast2(n_pn, read_config(config_file, "n_pn"));
  string file_delete_multi_prefix = read_config(config_file, "file_delete_multi_prefix");
  string f_out = file_delete_multi_prefix + "tmp" + int_to_string(n_pn);

  vector <string> pns;
  string pn;
  ifstream fin(read_config(config_file, "file_pn").c_str());
  int i = 0;
  while (fin >> pn) {
    if (i++ >= n_pn) break;
    pns.push_back(pn);
  }
  
  combine_multi(pns, file_delete_multi_prefix, f_out, "chr21");  
}
