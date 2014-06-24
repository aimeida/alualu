#define SEQAN_HAS_ZLIB 1
#include "utils.h"
#include "consensus.h"

typedef seqan::Dna5String TSeq;

void read_seq(string pn, string fn, map <string, vector<string> > & pos_seqs){
  ifstream fin( fn.c_str() );
  assert (fin);
  string line, pos, seq , readtype , correct_byRC , qName ;
  stringstream ss;
  getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pos >> seq >> readtype >> correct_byRC >> qName;
    pos_seqs[pos].push_back( pn + " " +  seq + " " + readtype+ " " + correct_byRC+ " " + qName );
  }
  fin.close();
}

int main( int argc, char* argv[] )
 {
   string config_file = argv[1];
   string opt =  argv[2];
   if (argc < 3) exit(1);
   
   boost::timer clocki;    
   clocki.restart();

   ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
   vector<string> chrns;  
   for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
   chrns.push_back("chrX");
   chrns.push_back("chrY");  

#ifdef DEBUG_MODE   // only look at chr1 for now 
  cerr << "!!!! DEBUG_MODE: only chr1 tested\n";
  chrns.clear();
  chrns.push_back("chr1");
#endif       

  string path_cons = cf_fh.get_conf("file_ins_cons");    
  std::set <string> pns_used;
  read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
  
  if (opt == "cons_reads_pns" ) {  // rewrite reads to another folder 
    
    for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
      string path2 = path_cons + *ci + "_pos/";
      check_folder_exists( path2 );
      map <string, vector<string> > pos_seqs;
      for (std::set <string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) 
	read_seq(*pi, path_cons + *ci + "/" + *pi, pos_seqs);
      if (pos_seqs.empty()) continue;      
      for (map <string,  vector<string> >::iterator pi = pos_seqs.begin(); pi != pos_seqs.end() ; pi++ ) {
	string file_cons1 = path2 + pi->first;
	ofstream fout1(file_cons1.c_str());
	for (  vector< string > ::iterator si = (pi->second).begin(); si != (pi->second).end(); si++ )
	  fout1 << *si << endl;
	sort_file_by_col<string> (file_cons1, 1, true); // sort by pn
	fout1.close();		
      }
    }
  } else if (opt == "cons_reads_build" ) {  // eg: /home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/cons/chr1_pos/159408
    
  }
  
  cout << "total time used " << clocki.elapsed() << endl;
  clocki.restart();  
  return 0;
 }

