#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include "consensus.h"

typedef seqan::Dna5String TSeq;

bool read_Tseq(string fn, seqan::StringSet <TSeq> & seqs, int & clipPos, AluconsHandler *alucons_fh, bool allow_duplicate, string file_seq) {
  ifstream fin( fn.c_str() );
  ofstream fseq;
  float sim_rate;
  if (file_seq!="") fseq.open(file_seq.c_str() );
  string line, seq, readtype, pn; 
  map <string, int> aluRead_cnt;
  map <string, int> clipRead_cnt;
  int p;
  bool sameRC;
  stringstream ss;
  if (!fin) return false;
  clipPos = 0;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> p >> seq >> readtype >> pn;
    if (!clipPos) clipPos = p;
    assert ( p == clipPos);
    if ( readtype == "clipRead" ) 
      addKey(clipRead_cnt, seq, 1);
    else {
      ss >> sameRC ;
      if (sameRC)  reverseComplement(seq); 
      addKey(aluRead_cnt, seq, 1);
    }
  }
  fin.close();

  for (map <string, int>::iterator si = aluRead_cnt.begin(); si != aluRead_cnt.end(); si++) {
    seq = si->first;
    if ( !align_alu_cons_call( seq, alucons_fh, sim_rate, ALUCONS_SIMILARITY) ) continue;
    fseq << p << " " << seq << " " << readtype << " " << pn << " " << sameRC << endl;
    if (allow_duplicate) 
      for ( int i = 0; i < si->second; i++) appendValue(seqs, (TSeq) seq );
    else
      appendValue(seqs, (TSeq) seq );
  }
  //cout << aluRead_cnt.size() << " " << length(seqs) << endl;  
  for (map <string, int>::iterator si = clipRead_cnt.begin(); si != clipRead_cnt.end(); si++) {
    seq = si->first;
    fseq << p << " " << seq << " " << readtype << " " << pn << " " << sameRC << endl;
    if (allow_duplicate) 
      for ( int i = 0; i < si->second; i++) appendValue(seqs, (TSeq) seq );
    else
      appendValue(seqs, (TSeq) seq );    
  }        
  return true;
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

   if (opt == "cons_reads_build" ) { 
     string file_input = argv[3];
     string file_seq = ""; //file_input + ".seqs";
     string file_output = file_input + ".cons";
     int clipPos;
     seqan::StringSet <TSeq>  pos_seqs;    
     AluconsHandler *alucons_fh = new AluconsHandler(cf_fh.get_conf("file_alu_cons"), cf_fh.get_conf("type_alu_cons"));
     bool allow_duplicate = true;
     cout << "allow_duplicate ? " << (allow_duplicate ? "yes\n" : "no\n");
     if ( !read_Tseq(file_input, pos_seqs, clipPos, alucons_fh, allow_duplicate, file_seq) ) return 1;
     if ( length(pos_seqs) == 0 ) return 1;
     cout << "seqs " << length(pos_seqs) << endl;
     seqan::StringSet<TSeq> consensus;
     int multi_align_fail = compute_consensus(consensus, pos_seqs);	
     cout <<  "a multiple read alignment failed ? "<< multi_align_fail << endl;
     ofstream fout(file_output.c_str());

     int max_len = length(consensus[0]);
     int max_id = 0;
     for ( size_t i = 1; i < length(consensus); ++i) 
       if (length(consensus[i]) > max_len) {
	 max_len = length(consensus[i]);
	 max_id = i;
       }

     fout << ">" << clipPos << "\n"<< consensus[max_id] << "\n";  

     for ( size_t i = 0; i < length(consensus); ++i)  
       fout << ">" << clipPos << "_" << i << "\n"<< consensus[i] << "\n";  // ok, i might want to debug, why so many contigs ?	
     fout.close();     
     delete alucons_fh;
   }
   
   cout << "total time used " << clocki.elapsed() << endl;
   clocki.restart();  
   return 0;
 }

