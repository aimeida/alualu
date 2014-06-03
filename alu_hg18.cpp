// print out some known alu sequence in hg18
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/graph_msa.h>
#include <seqan/align.h>
#include "common.h"
#include "utils.h"

#define N_ALU 100

int print_align( map<int, seqan::CharString> &fa_seqs, int i, int j){
  TAlign align;
  seqan::Score<int> scoringScheme(1, -3, -2, -5);
  resize(rows(align), 2);
  assignSource(row(align,0), fa_seqs[i]);
  assignSource(row(align,1), fa_seqs[j]); 
  
  int score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>()); 
  cout << align << endl;
}

void multi_align( map<int, seqan::CharString> &fa_seqs, size_t nseq, ofstream &fout){
  size_t i = 0;
  size_t n_align = min(nseq, fa_seqs.size());
  cerr << "n_align " << n_align << endl;
  map<int, seqan::CharString>::iterator fi;
  seqan::StringSet<seqan::CharString> seq;
  for (fi = fa_seqs.begin(); fi != fa_seqs.end(), i < n_align; fi++, i++) 
    appendValue(seq, fi->second);
  seqan::Graph<seqan::Alignment<seqan::StringSet<seqan::CharString, seqan::Dependent<> > > > aliG(seq);
  seqan::Score<int> scoringScheme(2, -3, -3, -7);
  globalMsaAlignment(aliG, scoringScheme);
  fout << aliG << std::endl;
}

int main( int argc, char* argv[] )
{
  string config_file = argv[1];
  string opt = argv[2];
  string align_alu_type = argv[3]; // AluY, AluSx, AluJo, AluJb
  string chrn = "chr1";

  ConfigFileHandler cf_fh = ConfigFileHandler(config_file);

  string file_alupos_prefix = cf_fh.get_conf("file_alupos_prefix"); 
  string file_fa = cf_fh.get_conf("file_fa_prefix") + chrn + ".fa";
  
  AluRefPosRead *alurefpos = new AluRefPosRead(file_alupos_prefix + chrn, 200);    
  map<int, seqan::CharString> fa_seqs; 
  int aluBegin, aluEnd;
  char chain;
  string alu_type;

  
  if (opt == "1") {
	  string file_alu = "/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/grocery/alu.seq."+align_alu_type;  
	  ofstream fout(file_alu.c_str());
	  for (int i=0; i < N_ALU; ) {
	    alurefpos->updatePos(aluBegin, aluEnd, chain, alu_type);
	    if (alu_type != align_alu_type ) continue;
	    seqan::CharString fa = FastaFileHandler::fasta_seq(file_fa, chrn, aluBegin, aluEnd);
	    if ( chain=='-' )  seqan::reverseComplement(fa);
	    fa_seqs[i] = fa;
	    fout << chrn << " " << aluBegin << " " << aluEnd << " " << fa << endl;
	    i++;
	  }  
	  
	//  multi_align(fa_seqs, 10, fout);
	//  print_align(fa_seqs, 0, 1);
	//  print_align(fa_seqs, 0, 2);
	//  print_align(fa_seqs, 1, 2);
	  fout.close();
 	  cerr << "write into " << file_alu << endl;
  } else if ( opt == "2" ) {
    string file_alu = "/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/grocery/alu.seq." + align_alu_type + ".fa";  
    ofstream fout(file_alu.c_str());
    int i = 0, skip_first_n = 50;
    while ( i < N_ALU + skip_first_n) {      
      alurefpos->updatePos(aluBegin, aluEnd, chain, alu_type);
      if (alu_type != align_alu_type ) continue;
      i++;
      if (i < skip_first_n) continue;
      seqan::CharString fa = FastaFileHandler::fasta_seq(file_fa, chrn, aluBegin, aluEnd);
      if ( chain=='-' )  seqan::reverseComplement(fa);
      string fa_name = chrn + "_" + int_to_string(aluBegin) + "_" + int_to_string(aluEnd);
      writeRecord(fout, fa_name.c_str(), fa, seqan::Fasta());      
    }      
    fout.close();
    cerr << "write into " << file_alu << endl;
  }	  
  
  return 0;  
}
