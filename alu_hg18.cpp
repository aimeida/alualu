// print out some known alu sequence in hg18
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "utils.h"

#define N_ALU 20

int print_align( map<int, seqan::CharString> &fa_seqs, int i, int j){
  TAlign align;
  seqan::Score<int> scoringScheme(1, -3, -2, -5);
  resize(rows(align), 2);
  assignSource(row(align,0), fa_seqs[i]);
  assignSource(row(align,1), fa_seqs[j]); 
  
  int score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>()); 
  cout << align << endl;
}

int main( int argc, char* argv[] )
{
  string chrn = "chr13";
  string config_file = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/config.properties";
  string file_alu = "/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/grocery/alu.seq";

  string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
  string file_fa_prefix = read_config(config_file, "file_fa_prefix");
  seqan::FaiIndex faiIndex;
  unsigned fa_idx;
  assert ( !read(faiIndex, (file_fa_prefix + chrn + ".fa").c_str()) );      
  assert ( getIdByName(faiIndex, chrn, fa_idx) );

  AluRefPosRead *alurefpos = new AluRefPosRead(file_alupos_prefix + chrn, 200);    
  map<int, seqan::CharString> fa_seqs; 
  int aluBegin, aluEnd;
  char chain;
  ofstream fout(file_alu.c_str());
  for (int i=0; i < N_ALU; i++) {
    alurefpos->updatePos(aluBegin, aluEnd, chain);
    seqan::CharString fa = fasta_seq(faiIndex, fa_idx, aluBegin, aluEnd, true); 
    if ( chain=='-' )  seqan::reverseComplement(fa);
    fa_seqs[i] = fa;
    fout << chrn << " " << aluBegin << " " << aluEnd << " " << fa << endl;
  }
  
  print_align(fa_seqs, 0, 1);
  print_align(fa_seqs, 0, 2);
  print_align(fa_seqs, 1, 2);

  fout.close();
  cerr << "write into " << file_alu << endl;
  return 0;  
}
