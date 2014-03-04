// some data struct to read and write data 
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "common.h"
#ifndef UTILS_H
#define UTILS_H

// Setup name store, cache, and BAM I/O context.
typedef seqan::StringSet<seqan::CharString> TNameStore;
typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
typedef seqan::BamIOContext<TNameStore>     TBamIOContext;

string get_pn(string pn_file, int idx_pn);
string read_config(string config_file, string key);
void print_cigar(seqan::String<seqan::CigarElement<> > &cigar, int len_cigar);
bool find_read(string &bam_input, string &bai_input, string &chrx, string &this_qName, int this_pos, seqan::BamAlignmentRecord &that_record, int flank_region);
void fasta_seq(string fa_input, string chrx, int beginPos, int endPos, seqan::CharString &seq);
void insertLen_of_nonUniq_mapping(vector<int> &starts_ends, vector<int> &reads_insert_len);

class AluRefPos
{
  queue<int> beginP, endP;
 public:
  AluRefPos(string file_alupos);
  int updatePos(int &beginPos, int &endPos);
  ~AluRefPos(void);
};

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
