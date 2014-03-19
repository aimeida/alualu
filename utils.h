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
typedef seqan::Align<seqan::CharString, seqan::ArrayGaps> TAlign;

struct MyUpper : public unary_function<char,char> {
  inline char operator()(char x) const  {
    if (('a' <= x) && (x <= 'z')) return (x + ('A' - 'a'));
    return x; 
  }
};

// without inline, compile error
inline bool left_read( seqan::BamAlignmentRecord &record){return (record.beginPos < record.pNext);};
inline bool QC_read( seqan::BamAlignmentRecord &record){  // if this read is counted in calculating coverage
  return ((not hasFlagQCNoPass(record) ) and hasFlagAllProper(record) and (not hasFlagDuplicate(record)) and hasFlagMultiple(record) );
};
inline bool has_soft_last(seqan::BamAlignmentRecord &record, unsigned min_bp){ 
  unsigned i = length(record.cigar) - 1;
  return (record.cigar[i].operation == 'S') and (record.cigar[i].count >= min_bp) ;
}; 
inline bool has_soft_first(seqan::BamAlignmentRecord &record, unsigned min_bp){ 
  return (record.cigar[0].operation == 'S') and (record.cigar[0].count >= min_bp); 
};

inline bool not_all_match(seqan::BamAlignmentRecord &record, int max_bp = 5){ 
  //return !(length(record.cigar) == 1 and record.cigar[0].operation == 'M');
  int non_match_len = 0;
  for (size_t i = 0; i < length(record.cigar); i++) 
    if (record.cigar[i].operation != 'M') non_match_len += record.cigar[i].count; 
  return non_match_len > max_bp;  // if < 5bp, consider as full match
};

string read_config(string config_file, string key);
void print_cigar(seqan::BamAlignmentRecord &record);
void print_cigar(seqan::BamAlignmentRecord &record, ofstream &fout);

void print_read(seqan::BamAlignmentRecord &record);
void print_read(seqan::BamAlignmentRecord &record, ofstream &fout);

bool find_read(string &bam_input, string &bai_input, string &chrx, string &this_qName, int this_pos, seqan::BamAlignmentRecord &that_record, int flank_region);
bool find_read(seqan::Stream<seqan::Bgzf> &inStream, TBamIOContext &context, int rID, int this_pos, string this_qName, seqan::BamAlignmentRecord &record);

seqan::CharString fasta_seq(string fa_input, string chrx, int beginPos, int endPos, bool upper = true);
seqan::CharString fasta_seq(seqan::FaiIndex &faiIndex, unsigned idx, int beginPos, int endPos, bool upper = true);

void insertLen_of_nonUniq_mapping(vector<int> &starts_ends, vector<int> &reads_insert_len);

void get_rID_chrx(string & bam_input, vector<string> &chrns, map<int, seqan::CharString> &rID_chrx);
char mappingType(seqan::BamAlignmentRecord &record);
int numOfBestHits(seqan::BamAlignmentRecord &record);

class AluRefPos
{
  queue<int> beginP, endP;
  vector<int> beginV, endV;
 public:
  int minP, maxP;
  AluRefPos(string file_alupos, bool use_vector = false);
  bool endOfChr(int p);
  int updatePos(int &beginPos, int &endPos);
  bool insideAlu(int beginPos, int endPos, int alu_min_overlap, int &len_overlap); // return 0 if not inside Alu
  bool insideAlu(int beginPos, int endPos, int alu_min_overlap, int &beginPos_match, int &endPos_match);
  ~AluRefPos(void);
};

class RepMaskPos
{
 public:
  vector<string> chrns;
  map<string, vector<int> > beginP, endP;
  RepMaskPos(string file_rmsk, int join_len = 20);
  void print_begin(int ni=10);
  //bool insideRep(string chrn, int pos);
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
