// some data struct to read and write data 
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/modifier.h>
#include "common.h"
#ifndef UTILS_H
#define UTILS_H

typedef seqan::StringSet<seqan::CharString> TNameStore;
typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
typedef seqan::Align<seqan::CharString, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow; 
typedef seqan::Iterator<TRow>::Type TRowIterator;

struct MyUpper : public unary_function<char,char> {
  inline char operator()(char x) const  {
    if (('a' <= x) && (x <= 'z')) return (x + ('A' - 'a'));
    return x; 
  }
};

class BamFileHandler{
 public:
  string fn_bai;
  string fn_bam_out;
  BamFileHandler(vector<string> &chrns, string bam_input, string bai_input, string bam_output="");
  ~BamFileHandler(void);
  bool jump_to_region(string chrn, int region_begin, int region_end);
  string fetch_a_read(string chrn, int region_begin, int region_end, seqan::BamAlignmentRecord & record);
  bool get_chrn(int query_rid, string & pairChrn);
  bool write_a_read(seqan::BamAlignmentRecord & record);  
 private:
  TNameStore  nameStore;
  TNameStoreCache nameStoreCache;
  TBamIOContext context;
  seqan::Stream<seqan::Bgzf> inStream ;
  map<seqan::CharString, int> chrn_rID ;
  map<int, string> rID_chrn ;
  seqan::BamIndex<seqan::Bai> baiIndex ;    
  seqan::Stream<seqan::Bgzf> outStream ;
};

class FastaFileHandler {
 public:
  seqan::FaiIndex faiIndex;
  FastaFileHandler(string fn_fa);
  void fetch_fasta_upper(string seq_name, int beginPos, int endPos, seqan::CharString & seq);
};

inline void check_folder_exists(string & path) {
  if ( access( path.c_str(), 0 ) != 0 ) system( ("mkdir " + path).c_str() );    
};
inline bool left_read( seqan::BamAlignmentRecord &record){return (record.beginPos < record.pNext);};

inline bool QC_delete_read( seqan::BamAlignmentRecord &record){  
  return ( (!hasFlagQCNoPass(record)) and (!hasFlagDuplicate(record)) and hasFlagMultiple(record) and hasFlagAllProper(record) and abs(record.tLen) <= 2000) and (!hasFlagSecondary(record));
  /// BAM_FLAG_ALL_PROPER: 0x0002 All fragments have been aligned properly.
};

inline bool QC_insert_read( seqan::BamAlignmentRecord &record){  
  return ( (!hasFlagQCNoPass(record)) and (!hasFlagDuplicate(record)) and hasFlagMultiple(record) and (!hasFlagUnmapped(record)) and (!hasFlagNextUnmapped(record)) and (!hasFlagSecondary(record)));
  // what about  ( ! hasFlagAllProper(record) )  ???
};

inline bool has_soft_last(seqan::BamAlignmentRecord &record, unsigned min_bp){ 
  unsigned i = length(record.cigar) - 1;
  return (record.cigar[i].operation == 'S') and (record.cigar[i].count >= min_bp) ;
};
inline bool has_soft_first(seqan::BamAlignmentRecord &record, unsigned min_bp){ 
  return (record.cigar[0].operation == 'S') and (record.cigar[0].count >= min_bp); 
};

inline bool not_all_match(seqan::BamAlignmentRecord &record, int max_err_bp = 5){ 
  int non_match_len = 0;
  for (size_t i = 0; i < length(record.cigar); i++) 
    if (record.cigar[i].operation != 'M') non_match_len += record.cigar[i].count; 
  return non_match_len > max_err_bp;  // if <= 5bp, consider as full match
};


class AluRefPosRead
{
  queue<int> beginP, endP;
  queue<char> strandP;
  queue<string> aluType;
 public:
  int minP, maxP;
  AluRefPosRead(string file_alupos, int minLen = 200);
  int updatePos(int &beginPos, int &endPos);   // alu_delete, alu_hg18 used
  int updatePos(int &beginPos, int &endPos, char & chain, string & alu_type);
};

class AluRefPos   // used by alu_insert, read alu in different way, update more efficient 
{
public:
  vector<int> beginV, endV;
  vector<string> typeV;
  AluRefPos(string file_alupos);   
  bool insideAlu(int beginPos, int endPos, int alu_min_overlap, int &len_overlap); 
  ~AluRefPos(void);
};


class RepMaskPos
{
 public:
  vector<string> chrns;
  map<string, vector<int> > beginP, endP;
  RepMaskPos(string file_rmsk, int join_len = 20);
  void print_begin(int ni=10);
};


void parse_reading_group(string file_rg, map<string, int> & rg_to_idx);
int get_rgIdx(map<string, int> & rg_to_idx, seqan::BamAlignmentRecord & record);

void parse_cigar(string cigar, list <char> & opts, list <int> & cnts);
string read_config(string config_file, string key);
string get_cigar(seqan::BamAlignmentRecord &record);
void print_read(seqan::BamAlignmentRecord &record, ostream & os = cout);
bool find_read(string &bam_input, string &bai_input, string &chrn, string &this_qName, int this_pos, seqan::BamAlignmentRecord &that_record, int flank_region);
seqan::CharString fasta_seq(string fa_input, string chrn, int beginPos, int endPos, bool upper = true);
seqan::CharString fasta_seq(seqan::FaiIndex &faiIndex, unsigned idx, int beginPos, int endPos, bool upper = true);
void get_rID_chrn(string & bam_input, vector<string> &chrns, map<int, seqan::CharString> &rID_chrn);
int numOfBestHits(seqan::BamAlignmentRecord &record);


#endif /*UTILS_H*/
