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

inline string get_name_rg(string prefix, string pn){ return prefix + "RG." + pn;}
inline string get_name_rg_pdf(string prefix, string pn, string rg, string pdf_param){ return prefix + pn + ".count." + rg + "." + pdf_param; }

inline void check_folder_exists(string & path) {
  if ( access( path.c_str(), 0 ) != 0 ) system( ("mkdir " + path).c_str() );    
};

inline bool left_read( seqan::BamAlignmentRecord &record){return (record.beginPos < record.pNext);};

inline bool QC_read( seqan::BamAlignmentRecord &record){  
  return (!hasFlagQCNoPass(record)) and (!hasFlagDuplicate(record)) and hasFlagMultiple(record) and (!hasFlagSecondary(record));
}

inline bool QC_delete_read( seqan::BamAlignmentRecord &record){  
  return QC_read(record) and hasFlagAllProper(record) and (abs(record.tLen) <= DISCORDANT_LEN) and (hasFlagRC(record) != hasFlagNextRC(record));
};

inline bool QC_insert_read( seqan::BamAlignmentRecord &record){  
  return QC_read(record) and (!hasFlagUnmapped(record)) and (!hasFlagNextUnmapped(record));
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

inline int count_non_match(seqan::BamAlignmentRecord &record){ 
  int non_match_len = 0;
  for (size_t i = 0; i < length(record.cigar); i++) 
    if (record.cigar[i].operation != 'M') non_match_len += record.cigar[i].count; 
  return non_match_len;
}

inline bool p00_is_dominant(float * log10_p, int min_log10p) { return  max( log10_p[2] - log10_p[0], log10_p[1] - log10_p[0]) <= min_log10p; }
inline bool p11_is_dominant(float * log10_p, int min_log10p) { return  max( log10_p[0] - log10_p[2], log10_p[1] - log10_p[2]) <= min_log10p; }
inline string phred_log (float p) { return p ? (int_to_string (-(int)(log10 (p) * 10)) ) : "255"; }

class EmpiricalPdf
{
  map <int, float> prob_vec; 
  int min_len, max_len, bin_width;
  float min_prob;
 public:
  EmpiricalPdf(string pdf_file);
  float pdf_obs(int insertlen);
  static void delete_map(map <int, EmpiricalPdf *> & epdf_rg); 
};

class BamFileHandler{
 public:
  string fn_bai;
  string fn_bam_out;
  map<int, string> rID_chrn ;
  BamFileHandler(vector<string> &chrns, string bam_input, string bai_input, string bam_output="");
  ~BamFileHandler(void);
  bool jump_to_region(int rid, int region_begin, int region_end);
  bool jump_to_region(string chrn, int region_begin, int region_end);
  bool fetch_a_read(seqan::BamAlignmentRecord & record);
  string fetch_a_read(string chrn, int region_begin, int region_end, seqan::BamAlignmentRecord & record);
  bool get_chrn(int query_rid, string & pairChrn);
  bool write_a_read(seqan::BamAlignmentRecord & record);  
 private:
  TNameStore  nameStore;
  TNameStoreCache nameStoreCache;
  TBamIOContext context;
  seqan::Stream<seqan::Bgzf> inStream ;
  map<seqan::CharString, int> chrn_rID ;
  seqan::BamIndex<seqan::Bai> baiIndex ;    
  seqan::Stream<seqan::Bgzf> outStream ;
};

class FastaFileHandler {
 public:
  seqan::FaiIndex faiIndex;
  unsigned idx;
  bool only_one_seq;
  FastaFileHandler(string fn_fa);
  FastaFileHandler(string fn_fa, string seq_name);
  void fetch_fasta_upper(int beginPos, int endPos, seqan::CharString & seq);
  void fetch_fasta_upper(int beginPos, int endPos, seqan::CharString & seq, string seq_name);
  static seqan::CharString fasta_seq(string fa_input, string seq_name,int beginPos, int endPos);
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
};

void read_pdf_pn( string prefix, string pn, string pdf_param,  map <int, EmpiricalPdf *> & pdf_rg);
void get_min_value(map <int, float> & m, float & min_val, int & min_key);
void log10P_to_P(float *log_gp, float *gp, int max_logp_dif);
string phred_scaled(float p0, float p1, float p2);
void parse_reading_group(string file_rg, map<string, int> & rg_to_idx);

int get_rgIdx(map<string, int> & rg_to_idx, seqan::BamAlignmentRecord & record);
void parse_cigar(string cigar, list <char> & opts, list <int> & cnts);
string get_cigar(seqan::BamAlignmentRecord &record);
void debug_print_read(seqan::BamAlignmentRecord &record, ostream & os = cout);
bool find_read(string &bam_input, string &bai_input, string &chrn, string &this_qName, int this_pos, seqan::BamAlignmentRecord &that_record, int flank_region);
int numOfBestHits(seqan::BamAlignmentRecord &record);
void read_first2col(string fn, vector < pair<int, int> > & insert_pos);


#endif /*UTILS_H*/
