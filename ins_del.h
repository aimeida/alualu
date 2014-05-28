//#define DEBUG_MODE  // test only chr1 for now
#define REF_EXT_LEN 150 // extend exact insert pos to both sides. This data, max read len = 135
#define ALU_FLANK 400  // 600 (alu_delete.cpp) - 200 (min alu length)
#define DISCORDANT_LEN 2000
#define CLIP_BP 10
#define ALUCONS_SIMILARITY 0.8
#define MIN_READS_CNT 3 // min interesting reads, in order to consider this pos
#define MIN_MATCH_LEN 60 // in order to use this read

enum I_READ {alu_read, aluPart_read, aluSkip_read, unknow_read, useless_read};

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/basic.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/modifier.h>

#include "common.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

class RecordInfo {
 public:
  string pairChrn;
  seqan::CharString qName;
  int insertBegin, insertEnd, thisBegin, thisEnd, pairBegin, pairEnd, rgIdx, insertLen;
  bool thisRC, pairRC;
 RecordInfo(string pc, seqan::CharString qn, int ib, int ie, int tb, int te, int pb, int pe, int idx, int il, bool tr, bool pr) 
   : pairChrn(pc), qName(qn), insertBegin(ib), insertEnd(ie), thisBegin(tb), thisEnd(te),
    pairBegin(pb), pairEnd(pe), rgIdx(idx), insertLen(il), thisRC(tr), pairRC(pr) {}
  static void delete_record_list(list <RecordInfo *> & records); 
  static bool sort_faPos(const RecordInfo* a, const RecordInfo* b);
  static bool sort_insertPos(const RecordInfo* a, const RecordInfo* b);
  static void debugprint_thisBegin(list <RecordInfo *> & records, ostream & os, size_t n = 0);
};

class qNameInfo {
 public:
  I_READ itype;
  int insertAluOri;
  seqan::CharString seq; 
 qNameInfo(I_READ t, int ori, seqan::CharString & sq) : itype(t), insertAluOri(ori), seq(sq) {}    
};

class ConfigFileHandler{
 public:
  map <string, string> configs;
  ConfigFileHandler(string config_file);
  string get_cf_val(string key);
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

class AluconsHandler : public FastaFileHandler {
 public:
  string seq_name;
  int seq_len;
  AluconsHandler(string fn_fa, string sn);
  void update_seq_name(string sn);
  seqan::CharString fetch_alucons(int key);
 private:
  map <int, seqan::CharString> seqs;
};
