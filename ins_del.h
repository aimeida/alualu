//#define DEBUG_MODE  // test only chr1 for now
#include "utils.h"

enum I_READ {alu_read, aluPart_read, aluSkip_read, unknow_read, useless_read}; // used for insertions

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
