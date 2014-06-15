//#define DEBUG_MODE  // test only chr1 for now
#include "utils.h"

// aluRead: completely mapped to alu (maybe another chr in the genome)
// aluClipRead: clip read, partly mapped to ref
// aluskipRead: align perfectly through alu insert positions (called mid read in alu_delete)
// unknowRead: might be aluskip_read or aluClip_read, use insert length info
enum I_READ {aluRead, aluClipRead, aluSkipRead, unknowRead, uselessRead}; // used for insertions

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
  int insertAluOri;  // take value [1,2,3,4]
  seqan::CharString seq; 
 qNameInfo(I_READ t, seqan::CharString & sq) : itype(t), seq(sq), insertAluOri(0) {}    
};

