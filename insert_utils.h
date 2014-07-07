#define SEQAN_HAS_ZLIB 1
#include "utils.h"
#ifndef INSERT_UTILS_H
#define INSERT_UTILS_H

typedef map<int, ofstream* > MapFO;
typedef map<string, ofstream* > MapSFO;
typedef seqan::Dna5String TSeq;

inline void close_fhs(MapFO & fileMap, map<int, string> & rID_chrn) {
  for (map<int, string>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++)   
    delete fileMap[rc->first];
}

class AlumateINFO {
 public:
  string qname;
  int len_read, rid1, rid2, pos1, pos2, rgIdx;
  bool RC1, RC2;
 AlumateINFO( string & qn, int lr, int id1, int id2, int p1, int p2, int rgIdx, bool r1, bool r2) 
   : qname(qn), len_read(lr), rid1(id1), rid2(id2), pos1(p1), pos2(p2), rgIdx(rgIdx), RC1(r1), RC2(r2) {}
  static bool sort_pos2(const AlumateINFO* a, const AlumateINFO* b);
  static void delete_list(list <AlumateINFO *> & alumate_list);
};

class READ_INFO {
 public:
  int beginPos, endPos;
  bool should_be_left; // inferred from RC flag
  string alu_type;
 READ_INFO(int p, int lr, bool sbl, string alu_type) : beginPos(p), endPos(p+lr-2), should_be_left(sbl), alu_type(alu_type) {}
};

class ALUREAD_INFO {
public:
  seqan::CharString qName;
  int clipLeft, pos, readLen;
  bool sameRC;
 ALUREAD_INFO(seqan::CharString & qn, int cl, int p, int rlen, bool t) : qName(qn), clipLeft(cl), pos(p), readLen(rlen), sameRC(t) {}
};


bool clipRight_move_left(seqan::CharString & read_seq, seqan::CharString & ref_fa, list <int> & cigar_cnts, int refBegin, int & clipPos, int & align_len);
bool clipLeft_move_right(seqan::CharString & read_seq, seqan::CharString & ref_fa, list <int> & cigar_cnts, int refBegin, int & clipPos, int & align_len);

bool read_match_clipLeft(string & line, int clipLeft, string & pn, string & qName);
bool read_match_clipLeft(string & line, int clipLeft, string & pn, string & qName, string & cigar, string & seq); 

bool parseline_del_tmp1(string &line, string & output_line, map <int, EmpiricalPdf *> & pdf_rg, int cnt_alumate, int insertLenPlus = 0);
bool global_align_insert(const int hasRCFlag, seqan::CharString & seq_read, seqan::CharString & seq_ref, int &score, int cutEnd, float th_score, bool verbose = false);
bool align_clip_to_ref(char left_right, int adj_clipPos,  int clipPos, int align_len, seqan::BamAlignmentRecord &record, FastaFileHandler *fasta_fh, ofstream &fout, string  header);
int align_clip_to_LongConsRef(string shortSeq, string longSeq, int & refBegin, int & refEnd, int clipLen);  // consensus sequence is quite long 
void align_clip_to_consRef(string shortSeq, string longSeq, int & refBegin, int & refEnd, int clipLen);
bool align_alu_to_consRef(const string & shortSeq, const string & longSeq, float dif_th, const string & atype) ;


// following, depreciated
string parse_alu_type(string alu_name);
void filter_outlier_pn(string path_input, string fn_suffix, map<int, string> &ID_pn, string chrn, string file_pn_used_output, float percentage_pn_used);

#endif /*INSERT_UTILS_H*/
