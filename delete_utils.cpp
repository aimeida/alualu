#include "delete_utils.h"

bool get_align_pos(int aluBegin, int aluEnd, int beginPos, int endPos, int &ref_a, int &ref_b, int &read_a, int &read_b, seqan::BamAlignmentRecord &record){
  unsigned nl = length(record.cigar) - 1; // ignore complicated alignment between S and M
  int len_read = length(record.seq);
  int len_clip_toAlign = 0;
  if ( abs(beginPos - aluEnd) <= BOUNDARY_OFFSET and record.cigar[0].operation == 'S') {  //eg. S41M60, has_soft_first 
    len_clip_toAlign = record.cigar[0].count;
    ref_a = aluBegin - len_clip_toAlign;
    ref_b = aluBegin;
    read_a = 0;
    read_b = len_clip_toAlign;
  } else if (abs(endPos - aluBegin) <= BOUNDARY_OFFSET and record.cigar[nl].operation == 'S' ) { // eg. M27S93, has_soft_last
    len_clip_toAlign = record.cigar[nl].count;
    ref_a = aluEnd;
    ref_b = aluEnd + len_clip_toAlign;
    read_a = len_read - len_clip_toAlign;
    read_b = len_read;
  }
  if ( !len_clip_toAlign ) return false;
  /* beginning of reads is bad quality. we cut a few bp before aligning
     eg: 10bp cut 1bp, 110bp(90bp) cut 10bp;*/
  int cut_bp = len_read > 110 ? round( 0.11 * len_clip_toAlign - 0.13 ) : round( 0.1 * len_clip_toAlign ); 
  if ( hasFlagRC(record) and read_a == 0) read_a += cut_bp;    
  if ( (!hasFlagRC(record)) and read_b == len_read) read_b -= cut_bp;
  return read_a < read_b;
}

bool split_global_align(seqan::CharString &fa_seq, seqan::BamAlignmentRecord &record, int read_a, int read_b){
  // no need to use seed or chain, since don't know exact position to start mapping
  seqan::Score<int> scoringScheme(1, -3, -2, -5); // match, mismatch, gap extend, gap open
  TAlign align;
  int align_len, score;
  int min_align_len = max(CLIP_BP, (int) (0.85 * (read_b - read_a - BOUNDARY_OFFSET) ) ); 
  resize(rows(align), 2);
  assignSource(row(align,0), infix(record.seq, read_a, read_b) );  // 2,3
  assignSource(row(align,1), fa_seq);   // 1,4

  // free gap both ends
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>()); 
  align_len = min(toViewPosition(row(align, 0), read_b - read_a), toViewPosition(row(align, 1), length(fa_seq)))
    - max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  //cout << clippedBeginPosition(row(align, 0)) << " " << clippedBeginPosition(row(align, 1)) << " " <<  endl;
  //cout << "score1 " << score << " " << align_len << " " << min_align_score(align_len) << " " << min_align_len << " " << endl;
  if ( align_len >= min_align_len and score >= min_align_score(align_len) ) return true; 
  
  // AlignConfig 2=true: free end gaps for record.seq
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<false, true, true, false>()); 
  align_len = min(toViewPosition(row(align, 0), read_b - read_a), toViewPosition(row(align, 1), length(fa_seq)))
    - max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  return ( align_len >= min_align_len and score >= min_align_score(align_len) );
}

T_READ classify_read(seqan::BamAlignmentRecord & record, int align_len, int aluBegin, int aluEnd, FastaFileHandler *fasta_fh, bool read_is_left){

  ////return unknow_read;  // just try it for fun! not serious 
  int beginPos = record.beginPos;
  int endPos = record.beginPos + align_len;  
  // 
  if ( (!read_is_left) and beginPos > aluEnd - BOUNDARY_OFFSET and  
       (!has_soft_first(record, CLIP_BP)) and record.pNext + DEFAULT_READ_LEN < aluBegin + BOUNDARY_OFFSET )
    return unknow_read;
  // left read, the pair cross alu, no need to check its right pair afterwards
  if ( read_is_left and endPos < aluBegin - BOUNDARY_OFFSET and  
       (!has_soft_last(record, CLIP_BP)) and record.pNext >=  aluEnd + ALU_FLANK )
    return unknow_read;
    
  if ( (has_soft_first(record, CLIP_BP) and abs(beginPos - aluEnd) <= BOUNDARY_OFFSET ) or 
       ( has_soft_last(record, CLIP_BP) and abs(endPos - aluBegin) <= BOUNDARY_OFFSET ) ) {    
    int ref_a, ref_b, read_a, read_b;
    if ( get_align_pos(aluBegin, aluEnd, beginPos, endPos, ref_a, ref_b, read_a, read_b, record)) {
      seqan::CharString fa_seq;
      fasta_fh->fetch_fasta_upper(ref_a - BOUNDARY_OFFSET, ref_b + BOUNDARY_OFFSET, fa_seq);          if (split_global_align(fa_seq, record, read_a, read_b)) return clip_read;
    }
    if (read_is_left) return useless_read; 
  }
  // only consider as mid_read if very certain, otherwise classify as unknown read
  if ( endPos > aluBegin + BOUNDARY_OFFSET and beginPos < aluEnd - BOUNDARY_OFFSET) 
    if (!not_all_match(record))  return mid_read;  
  
  /* no big difference if ignore this part
    if ( (!read_is_left) and beginPos > aluEnd - BOUNDARY_OFFSET and 
    record.pNext + DEFAULT_READ_LEN < aluBegin + BOUNDARY_OFFSET )
    return unknow_read;
  */  
  return useless_read;  // alu_flank is too large, we have a lot reads not useful 
}

