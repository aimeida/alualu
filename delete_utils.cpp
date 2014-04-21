#include "delete_utils.h"

string phred_scaled(float p0, float p1, float p2) {
  float pmax = max(p0, max(p1, p2));
  string s0 =  ( p0 == pmax ) ?  "0" : phred_log(p0/pmax);
  string s1 =  ( p1 == pmax ) ?  "0" : phred_log(p1/pmax);
  string s2 =  ( p2 == pmax ) ?  "0" : phred_log(p2/pmax);
  return s0  + "," + s1 + "," + s2; 
}

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

string enum_to_string(T_READ x){
  if ( x == useless_read ) return "useless_read";
  if ( x == unknow_read ) return "unknow_read";
  if ( x == mid_read ) return "mid_read";
  if ( x == clip_read ) return "clip_read";
  return "ERROR!" ;
};

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
  //cout << align << endl;
  if ( align_len >= min_align_len and score >= min_align_score(align_len) ) return true; 
  
  // AlignConfig 2=true: free end gaps for record.seq
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<false, true, true, false>()); 
  align_len = min(toViewPosition(row(align, 0), read_b - read_a), toViewPosition(row(align, 1), length(fa_seq)))
    - max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  return ( align_len >= min_align_len and score >= min_align_score(align_len) );
}


T_READ classify_read(seqan::BamAlignmentRecord & record, int align_len, int aluBegin, int aluEnd, seqan::FaiIndex &faiIndex, unsigned fa_idx, bool read_is_left){  
  int beginPos = record.beginPos;
  int endPos = record.beginPos + align_len;  
  // more unknow_read, more information used 
  if ( (!read_is_left) and beginPos > aluEnd - BOUNDARY_OFFSET and  
       (!has_soft_first(record, CLIP_BP)) and record.pNext + DEFAULT_READ_LEN < aluBegin + BOUNDARY_OFFSET )
    return unknow_read;
  // left reads, and we won't check its right pair afterwards
  if ( read_is_left and endPos < aluBegin - BOUNDARY_OFFSET and  
       (!has_soft_last(record, CLIP_BP)) and record.pNext >=  aluEnd + ALU_FLANK )
    return unknow_read;
    
  if ( (has_soft_first(record, CLIP_BP) and abs(beginPos - aluEnd) <= BOUNDARY_OFFSET ) or 
       ( has_soft_last(record, CLIP_BP) and abs(endPos - aluBegin) <= BOUNDARY_OFFSET ) ) {    
    int ref_a, ref_b, read_a, read_b;
    if ( get_align_pos(aluBegin, aluEnd, beginPos, endPos, ref_a, ref_b, read_a, read_b, record)) {
      seqan::CharString fa_seq = fasta_seq(faiIndex, fa_idx, ref_a - BOUNDARY_OFFSET, ref_b + BOUNDARY_OFFSET, true);
      if (split_global_align(fa_seq, record, read_a, read_b)) return clip_read;
    }
    return useless_read; // if left reads, we don't know the type yet
  }
  // only consider as mid_read if very certain, otherwise look at its pair later
  if ( endPos > aluBegin + BOUNDARY_OFFSET and beginPos < aluEnd - BOUNDARY_OFFSET) 
    if (!not_all_match(record))  return mid_read;  
  return useless_read;  // alu_flank is too larger, we have a lot reads not useful 
}

void get_min_value(map <int, float> & m, float & min_val, int & min_key) {
  map <int, float>::iterator mi = m.begin();
  min_val = mi->second;
  min_key = mi->first ;
  mi++;
  while ( mi != m.end() ) {
    if ( min_val > mi->second) {
      min_val = mi->second;
      min_key = mi->first;
    }
    mi++;
  }
}

void log10P_to_P(float *log_gp, float *gp, int max_logp_dif){  
  assert( max_logp_dif > 0);
  map <int, float> idx_logp;
  int i, min_idx;
  for ( i = 0; i < 3; i++) idx_logp[i] = log_gp[i];
  float min_logp;
  get_min_value(idx_logp, min_logp, min_idx);
  for ( i = 0; i < 3; i++) {
    if ( idx_logp[i] - min_logp > max_logp_dif ) {
      idx_logp[min_idx] = 1;   // set p[i] = 0 afterwards
      min_logp = 1; 
      break;
    }
  }  
  if (min_logp > 0) { // check the two elements left
    get_min_value(idx_logp, min_logp, min_idx);
    for ( i = 0; i < 3; i++) { 
      if ( idx_logp[i] < 1 and idx_logp[i] - min_logp > max_logp_dif ) {
	idx_logp[min_idx] = 1;
	break;
      }
    }
  }
  get_min_value(idx_logp, min_logp, min_idx);
  //cout << "log " << min_logp << " " << max_logp_dif << " " << endl;
  
  float ratio_sum = 0;
  for ( i = 0; i < 3; i++) 
    if (idx_logp[i] < 1) ratio_sum += pow(10, (idx_logp[i] - min_logp));

  for ( i = 0; i < 3; i++) {
    if (idx_logp[i] > 0) gp[i] = 0;
    else gp[i] = pow(10, (idx_logp[i] - min_logp)) / ratio_sum;}    
}

bool check_delete_region(string const & bam_input, string const &bai_input, string const & fa_input,  string chrn, int beginPos, int endPos ){
  seqan::FaiIndex faiIndex;
  unsigned fa_idx = 0;
  assert (!read(faiIndex, fa_input.c_str()) );      
  assert (getIdByName(faiIndex, chrn, fa_idx));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::Stream<seqan::Bgzf> inStream;  
  open(inStream, bam_input.c_str(), "r");
  seqan::BamIndex<seqan::Bai> baiIndex;
  assert( !read(baiIndex, bai_input.c_str())) ;
  assert ( !readRecord(header, context, inStream, seqan::Bam()) );
  int rID = 0;
  assert ( getIdByName(nameStore, chrn, rID, nameStoreCache) ); // change context ??
  seqan::BamAlignmentRecord record;  
  bool hasAlignments = false;
  if (!jumpToRegion(inStream, hasAlignments, context, rID, beginPos, endPos, baiIndex)) return false;
  if (!hasAlignments) return false;

  int aluBegin = beginPos + 600;
  int aluEnd = endPos - 600;

  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != rID || record.beginPos >= endPos) break;
    if (record.beginPos < beginPos) continue;            
    if ( !QC_delete_read(record) ) continue;
    int align_len = getAlignmentLengthInRef(record);    
    bool read_is_left = left_read(record);
    T_READ rt_val = classify_read( record, align_len, aluBegin, aluEnd, faiIndex, fa_idx, read_is_left);
    if (rt_val != useless_read) {
      cout << aluBegin << " " << aluEnd << " " ;
      print_read(record, cout);
    }
  } 
  return true;
}


EmpiricalPdf::EmpiricalPdf(string pdf_file){ 
  int pos;
  float posp;
  bin_width = 0;
  ifstream fin( pdf_file.c_str());
  assert (fin);
  fin >> pos >> posp;
  prob_vec[pos] = posp;
  min_len = pos;
  min_prob = posp;
  while (fin >> pos >> posp) {
    prob_vec[pos] = posp;
    min_prob = min_prob > posp ? posp : min_prob;      
    if (!bin_width) bin_width = pos - min_len;
  }
  max_len = pos;
  fin.close();
  cerr << "read pdf dist " << pdf_file << endl;
  cerr << min_len << " " << max_len << " " << bin_width << " " << min_prob << endl;
}

float EmpiricalPdf::pdf_obs(int insertlen) {
  if (insertlen >= max_len || insertlen <= min_len) return min_prob;
  int nearby_pos = min_len + (insertlen-min_len)/bin_width*bin_width;
  map<int, float >::iterator it = prob_vec.find(nearby_pos);
  if ( it != prob_vec.end())  return it->second;
  return min_prob;
}   

