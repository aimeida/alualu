#include "delete_utils.h"

string debug_print_tread(T_READ td){
  if ( td == unknow_read ) return "unknow_read";
  if ( td == mid_read ) return "mid_read";
  if ( td == clip_read ) return "clip_read";
  return "";
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

T_READ classify_read(seqan::BamAlignmentRecord & record, int aluBegin, int aluEnd, FastaFileHandler *fasta_fh, bool only_tLen_info){
  bool read_is_left = left_read(record);
  int beginPos = record.beginPos;
  int endPos = record.beginPos + getAlignmentLengthInRef(record);    
  int pair_begin = read_is_left ? beginPos : record.pNext;
  int pair_end = pair_begin + abs(record.tLen);

  if (only_tLen_info ) {// silly to use only insert length info. (1) mid reads with small insert length will be wrongly interpreted. (2) many clip reads are ignored 
    if ( pair_begin < aluBegin + BOUNDARY_OFFSET and pair_end > aluEnd - BOUNDARY_OFFSET ) return unknow_read;
    else return useless_read;
  }

  if ( (has_soft_first(record, CLIP_BP) and abs(beginPos - aluEnd) <= BOUNDARY_OFFSET ) or 
       ( has_soft_last(record, CLIP_BP) and abs(endPos - aluBegin) <= BOUNDARY_OFFSET ) ) {    
    int ref_a, ref_b, read_a, read_b;
    if ( get_align_pos(aluBegin, aluEnd, beginPos, endPos, ref_a, ref_b, read_a, read_b, record)) {
      seqan::CharString fa_seq;
      fasta_fh->fetch_fasta_upper(ref_a - BOUNDARY_OFFSET, ref_b + BOUNDARY_OFFSET, fa_seq);          
      if (split_global_align(fa_seq, record, read_a, read_b)) 
	return clip_read;
    }
  }  
  // only consider as mid_read if very certain, otherwise classify as unknown read
  if ( beginPos < aluEnd - BOUNDARY_OFFSET and endPos > aluBegin + BOUNDARY_OFFSET and count_non_match(record) <= 5)       
    return mid_read;    
  if ( pair_begin < aluBegin + BOUNDARY_OFFSET and pair_end > aluEnd - BOUNDARY_OFFSET )
    return unknow_read;  
  return useless_read;  // alu_flank is too large, we have a lot reads not useful 
}

void filter_by_llh_noPrivate(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> &chrns, int col_00) {
  stringstream ss;
  string line, chrn, tmpfield;
  int aluBegin, flag;
  float p0, p1, p2;

  map < int, int > pos_altCnt; // count of alternative alleles 
  map < int, int > pos_pnCnt;  // count of pns having polymorphism at this loci
  map < int, map<string, GENO_PROB * > > pos_pnProb;  
  map < int, map<string, GENO_PROB * > >::iterator pp;
  map < string, GENO_PROB * >::iterator ppi;

  ofstream fout(f_out.c_str());  
  fout <<  "chrn aluBegin pnCnt alleleFreq Log_LikelihoodRatio\n"; 
  int alleleCnt = pns.size() * 2; 
  for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
      ifstream fin( (path0 + *pi + f_in_suffix).c_str() );
      assert(fin);
      getline(fin, line); 
      flag = 1;
      while ( getline(fin, line) ) {
	ss.clear(); ss.str( line );
	ss >> chrn >>  aluBegin;
	for (int ti = 0; ti < col_00 - 3; ti++) ss >> tmpfield;

	ss >> p0 >> p1 >> p2 ;
	  
	if (chrn != *ci) {
	  if (flag) continue;
	  else break;
	}
	flag = 0;
	pos_pnProb[ aluBegin ][*pi] = new GENO_PROB( p0, p1, p2);	  
	if (p1 > p0 or p2 > p0) {
	  addKey(pos_altCnt, aluBegin, p1 > p2 ? 1 : 2);	
	  addKey(pos_pnCnt, aluBegin, 1);
	}
      }
      fin.close();
    }
    
    map<int, int>::iterator pan = pos_pnCnt.begin();
    for ( map<int, int>::iterator pa = pos_altCnt.begin(); pa != pos_altCnt.end(); pa++, pan++) {
      if ( pan->second <= 1) continue; // no private 
      float altFreq = (pa->second) / (float)alleleCnt;
      float freq0 = (1 - altFreq) * (1 - altFreq);
      float freq1 = 2 * altFreq * (1 - altFreq);
      float freq2 = altFreq * altFreq;
      float llh_alt = 0, llh_00 = 0;
      for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
	if ( (ppi = pos_pnProb[pa->first].find(*pi)) != pos_pnProb[pa->first].end() ) {
	  llh_00 +=  ( (ppi->second)->g0 > 0 ? log( (ppi->second)->g0 ) : ( - LOG10_GENO_PROB ) ) ;
	  llh_alt += log (freq0 * (ppi->second)->g0 + freq1 * (ppi->second)->g1 + freq2 * (ppi->second)->g2);  
	} else
	  llh_alt += log (freq0);  // g0 = 1	
      }
      fout << *ci << " " << pa->first << " " << pan->second << " " << altFreq << " " << llh_alt - llh_00 << endl;
    }
    cout << "done with " << *ci << endl;

    pos_altCnt.clear(); 
    pos_pnCnt.clear();
    for (pp = pos_pnProb.begin(); pp != pos_pnProb.end(); pp++) {
      for ( ppi = (pp->second).begin(); ppi != (pp->second).end(); ppi++) 
	delete ppi->second;
    }
    pos_pnProb.clear();
  } 
  fout.close();
}

bool combine_pns_vcf_noPrivate(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> & chrns, int col_00) {
  ifstream fin;
  stringstream ss;
  string line, chrn, tmp1, tmp2;
  int aluBegin, aluEnd, cii, flag;
  map < pair<int, int>, map<string, string> > pos_pnProb;
  map < pair<int, int>, map<string, string> > pos_pnProb_all;
  map < pair<int, int>, map<string, string> >::iterator ppp;
  map<string, string>::iterator pp;
  vector <string>::iterator ci, pi;
  //seqan::VcfStream vcfout("-", seqan::VcfStream::WRITE);
  // seqan::VcfStream vcfout(seqan::toCString(f_out), seqan::VcfStream::WRITE);
  seqan::VcfStream vcfout(f_out.c_str(), seqan::VcfStream::WRITE); //BJARNI INS
  for ( ci = chrns.begin(); ci != chrns.end(); ci++)
    appendValue(vcfout.header.sequenceNames, *ci);
  for ( pi = pns.begin(); pi != pns.end(); pi++) 
    appendValue(vcfout.header.sampleNames, *pi);
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("fileformat", "VCFv4.1"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("fileDate", "201402"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("reference", "hg18"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FILTER", "<ID=PL,Number=3,Type=Integer, Description=\"Phred-scaled likelihoods for genotypes\">"));
  seqan::VcfRecord record;    
  record.ref = ".";
  record.alt = "1";
  record.qual = 0;
  record.filter = ".";
  record.info = ".";
  record.format = ".";
  float p0, p1, p2;      
  string tmpfield;
  for ( cii = 0, ci = chrns.begin(); ci != chrns.end(); ci++, cii++) {
    pos_pnProb.clear();
    pos_pnProb_all.clear();
    for ( pi = pns.begin(); pi != pns.end(); pi++) {
      fin.open( (path0 + *pi + f_in_suffix).c_str() );
      assert(fin);
      getline(fin, line); 
      flag = 1;
      while ( getline(fin, line) ) {
	ss.clear(); ss.str( line );
	ss >> chrn >>  aluBegin >> aluEnd;
	for (int ti = 0; ti < col_00 - 4; ti++) ss >> tmpfield;
	ss >> p0 >> p1 >> p2 ;
	if (chrn != *ci) {
	  if (flag) continue;
	  else break;
	  }
	flag = 0;
	string phred_scaled_str = phred_scaled(p0, p1, p2);
	pos_pnProb_all[ make_pair(aluBegin, aluEnd) ][*pi] = phred_scaled_str;
	if (p1 > p0 or p2 > p0) pos_pnProb[ make_pair(aluBegin, aluEnd) ][*pi] = phred_scaled_str;
      }
      fin.close();
    }
    cout << *ci << " size " << pos_pnProb.size() << endl;
    // Write out the records.
    for (ppp = pos_pnProb.begin(); ppp != pos_pnProb.end(); ppp++) {
      record.beginPos = (ppp->first).first;
      record.rID = cii;
      record.id = int_to_string((ppp->first).second);
      int n_pn = 0;
      for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
	if ( (pp = ppp->second.find(*pi)) != ppp->second.end() ) { 
	  appendValue(record.genotypeInfos, pp->second);
	  n_pn++ ;
	} else if ((pp = pos_pnProb_all[ppp->first].find(*pi)) !=  pos_pnProb_all[ppp->first].end() ) {
	  appendValue(record.genotypeInfos,  pos_pnProb_all[ppp->first][*pi]); 	  
	} else {
	  appendValue(record.genotypeInfos, "0,255,255"); // otherwise vcf consider it as missing
	}
      }
      if (n_pn > 1) writeRecord(vcfout, record);  // don't output private loci
      clear(record.genotypeInfos);
    }        
    cout << "done with " << *ci << endl;
  }  // chrn finished
  clear(record);
  seqan::close(vcfout);
  return true;
}

