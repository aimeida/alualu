#define SEQAN_HAS_ZLIB 1
#include "ins_del.h"

void RecordInfo::delete_record_list(list <RecordInfo *> & records){
  for (list <RecordInfo *>::iterator ri = records.begin(); ri != records.end(); ri++) 
    delete *ri;
  records.clear();
}

bool RecordInfo::sort_faPos (const RecordInfo* a, const RecordInfo* b) {
  return (a->pairChrn == b->pairChrn) ? (a->pairBegin < b->pairBegin) : (a->pairChrn < b->pairChrn) ;
}

bool RecordInfo::sort_insertPos (const RecordInfo* a, const RecordInfo* b) {
  return (a->insertBegin == b->insertBegin) ? (a->thisBegin < b->thisBegin) : (a->insertBegin < b->insertBegin);
}

void RecordInfo::debugprint_thisBegin(list <RecordInfo *> & records, ostream & os, size_t n) {
  size_t ni = 0;
  size_t n_max = ( (n == 0) ? records.size() : min(n, records.size()) );
  for (list <RecordInfo *>::iterator ri = records.begin(); ri != records.end(),ni < n_max; ri++,ni++) 
    os << (*ri)->thisBegin << " ";
  os << endl;
}

AluconsHandler::AluconsHandler(string fn_fa, string sn):
  FastaFileHandler(fn_fa) {
  update_seq_name(sn);
}

void AluconsHandler::update_seq_name(string sn) {
  seq_name = sn;
  seqs.clear();  
  unsigned idx = 0;
  seqan::CharString seq;
  assert (getIdByName(faiIndex, seq_name, idx));
  readSequence(seq, faiIndex, idx);
  seqan::toUpper(seq);    
  seq_len = length(seq);
  seqs[1] = seq;
  seqan::ModifiedString <seqan::CharString, seqan::ModReverse > myRev(seq);
  seqan::CharString rev_seq = myRev;
  seqs[4] = rev_seq;  // 1,4 are reversed
  reverseComplement(seq);
  seqs[2] = seq;
  rev_seq = myRev;
  seqs[3] = rev_seq; // 2,3 are reversed
}

seqan::CharString AluconsHandler::fetch_alucons(int key) {
  assert ( key >=1 and key <= 4) ;
  return seqs[key];
} 

bool align_alu_cons(seqan::CharString &ref_fa, seqan::CharString & alucons, float & sim_rate,
		    int & align_consBegin, int & align_len, float sim_th){
  TAlign align;
  seqan::Score<int> scoringScheme(0, -1, -1, -2); 
  resize(rows(align), 2);
  assignSource(row(align,0), ref_fa);  // 2,3
  assignSource(row(align,1), alucons);   // 1,4, free gap at end

  globalAlignment(align, scoringScheme, seqan::AlignConfig<true, false, false, true>());
  align_consBegin = toViewPosition(row(align, 0), 0) - toViewPosition(row(align, 1), 0);
  int align_start = max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  int align_end = min(toViewPosition(row(align, 0), length(ref_fa)), toViewPosition(row(align, 1), length(alucons)));
  /////align_end = min(toViewPosition(row(align, 0), ref_fa.size()), toViewPosition(row(align, 1), alucons.size()));
  sim_rate = 0;
  align_len = align_end - align_start;
  if ( align_len <= CLIP_BP or align_len <= sim_th * length(ref_fa) )
    return false;
  TRow &row0 = row(align,0);
  TRowIterator it0 = begin(row0);
  TRow &row1 = row(align,1);
  TRowIterator it1 = begin(row1);
  int i = 0, dif = 0;
  while ( i++ < align_start ) {  it0++; it1++; }
  while ( i++ <= align_end) {
    if ( (*it0) != (*it1) ) dif++;     ////if(isGap(it1))
    it0++; it1++;
  }
  sim_rate = 1 - dif / (float) align_len;
  return  sim_rate >= sim_th;
}

int align_alu_cons_call(seqan::CharString &ref_fa, AluconsHandler *alucons_fh, float & sim_rate, 
			int & align_consBegin, int & align_len, float sim_th){
  for (int k = 1; k <= 4; k++) {
    seqan::CharString alucons = alucons_fh->fetch_alucons(k);
    /// fixme: static_cast<seqan::CharString &> (alucons_fh->fetch_alucons(k)) 
    if (align_alu_cons(ref_fa, alucons, sim_rate, align_consBegin, align_len, sim_th))
      return k;
  }
  return 0;
}

void filter_aluRead(list <RecordInfo *> & records, AluconsHandler *alucons_fh, string file_fa, map <seqan::CharString, qNameInfo > & qName_info) {
  records.sort( RecordInfo::sort_faPos );      
  FastaFileHandler *fasta_fh;
  string chrn = "";
  seqan::CharString ref_fa;     
  list <RecordInfo *>::iterator ri = records.begin();
  while (ri != records.end()){
    map <seqan::CharString, qNameInfo >::iterator qItr = qName_info.find( (*ri)->qName);
    assert (qItr != qName_info.end() );
    if ( (qItr->second).itype != aluRead ) { 
      ri++ ;
      continue;
    }
    if ( (*ri)->pairChrn != chrn ) {
      if ( chrn != "") delete fasta_fh;
      chrn = (*ri)->pairChrn;
      fasta_fh = new FastaFileHandler(file_fa + chrn + ".fa", chrn);      
    }
    fasta_fh -> fetch_fasta_upper((*ri)->pairBegin, (*ri)->pairEnd, ref_fa);
    float sim_rate;    
    int align_consBegin, align_len;    
    int cons_ori = align_alu_cons_call(ref_fa, alucons_fh, sim_rate, align_consBegin, align_len, ALUCONS_SIMILARITY);

    //cout << "check1 " << records.size() << " " << (*ri)->pairChrn << endl;    
    //cout << (*ri)->insertBegin << " " << chrn << " " << sim_rate << endl;
    int flag = 0;
    int insertLen;
    if ( cons_ori > 0 ) {
      if ( (*ri)->thisBegin <= (*ri)->insertBegin - MIN_MATCH_LEN) { // this read is left read
	if ( (*ri)->thisRC == (*ri)->pairRC) {
	  (qItr->second).insertAluOri = 5 - cons_ori;
	  insertLen = alucons_fh->seq_len - align_consBegin;
	} else { 
	  (qItr->second).insertAluOri = cons_ori;
	  insertLen = align_consBegin + align_len;
	}
	insertLen += ( (*ri)->insertBegin - (*ri)->thisBegin );	
	flag = 1;
      } else if ( (*ri)->thisEnd >= (*ri)->insertEnd + MIN_MATCH_LEN) { // this read is right read
	if ( (*ri)->thisBegin <= (*ri)->insertBegin - MIN_MATCH_LEN) { // this read is left read
	  if ( (*ri)->thisRC == (*ri)->pairRC) {
	    (qItr->second).insertAluOri = 5 - cons_ori;
	    insertLen = align_consBegin + align_len;
	  } else { 
	    (qItr->second).insertAluOri = cons_ori;
	    insertLen = alucons_fh->seq_len - align_consBegin;
	  }
	  insertLen += ( (*ri)->thisEnd - (*ri)->insertEnd );
	  flag = 1;
	}
      }
    }
    if (!flag) {
      qName_info.erase( (*ri)->qName );
      delete *ri;
      ri = records.erase(ri);
    } else {
      (*ri)->insertLen = insertLen;
      ri++;
    }
  }
  delete fasta_fh; 
}

int major_element(vector <int> & ints, float & major_ratio) {
  map <int, int>  key_cnt;
  for (vector <int>::iterator it = ints.begin(); it != ints.end(); it++)
    addKey(key_cnt, *it);
  int max_key = (key_cnt.begin())->first;
  int max_val = (key_cnt.begin())->second;
  for (map <int, int>::iterator ki = key_cnt.begin(); ki != key_cnt.end(); ki++) {
    if ( ki->second > max_val ) {
      max_key = ki->first;
      max_val = ki->second;
    }
  }
  major_ratio = max_val / (float)(ints.size());
  return max_key;
}

void filter_aluRead_aluOri(list <RecordInfo *> & records, map <seqan::CharString, qNameInfo > & qName_info) {  
  records.sort( RecordInfo::sort_insertPos );      
  map <seqan::CharString, qNameInfo >::iterator qItr ;
  map <int, vector <int> > insertBegin_ori;
  list <RecordInfo *>::iterator ri ;
  for ( ri = records.begin(); ri != records.end(); ri++ ) {
    assert( (qItr = qName_info.find( (*ri)->qName)) != qName_info.end() );
    if ( (qItr->second).itype != aluRead ) continue;
    insertBegin_ori[(*ri)->insertBegin].push_back((qItr->second).insertAluOri);
  }  

  map <int, vector <int> > :: iterator oi;
  map <int, int > insertBegin_majorOri;
  for ( oi = insertBegin_ori.begin(); oi != insertBegin_ori.end(); oi++) {
    float major_freq;
    int major_key = major_element(oi->second, major_freq);
    if (major_freq < 0.9) {
      //cout << major_freq << " " ;
      insertBegin_majorOri[oi->first] = 0; // majorOri valid value: 1,2,3,4
    } else
      insertBegin_majorOri[oi->first] = major_key;
  }
  
  ri = records.begin();
  while ( ri != records.end() ) {
    assert( (qItr = qName_info.find( (*ri)->qName)) != qName_info.end() );
    if ( (qItr->second).itype != aluRead ) {
      ri++;
      continue;
    }
    if ( (qItr->second).insertAluOri != insertBegin_majorOri[(*ri)->insertBegin] ) {
      qName_info.erase( (*ri)->qName );
      delete *ri;
      ri = records.erase(ri);
    } else ri++;
  }  
}

void filter_pos_by_cnt(list <RecordInfo *> & records, map <seqan::CharString, qNameInfo > & qName_info, int num_informative_read) {
  records.sort( RecordInfo::sort_insertPos );      
  map <int, int> insertBegin_cnt;
  set <int> insertPos_all; 
  set <int> insertPos_rm; 
  set <int> insertPos_keep; 
  list <RecordInfo *>::iterator ri ;
  for ( ri = records.begin(); ri != records.end(); ri++ ) {
    insertPos_all.insert( (*ri)->insertBegin );
    map <seqan::CharString, qNameInfo >::iterator qItr = qName_info.find( (*ri)->qName);
    if ( (qItr->second).itype == aluRead or (qItr->second).itype == aluClipRead )  
      addKey(insertBegin_cnt, (*ri)->insertBegin);
  }  
  for (map <int, int>::iterator ii = insertBegin_cnt.begin(); ii != insertBegin_cnt.end(); ii++) {
    if ( ii->second >= num_informative_read ) 
      insertPos_keep.insert(ii->first);
    //cout << "insertBegin " << ii->first << " " << ii->second << endl;
  }
  for (set <int>::iterator ip = insertPos_all.begin(); ip != insertPos_all.end(); ip++ ) {
    if ( insertPos_keep.find(*ip) == insertPos_keep.end())
      insertPos_rm.insert(*ip);
  }
  
  ri = records.begin();
  while (ri != records.end()) {
    if (insertPos_rm.find((*ri)->insertBegin) != insertPos_rm.end() ) {
      qName_info.erase( (*ri)->qName );
      delete *ri;
      ri = records.erase(ri);
    } else ri++;
  }
}

void filter_aluClipRead(list <RecordInfo *> & records, map <seqan::CharString, qNameInfo > & qName_info, AluconsHandler *alucons_fh, string file_fa, map <int, int> & insertBegin_consLen) {
  // if aluClipRead, find cons length
  // if not aluClipRead, assign as unknowRead, use insert_len
  cout << "fixme: filter_aluClipRead() to be implemented \n";
}

void filter_aluSkipRead(list <RecordInfo *> & records, map <seqan::CharString, qNameInfo > & qName_info, AluconsHandler *alucons_fh, string file_fa) {
  // if aluClipRead, find cons length
  // if not assign as unknowRead, use insert_len
  cout << "fixme: filter_aluSkipRead() to be implemented \n";
}

bool writeReads_type(RecordInfo * ii, int alu_ins_len, map <seqan::CharString, qNameInfo > & qName_info, int & midCnt, int &clipCnt, int &unknowCnt, stringstream & ss) {
  map <seqan::CharString, qNameInfo >::iterator qItr = qName_info.find( ii->qName );
  assert (qItr != qName_info.end() );
  if ( (qItr->second).itype == aluSkipRead ) { 
    clipCnt++;
    return true;
  } 
  if ( (qItr->second).itype == aluRead or (qItr->second).itype == aluClipRead) { 
    midCnt++;
    return true;
  } 
  if ( (qItr->second).itype == unknowRead ) { 
    unknowCnt++;
    //cout << "alu_ins_len " << alu_ins_len << " " <<  ii->insertLen << endl;
    ss << " " << ii->rgIdx << ":" << abs(ii->insertLen) + alu_ins_len;
    return true;
  }
  return false;
}

/*
bool writeReads_insertLen(RecordInfo * ii, int alu_ins_len, map <seqan::CharString, qNameInfo > & qName_info, int & midCnt, int &clipCnt, int &unknowCnt, stringstream & ss) {
  map <seqan::CharString, qNameInfo >::iterator qItr = qName_info.find( ii->qName );
  assert (qItr != qName_info.end() );
  cout << "writeReads_insertLen() to be implemented\n";
  return false;
}
*/

bool align_aluSkip(seqan::CharString & seq_read, seqan::CharString & seq_ref, int readLen, float th_score = 0.8, float th_sim = 0.9) {
  seqan::Score<int> scoringScheme(1, -1, -1, -2); // match, mismatch, gap extend, gap open
  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align,0), seq_read); // 2,3, true means free gap
  assignSource(row(align,1), seq_ref);   // 1,4
  int score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>()); 
  int align_start = max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  int align_end = min(toViewPosition(row(align, 0), length(seq_read)), toViewPosition(row(align, 1), length(seq_ref)));
  bool align_pass = (score >= round(th_score * readLen)) and align_end - align_start >= th_sim * readLen ;
  /*
    if (!align_pass) {
    cout << score << " " << align_end - align_start << " " << readLen << endl ;
    cout << seq_read << endl;
    cout << seq_ref << endl;
    cout << align << endl;
    }
  */
  return align_pass;
}


I_READ classify_read(seqan::BamAlignmentRecord &record, seqan::CharString & ref_fa, int insertBegin, int insertEnd, int insertFlank) {
  if (record.rID != record.rNextId or abs(record.tLen) >= DISCORDANT_LEN) 
    return aluRead;
  if ( ! (record.rID == record.rNextId and hasFlagRC(record) != hasFlagNextRC(record) and
	  abs(record.tLen) < DISCORDANT_LEN ) )
    return uselessRead;
  int thisEnd = record.beginPos + (int)getAlignmentLengthInRef(record);   
  if ( ( has_soft_first(record, CLIP_BP) and abs( record.beginPos - insertEnd) < 30 ) or  	    
       ( has_soft_last(record, CLIP_BP) and abs( thisEnd - insertBegin) < 30 ) )
    return aluClipRead;
  int readLen = length(record.seq);
  if ( ( left_read(record) and record.beginPos >= insertBegin - readLen and record.beginPos < insertEnd + CLIP_BP ) or 
       ( !left_read(record) and thisEnd > insertBegin - CLIP_BP and thisEnd <= insertEnd + readLen) ) {
    int non_match_bp = count_non_match(record);
    int ref_fa_begin = record.beginPos - ( insertBegin - insertFlank);
    seqan::CharString ref_fa_toalign = infix(ref_fa, ref_fa_begin - non_match_bp, ref_fa_begin + readLen + non_match_bp);
    if (  align_aluSkip(record.seq, ref_fa_toalign, readLen) )  
      return aluSkipRead;
    return unknowRead;
  }
  return uselessRead;
}

void insert_flanking_reads(string chrn, FastaFileHandler * fasta_fh, vector< pair<int, int> > & insert_pos, map<string, int> & rg_to_idx,  list <RecordInfo *> & aluinfo, BamFileHandler* bam_fh, map <seqan::CharString, qNameInfo > & qName_info, string path_cons ){
  ofstream fout;
  if (path_cons !="") fout.open(path_cons.c_str());    
  seqan::BamAlignmentRecord record;
  string read_status;
  for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
    int region_begin = (*pi).first - ALU_INSERT_FLANK;
    int region_end = (*pi).second + ALU_INSERT_FLANK;	
    if (! bam_fh->jump_to_region(chrn, region_begin, region_end)) 
      continue;    
    seqan::CharString ref_fa;
    fasta_fh->fetch_fasta_upper(region_begin, region_end, ref_fa);
    while ( true ) {
      read_status = bam_fh->fetch_a_read(chrn, region_begin, region_end, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_insert_read(record)) continue;
      string pairChrn;
      if ( (bam_fh->get_chrn)(record.rNextId, pairChrn) == false) continue;	       
      map <seqan::CharString, qNameInfo >::iterator qItr = qName_info.find(record.qName);
      if ( qItr != qName_info.end() and (qItr->second).itype != unknowRead) 
	continue;

      I_READ iread = classify_read(record, ref_fa, (*pi).first, (*pi).second, ALU_INSERT_FLANK); //aluRead needs ReCheck afterwards 
      if ( iread == uselessRead) continue;      
      if ( qItr == qName_info.end() ) { // not exists
	qNameInfo  one_qi = qNameInfo(iread, record.seq);
	qName_info.insert( pair<seqan::CharString, qNameInfo> (record.qName, one_qi) );
	int rgIdx = get_rgIdx(rg_to_idx, record);
	int thisEnd = record.beginPos + (int)getAlignmentLengthInRef(record);   
	RecordInfo * one_ri = new RecordInfo(pairChrn, record.qName, (*pi).first, (*pi).second, record.beginPos,
					     thisEnd, record.pNext, record.pNext + length(record.seq), rgIdx,
					     record.tLen, hasFlagRC(record), hasFlagNextRC(record));	
	if ( iread == aluRead) 	one_ri->insertLen = 0;  
	aluinfo.push_back(one_ri);
      } else  if ( (qItr->second).itype == unknowRead and iread != unknowRead) {
	(qItr->second).itype = iread;
	(qItr->second).seq = record.seq;
      }      
    }
    cout << "read " << (*pi).first << endl;
  }
  if (path_cons !="")  fout.close();
}

void write_cnt0(string chrn, list <RecordInfo *> & records, map <seqan::CharString, qNameInfo > & qName_info,map <int, int> & insertBegin_consLen, int alu_ins_len_fixed, ofstream &fout, int num_informative_read) {
  bool use_fixed_consLen = insertBegin_consLen.empty() ? true : false;
  int alu_ins_len;
  records.sort( RecordInfo::sort_insertPos );      
  list <RecordInfo *>::iterator ri;  
  int midCnt, clipCnt, unknowCnt;
  int insertBegin = 0, insertEnd;
  stringstream ss;
  for (list <RecordInfo *>::iterator ri = records.begin(); ri != records.end(); ri++) {
    if ( (*ri)->insertBegin != insertBegin ) {
      if (insertBegin != 0) 
	fout << chrn << " " << (*ri)->insertBegin << " " << (*ri)->insertEnd << " " << alu_ins_len 
	     << " " << midCnt << " " << clipCnt << " " << unknowCnt << " " << ss.str() << endl;
      insertBegin = (*ri)->insertBegin;
      alu_ins_len = use_fixed_consLen ? alu_ins_len_fixed : insertBegin_consLen[insertBegin];
      insertEnd = (*ri)->insertEnd;
      midCnt = 0, clipCnt = 0, unknowCnt = 0;
      ss.clear(); ss.str("");
    }
    writeReads_type(*ri, alu_ins_len, qName_info, midCnt, clipCnt, unknowCnt, ss);  
    //writeReads_insertLen(*ri, alu_ins_len, qName_info, midCnt, clipCnt, unknowCnt, ss);  
  }
  fout << chrn << " " << insertBegin << " " << insertEnd << " -1 "
       << midCnt << " " << clipCnt << " " << unknowCnt << ss.str() << endl;
}

bool parseline_del_tmp1(string &line, string & output_line, map <int, EmpiricalPdf *> & pdf_rg){
  float *log10_gp = new float[3];
  stringstream ss, ss_out;
  string chrn, insertBegin, insertEnd;
  int estiAluLen, midCnt, clipCnt, unknowCnt;
  ss.clear(); ss.str(line); 
  ss >> chrn >> insertBegin >> insertEnd >> estiAluLen >> midCnt >> clipCnt >> unknowCnt ;
  float prob_ub = pow(10, -LOG10_RATIO_UB);
  float prob_known = (midCnt+clipCnt)/(float)(midCnt + clipCnt + unknowCnt);
  for (int i = 0; i < 3; i++) log10_gp[i] = 0;
  if (prob_known) {
    log10_gp[0] = clipCnt * log10 ( prob_known * prob_ub ) + midCnt * log10 ( prob_known * (1 - prob_ub) );
    log10_gp[1] = (midCnt + clipCnt) * log10 (prob_known * 0.5) ; 
    log10_gp[2] = midCnt * log10 ( prob_known * prob_ub ) + clipCnt * log10 ( prob_known * (1 - prob_ub) );
  }
  if (unknowCnt) { 
    int insert_len, idx;
    string token;
    for (int i = 0; i < unknowCnt; i++) {
      getline(ss, token, ':');
      seqan::lexicalCast2(idx, token);
      getline(ss, token, ' ');
      seqan::lexicalCast2(insert_len, token);      
      float p_y = pdf_rg[idx]->pdf_obs(insert_len);
      float p_z = pdf_rg[idx]->pdf_obs(insert_len - estiAluLen);
      //float freq0 = 0.67;  // high FP ratio      
      float freq0 = ( midCnt + 1 )/(float)(midCnt + clipCnt + 2); // 1 and 2 are psudo count
      log10_gp[0] += log10 (p_y * (1 - prob_known));
      log10_gp[1] += log10 ((freq0 * p_y + (1 - freq0) * p_z) * (1 - prob_known) ) ;
      log10_gp[2] += log10 (p_z * (1 - prob_known));
    }
  }
  bool use_this_line = false;
  if ( !p11_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    ss_out << chrn << " " << insertBegin << " " << insertEnd << " " << midCnt << " " << clipCnt << " " << unknowCnt ;
    float *gp = new float[3];
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
    ss_out << " " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2]; 
    delete gp;    
    output_line = ss_out.str();
    use_this_line = true;
  }
  delete log10_gp;
  return use_this_line;
}

int main( int argc, char* argv[] )
 {
   string config_file = argv[1];
   string opt =  argv[2];
   if (argc < 3) exit(1);
   
   boost::timer clocki;    
   clocki.restart();

   ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
   map<int, string> ID_pn;
   get_pn(cf_fh.get_conf("file_pn"), ID_pn);      
   vector<string> chrns;  
   for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
   chrns.push_back("chrX");
   chrns.push_back("chrY");  
   string file_dist_prefix = cf_fh.get_conf("file_dist_prefix");

   if (opt == "write_tmp1") { 
     int idx_pn;
     seqan::lexicalCast2(idx_pn, argv[3]);
     assert (argc == 4);
     string pn = ID_pn[idx_pn];      
     string bam_input = cf_fh.get_conf("file_bam_prefix") + pn + ".bam";

#ifdef DEBUG_MODE   
     bam_input = cf_fh.get_conf("file_test_bam_prefix") + pn + ".bam";
#endif

     string bai_input = bam_input + ".bai";  
     BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bai_input); 

#ifdef DEBUG_MODE
     cerr << "!!!! DEBUG_MODE, using test bam input and chr1 \n";
     chrns.clear();
     chrns.push_back("chr1");  
#endif       

     string file_rg = get_name_rg(file_dist_prefix, pn);
     map<string, int> rg_to_idx;
     parse_reading_group( file_rg, rg_to_idx );
     AluconsHandler *alucons_fh = new AluconsHandler(cf_fh.get_conf("file_alu_cons"), cf_fh.get_conf("type_alu_cons"));
     string file_fa = cf_fh.get_conf("file_fa_prefix");         
     string path0 = cf_fh.get_conf("file_ins_del0");
     check_folder_exists( path0);
     string path_cons = cf_fh.get_conf("file_ins_cons");  // for consensus sequence
     check_folder_exists( path_cons);


     ofstream fout1 ( (path0 + pn + ".tmp1").c_str());     
     fout1 << "chr insertBegin insertEnd estimatedAluLen midCnt clipCnt unknowCnt unknowStr\n";    
     for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
       string file_cons = path_cons + *ci + "_" + pn;
       FastaFileHandler * fasta_fh = new FastaFileHandler(cf_fh.get_conf("file_fa_prefix") + *ci + ".fa", *ci);
       vector< pair<int, int> > insert_pos;
       string posfile = cf_fh.get_conf("file_clip_reads") + *ci + "_" + cf_fh.get_conf("freq_min") + "_" + cf_fh.get_conf("freq_max");
       read_first2col(posfile, insert_pos);
       list <RecordInfo *> aluinfo; 
       map <seqan::CharString, qNameInfo > qName_info;
       insert_flanking_reads(*ci, fasta_fh, insert_pos, rg_to_idx, aluinfo, bam_fh, qName_info, file_cons);
       
//       filter_aluRead(aluinfo, alucons_fh, file_fa, qName_info);  // for aluRead
//       filter_aluRead_aluOri(aluinfo, qName_info);
//       filter_pos_by_cnt(aluinfo, qName_info, MIN_READS_CNT);
//       //cout << "size " << aluinfo.size() << endl;
//       map <int, int> insertBegin_consLen;
//       //filter_aluClipRead(aluinfo, qName_info, alucons_fh, file_fa, insertBegin_consLen); 
//       //filter_aluSkipRead(aluinfo, qName_info, alucons_fh, file_fa); // re_align
//       write_cnt0(*ci, aluinfo, qName_info, insertBegin_consLen, alucons_fh->seq_len, fout1, MIN_READS_CNT);       
       qName_info.clear();
       delete fasta_fh;
       RecordInfo::delete_record_list(aluinfo);       
     }
     fout1.close();
     delete bam_fh;    
     delete alucons_fh;

   } else if (opt == "write_tmp2") { 
     int idx_pn = seqan::lexicalCast <int> (argv[3]);
     string pn = ID_pn[idx_pn];      
     map <int, EmpiricalPdf *> pdf_rg;
     string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5                                   
     read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
     string path0 = cf_fh.get_conf("file_ins_del0");
     string fn_tmp1 = path0 + pn + ".tmp1";
     string fn_tmp2 = path0 + pn + ".tmp2";

     ofstream fout(fn_tmp2.c_str());
     fout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";
     string line, output_line;
     ifstream fin(fn_tmp1.c_str());
     assert(fin);
     getline(fin, line); // read header                                                               
     while (getline(fin, line))
       if (parseline_del_tmp1(line, output_line, pdf_rg))
	 fout << output_line << endl;
     fin.close();
     fout.close();
     
   } else if ( opt == "debug") {
     string ref_fa = "CAGCACTTTGGGAGGCCGAGGCGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAATACGGTGAAACCCCGTCTCTACTAAAAATACAAAAA";
     //seqan::String<seqan::Dna> ref_fa = "CAGCACTTTGGGAGGCCGAGGCGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAATACGGTGAAACCCCGTCTCTACTAAAAATACAAAAA";     
   } else if ( opt == "write_testbam") { // create tmp bam file, for testing 
     int idx_pn;
     seqan::lexicalCast2(idx_pn, argv[3]);
     string pn = ID_pn[idx_pn];      
     string bam_input = cf_fh.get_conf("file_bam_prefix") + pn + ".bam";
     string bai_input = bam_input + ".bai";  
     string bam_output = cf_fh.get_conf("file_test_bam_prefix") + pn + ".bam";
     BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bai_input, bam_output);      
     vector< pair<int, int> > insert_pos;
     seqan::BamAlignmentRecord record;
     string read_status;
     chrns.clear();
     chrns.push_back("chr1");
     
     for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
       // posfile manually created, avoid one read appear multiple times if pos are close to each other
       string posfile = cf_fh.get_conf("file_clip_reads") + *ci + "_" + cf_fh.get_conf("freq_min") + "_" + cf_fh.get_conf("freq_max") + ".testbam";;
       read_first2col(posfile, insert_pos);
       for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
	 int region_begin = (*pi).first - ALU_INSERT_FLANK;
	 int region_end = (*pi).second + ALU_INSERT_FLANK;	
	 if (! bam_fh->jump_to_region(*ci, region_begin, region_end)) 
	   continue;
	 int ni = 0;
	 while ( true ) {
	   read_status = bam_fh->fetch_a_read(*ci, region_begin, region_end, record);
	   if (read_status == "stop" ) break;
	   if (read_status == "skip") continue;
	   string pairChrn;
	   if ( (bam_fh->get_chrn)(record.rNextId, pairChrn) == false) continue;	    
	   if (! (bam_fh->write_a_read(record) )) cout << "err1\n";
	   ni ++;
	 }  
	 cout  << "checked " << region_begin << " to " << region_end << " with " << ni << " counts\n";
       }
     }
     delete bam_fh;
   }
   
   cout << "total time used " << clocki.elapsed() << endl;
   clocki.restart();  
   return 0;
 }

