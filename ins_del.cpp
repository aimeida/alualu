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

ConfigFileHandler::ConfigFileHandler(string config_file) {
  ifstream fin(config_file.c_str());
  if (!fin) 
    try {
      throw 1;    
    } catch(int e) {
      cerr << "#### ERROR #### file: "<< config_file << " not exists!!!" << endl;
    }   
  stringstream ss;
  string line, key, value;
  while (getline(fin, line)) {
    if (line[0] == '#')  continue;
    ss.clear(); ss.str( line );  
    ss >> key >> value;
    configs[key] = value;
  }
  fin.close();
}

string ConfigFileHandler::get_cf_val(string key) {
  map <string, string>::iterator ci;
  if ( (ci = configs.find(key)) == configs.end() ) {
    cerr << "#### ERROR #### key: " << key <<" doesn't exist\n";
    exit(1);
  }
  return ci->second;
}

BamFileHandler::BamFileHandler(vector<string> &chrns, string bam_input, string bai_input, string bam_output) 
  : fn_bai(bai_input)
  , fn_bam_out(bam_output)
  , nameStoreCache(nameStore)
  , context(nameStore, nameStoreCache)
{
  if (fn_bai != "" and read(baiIndex, fn_bai.c_str()) != 0){
    cerr << "ERROR: Could not read BAI index file " << fn_bai << endl;
    exit(1);
  }
  if (!open(inStream, bam_input.c_str(), "r")) {
    std::cerr << "ERROR: Could not open " << bam_input << " for reading.\n";
    exit(1);
  }
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );
  if (!fn_bam_out.empty()) {
    assert(open(outStream, bam_output.c_str(), "w") );
    assert(write2(outStream, header, context, seqan::Bam()) == 0);
  }
  int rID;
  for ( vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    assert ( (*ci).substr(0,3) == "chr");
    if ( !getIdByName(nameStore, *ci, rID, nameStoreCache) )
      if (!getIdByName(nameStore, (*ci).substr(3), rID, nameStoreCache)) {
	cerr << "ERROR: Reference sequence named "<< *ci << " not known.\n";
	exit(1);
      }
    //cout << *ci << " " << rID  << endl;
    chrn_rID[*ci] = rID;
    rID_chrn[rID] = *ci;
  }    
} 

BamFileHandler::~BamFileHandler(void){
  seqan::close(inStream); 
  if (!fn_bam_out.empty())
    seqan::close(outStream); 
}

bool BamFileHandler::jump_to_region(string chrn, int region_begin, int region_end) {
  if ( fn_bai.empty())  // bai file not available
    return false;
  bool hasAlignments = false;	
  if (!jumpToRegion(inStream, hasAlignments, context, chrn_rID[chrn], region_begin, region_end, baiIndex))
    return false;
  return hasAlignments;
}

string BamFileHandler::fetch_a_read(string chrn, int region_begin, int region_end, seqan::BamAlignmentRecord & record) {
  if (atEnd(inStream)) return "stop";
  assert (!readRecord(record, context, inStream, seqan::Bam())); 
  if (record.rID != chrn_rID[chrn] || record.beginPos >= region_end) return "stop";
  if (record.beginPos + (int)getAlignmentLengthInRef(record) < region_begin) return "skip";
  if (!QC_insert_read(record) ) return "skip";
  return "record";
}

bool BamFileHandler::get_chrn(int query_rid, string & pairChrn) {
  map <int, string>::iterator ri;
  if ( (ri = rID_chrn.find(query_rid)) == rID_chrn.end() )
    return false;
  pairChrn = ri->second;
  return true;
}

bool BamFileHandler::write_a_read(seqan::BamAlignmentRecord & record) {
  if (fn_bam_out.empty())  return false;
  return write2(outStream, record, context, seqan::Bam()) == 0;
}

FastaFileHandler::FastaFileHandler(string fn_fa) {
  string fn_fai = fn_fa + ".fai";
  if ( read(faiIndex, fn_fa.c_str()) ) {
    build(faiIndex, fn_fa.c_str() ); 
    write(faiIndex, fn_fai.c_str() ); 
  }      
}

void FastaFileHandler::fetch_fasta_upper(string seq_name, int beginPos, int endPos, seqan::CharString &seq){
  unsigned idx = 0;
  assert (getIdByName(faiIndex, seq_name, idx));
  //cout << beginPos << " " << endPos << endl;
  assert (!readRegion(seq, faiIndex, idx, beginPos, endPos));
  seqan::toUpper(seq);
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

void read_first2col(string fn, vector < pair<int, int> > & insert_pos) {
  insert_pos.clear();
  ifstream fin(fn.c_str());
  assert(fin);
  stringstream ss;
  string line;
  getline(fin, line); // skip header!
  int beginPos, endPos;
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );  
    ss >> beginPos >> endPos;
    insert_pos.push_back( make_pair(beginPos, endPos) );
#ifdef DEBUG_MODE
    if (insert_pos.size() > 10) break;
#endif
  }
  fin.close();
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

void pair_is_alu(list <RecordInfo *> & records, AluconsHandler *alucons_fh, string file_fa) {
  records.sort( RecordInfo::sort_faPos );      
  FastaFileHandler *fasta_fh;
  string chrn = "";
  seqan::CharString ref_fa;     
  list <RecordInfo *>::iterator ri = records.begin();
  while (ri != records.end()){
    if ( (*ri) -> maybe_clip_read() ) {
      ri++ ;
      continue;
    }
    //cout << "debug " <<  (*ri)->insertBegin << " " << (*ri)->seq << " " << (*ri)->insertLen << endl;
    if ( chrn == "" or (*ri)->pairChrn != chrn ) {
      if ( chrn != "") delete fasta_fh;
      chrn = (*ri)->pairChrn;
      fasta_fh = new FastaFileHandler(file_fa + chrn + ".fa");      
    }
    fasta_fh -> fetch_fasta_upper(chrn, (*ri)->pairBegin, (*ri)->pairEnd, ref_fa);
    float sim_rate;    
    int align_consBegin, align_len;    
    int cons_ori = align_alu_cons_call(ref_fa, alucons_fh, sim_rate, align_consBegin, align_len, ALUCONS_SIMILARITY);

    //cout << "check1 " << records.size() << " " << (*ri)->pairChrn << endl;    
    //cout << (*ri)->insertBegin << " " << chrn << " " << sim_rate << endl;
    int flag = 0;
    if ( cons_ori > 0 ) {
      if ( (*ri)->thisBegin <= (*ri)->insertBegin - MIN_MATCH_LEN) { // this read is left read
	if ( (*ri)->thisRC == (*ri)->pairRC) {
	  (*ri)->insertAluOri = 5 - cons_ori; 
	  (*ri)->insertLen = alucons_fh->seq_len - align_consBegin;
	} else { 
	  (*ri)->insertAluOri = cons_ori;
	  (*ri)->insertLen = align_consBegin + align_len;
	}
	(*ri)->insertLen += ( (*ri)->insertBegin - (*ri)->thisBegin );	
	flag = 1;
      } else if ( (*ri)->thisEnd >= (*ri)->insertEnd + MIN_MATCH_LEN) { // this read is right read
	if ( (*ri)->thisBegin <= (*ri)->insertBegin - MIN_MATCH_LEN) { // this read is left read
	  if ( (*ri)->thisRC == (*ri)->pairRC) {
	    (*ri)->insertAluOri = 5 - cons_ori;
	    (*ri)->insertLen = align_consBegin + align_len;
	  } else { 
	    (*ri)->insertAluOri = cons_ori;
	    (*ri)->insertLen = alucons_fh->seq_len - align_consBegin;
	  }
	  (*ri)->insertLen += ( (*ri)->thisEnd - (*ri)->insertEnd );
	  flag = 1;
	}
      }
    }
    if (!flag) {
      delete *ri;
      ri = records.erase(ri);
    } else ri++;
  }
  delete fasta_fh; 
}

void keep_pos_of_enough_reads(list <RecordInfo *> & records, int num_informative_read, map <int, int> & insertPos_cnt) {
  records.sort( RecordInfo::sort_insertPos );      
  insertPos_cnt.clear();
  set <int> insertPos_rm; // only reads @ this pos are kept
  list <RecordInfo *>::iterator ri = records.begin();
  // fixme: read tmp1, rm reads that align cons at wrong direction 
  while (ri != records.end()) {
    addKey(insertPos_cnt, (*ri)->insertBegin);
    ri++;
  }
  for (map <int, int>::iterator ii = insertPos_cnt.begin(); ii != insertPos_cnt.end(); ii++) {
    if (ii->second < num_informative_read) 
      insertPos_rm.insert(ii->first);
  }
  ri = records.begin();
  while (ri != records.end()) {
    if (insertPos_rm.find((*ri)->insertBegin) != insertPos_rm.end() ) {
      delete *ri;
      ri = records.erase(ri);
    } else ri++;
  }
  for (set <int>::iterator im = insertPos_rm.begin(); im != insertPos_rm.end(); im++) 
    insertPos_cnt.erase(*im);
}

void filter_clip_read(list <RecordInfo *> & records, AluconsHandler *alucons_fh, string file_fa) {
  // if not a ok clip_read, assign as normal read pair with known insert_len
  cout << "fixme: not implemented \n";
  
}

void write_cnt0(list <RecordInfo *> & records, string fn, int num_informative_read ) {
  records.sort( RecordInfo::sort_insertPos );      
  list <RecordInfo *>::iterator ri = records.begin();
  cout << records.size() << endl;
  while (ri != records.end()) {
    //fixme: to continue
    ri++;
  }
  ofstream fout( fn.c_str() );
  fout << "chr aluBegin aluEnd mean_coverage midCnt clipCnt unknowCnt unknowStr\n";    
  fout.close();
}

void parse_reading_group(string file_rg, map<string, int> & rg_to_idx){
  ifstream fin (file_rg.c_str() );
  int idx = 0;
  assert(fin);
  string rg;
  while (fin >> rg) rg_to_idx[rg] = idx++;
  fin.close();
  //cout << rg_to_idx.size() << " reading groups exist\n";
}

void search_chrn(string chrn, vector< pair<int, int> > & insert_pos, map<string, int> & rg_to_idx,  list <RecordInfo *> & aluinfo, BamFileHandler* bam_fh) {
  seqan::BamAlignmentRecord record;
  int region_begin, region_end;
  string read_status;
  RecordInfo * one_read;
  for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
    region_begin = (*pi).first - ALU_FLANK;
    region_end = (*pi).second + ALU_FLANK;	
    if (! bam_fh->jump_to_region(chrn, region_begin, region_end)) 
      continue;
    size_t cnt = 0;
    while ( true ) {
      read_status = bam_fh->fetch_a_read(chrn, region_begin, region_end, record);
      int flag = -1;
      if (read_status == "stop" ) break;
      cnt += 1;
      if (read_status == "skip") continue;
      string pairChrn;
      if ( (bam_fh->get_chrn)(record.rNextId, pairChrn) == false) continue;	    
      int thisEnd = record.beginPos + (int)getAlignmentLengthInRef(record); 
      if (record.rID != record.rNextId or abs(record.tLen) >= DISCORDANT_LEN) { 
	flag = 0;
      } else if (record.rID == record.rNextId and hasFlagRC(record) != hasFlagNextRC(record) and
		 abs(record.tLen) < DISCORDANT_LEN )  {  // proper mapped read pair
	if ( ( has_soft_first(record, CLIP_BP) and abs( record.beginPos - (*pi).second) < 30 ) or  	    
	     ( has_soft_last(record, CLIP_BP) and abs( thisEnd - (*pi).first) < 30 ) )
	  flag = 1;
      }
      if (flag < 0) continue;
      
      seqan::BamTagsDict tags(record.tags);
      unsigned idx_rg1;  // idx in bam file
      assert (findTagKey(idx_rg1, tags, "RG"));
      string rg_name = seqan::toCString(getTagValue(tags, idx_rg1));

      int idx_rg2 = 0; // idx in my file 
      map <string, int>::iterator rgi = rg_to_idx.find(rg_name);
      if (  rgi != rg_to_idx.end() )
	idx_rg2 = rgi->second;
      one_read = new RecordInfo(pairChrn, (*pi).first, (*pi).second, record.beginPos, thisEnd, record.pNext, 
				record.pNext + length(record.seq), idx_rg2, "", hasFlagRC(record), hasFlagNextRC(record));
      if (flag == 1) 
	one_read->seq = seqan::toCString(record.seq);
      aluinfo.push_back(one_read);	      
    }
    cout  << "read " << (*pi).first << " with " << cnt << " counts\n";
  }
}

int main( int argc, char* argv[] )
 {
   if (argc < 2) exit(1);
   string opt = argv[1];
   string config_file = argv[2];

   boost::timer clocki;    
   clocki.restart();
   ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
   map<int, string> ID_pn;
   get_pn(cf_fh.get_cf_val("file_pn"), ID_pn);      


   vector<string> chrns;  
   for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
   chrns.push_back("chrX");
   chrns.push_back("chrY");  

   if (opt == "1") { 
     int idx_pn;
     seqan::lexicalCast2(idx_pn, argv[3]);
     string pn = ID_pn[idx_pn];      
     string posfile_prefix = cf_fh.get_cf_val("file_clip_reads");
     string posfile_suffix = cf_fh.get_cf_val("freq_min") + "_" + cf_fh.get_cf_val("freq_max");

     string bam_input = cf_fh.get_cf_val("file_bam_prefix") + pn + ".bam";

#ifdef DEBUG_MODE   
     cerr << "!!!! DEBUG_MODE, using test bam input and chr1 \n";
     bam_input = cf_fh.get_cf_val("file_test_bam_prefix") + pn + ".bam";
#endif
     string bai_input = bam_input + ".bai";  
     BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bai_input); 

#ifdef DEBUG_MODE   // only look at chr1 for now 
     cerr << "!!!! please recompile if you want to run on all chromosomes\n";
     chrns.clear();
     chrns.push_back("chr1");  
 #endif       

     // reading groups 
     string file_rg = cf_fh.get_cf_val("file_dist_prefix") + "RG." + pn;
     map<string, int> rg_to_idx;
     parse_reading_group( file_rg, rg_to_idx );
     AluconsHandler *alucons_fh = new AluconsHandler(cf_fh.get_cf_val("file_alu_cons"), cf_fh.get_cf_val("type_alu_cons"));
     string file_fa = cf_fh.get_cf_val("file_fa_prefix");         
     string fout_path0 = cf_fh.get_cf_val("file_ins_del0");
     check_folder_exists( fout_path0);
     
     for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
       string posfile = posfile_prefix + *ci + "_" + posfile_suffix;
       vector< pair<int, int> > insert_pos;
       read_first2col(posfile, insert_pos);
       list <RecordInfo *> aluinfo; 
       search_chrn(*ci, insert_pos, rg_to_idx, aluinfo, bam_fh);
       cout << "size 0 " << aluinfo.size() << endl;
       pair_is_alu(aluinfo, alucons_fh, file_fa);
       map <int, int> insertPos_cnt;
       keep_pos_of_enough_reads(aluinfo, LEFT_PLUS_RIGHT, insertPos_cnt);
       cout << "size 2 " << aluinfo.size() << endl;             
       filter_clip_read(aluinfo, alucons_fh, file_fa); // can be improved       
       write_cnt0(aluinfo, fout_path0 + pn + ".tmp1", LEFT_PLUS_RIGHT);
       RecordInfo::delete_record_list(aluinfo);
     }
     delete bam_fh;    
     delete alucons_fh;
   }  else if ( opt == "debug") {
     string ref_fa = "CAGCACTTTGGGAGGCCGAGGCGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAATACGGTGAAACCCCGTCTCTACTAAAAATACAAAAA";
     //seqan::String<seqan::Dna> ref_fa = "CAGCACTTTGGGAGGCCGAGGCGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAATACGGTGAAACCCCGTCTCTACTAAAAATACAAAAA";
     string alucons = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAA";     
     
   } else if ( opt == "write_testbam") { // create tmp bam file, for testing 
     int idx_pn;
     seqan::lexicalCast2(idx_pn, argv[3]);
     string pn = ID_pn[idx_pn];      
     string posfile_prefix = cf_fh.get_cf_val("file_clip_reads");
     string posfile_suffix = cf_fh.get_cf_val("freq_min") + "_" + cf_fh.get_cf_val("freq_max");

     string bam_input = cf_fh.get_cf_val("file_bam_prefix") + pn + ".bam";
     string bai_input = bam_input + ".bai";  
     string bam_output = cf_fh.get_cf_val("file_test_bam_prefix") + pn + ".bam";
     BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bai_input, bam_output); 
     
     vector< pair<int, int> > insert_pos;
     seqan::BamAlignmentRecord record;
     int region_begin, region_end;
     string read_status;
     chrns.clear();
     chrns.push_back("chr1");
     chrns.push_back("chr2");
     
     for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
       // posfile manually created, avoid one read appear multiple times if pos are close to each other
       string posfile = posfile_prefix + *ci + "_" + posfile_suffix + ".testbam";
       read_first2col(posfile, insert_pos);

       for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
	 region_begin = (*pi).first - ALU_FLANK;
	 region_end = (*pi).second + ALU_FLANK;	
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

