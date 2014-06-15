#include "utils.h"

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
  //cerr << "read pdf dist " << pdf_file << endl;
  cerr << "min_len max_len bin_width min_prob\n";
  cerr << min_len << " " << max_len << " " << bin_width << " " << min_prob << endl;
}

float EmpiricalPdf::pdf_obs(int insertlen) {
  if (insertlen >= max_len || insertlen <= min_len) return min_prob;
  int nearby_pos = min_len + (insertlen-min_len)/bin_width*bin_width;
  map<int, float >::iterator it = prob_vec.find(nearby_pos);
  if ( it != prob_vec.end())  return it->second;
  return min_prob;
}   

void EmpiricalPdf::delete_map(map <int, EmpiricalPdf *> & epdf_rg){
  for (map <int, EmpiricalPdf *>::iterator ri = epdf_rg.begin(); ri != epdf_rg.end(); ri++) 
    delete ri->second;
  epdf_rg.clear();
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

bool BamFileHandler::jump_to_region(int rid, int region_begin, int region_end) {
  assert ( !fn_bai.empty()); // need bai file if want to jump 
  bool hasAlignments = false;	
  if (!jumpToRegion(inStream, hasAlignments, context, rid, region_begin, region_end, baiIndex))
    return false;
  return hasAlignments;
}

bool BamFileHandler::jump_to_region(string chrn, int region_begin, int region_end) {
  assert ( !fn_bai.empty()); // need bai file if want to jump 
  bool hasAlignments = false;	
  if (!jumpToRegion(inStream, hasAlignments, context, chrn_rID[chrn], region_begin, region_end, baiIndex))
    return false;
  return hasAlignments;
}

bool BamFileHandler::fetch_a_read(seqan::BamAlignmentRecord & record) {
  if (atEnd(inStream)) return false;
  assert (!readRecord(record, context, inStream, seqan::Bam())); 
  return true;
}

string BamFileHandler::fetch_a_read(string chrn, int region_begin, int region_end, seqan::BamAlignmentRecord & record) {
  if (atEnd(inStream)) return "stop";
  assert (!readRecord(record, context, inStream, seqan::Bam())); 
  if (record.rID != chrn_rID[chrn] || record.beginPos >= region_end) return "stop";
  if (record.beginPos + (int)getAlignmentLengthInRef(record) < region_begin) return "skip";
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
  if ( read(faiIndex, fn_fa.c_str()) ) {
    build(faiIndex, fn_fa.c_str() ); 
    write(faiIndex, (fn_fa+".fai").c_str() ); 
  }      
  only_one_seq = false;
}

FastaFileHandler::FastaFileHandler(string fn_fa, string seq_name) {
  if ( read(faiIndex, fn_fa.c_str()) ) {
    build(faiIndex, fn_fa.c_str() ); 
    write(faiIndex, (fn_fa+".fai").c_str() ); 
  }      
  idx = 0;
  assert (getIdByName(faiIndex, seq_name, idx));
  only_one_seq = true;
}

void FastaFileHandler::fetch_fasta_upper(int beginPos, int endPos, seqan::CharString &seq){
  assert (only_one_seq);
  assert (!readRegion(seq, faiIndex, idx, beginPos, endPos));
  seqan::toUpper(seq);
}

void FastaFileHandler::fetch_fasta_upper(int beginPos, int endPos, seqan::CharString &seq, string seq_name){
  assert (!only_one_seq);
  assert (getIdByName(faiIndex, seq_name, idx));
  assert (!readRegion(seq, faiIndex, idx, beginPos, endPos));
  seqan::toUpper(seq);
}

seqan::CharString FastaFileHandler::fasta_seq(string fa_input, string seq_name,int beginPos, int endPos) {
  FastaFileHandler * fasta_fh = new FastaFileHandler(fa_input, seq_name);
  seqan::CharString fa_seq;
  fasta_fh -> fetch_fasta_upper(beginPos, endPos, fa_seq);
  delete fasta_fh;
  return fa_seq;
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


void parse_reading_group(string file_rg, map<string, int> & rg_to_idx){
  ifstream fin (file_rg.c_str() );
  int idx = 0;
  assert(fin);
  string rg;
  while (fin >> rg) rg_to_idx[rg] = idx++;
  fin.close();
  //cout << rg_to_idx.size() << " reading groups exist\n";
}

int get_rgIdx(map<string, int> & rg_to_idx, seqan::BamAlignmentRecord & record) {
  seqan::BamTagsDict tags(record.tags);
  unsigned idx_rg1;  // idx in bam file
  assert (findTagKey(idx_rg1, tags, "RG"));
  string rg_name = seqan::toCString(getTagValue(tags, idx_rg1));
  int idx_rg2 = 0; // idx in my file 
  map <string, int>::iterator rgi = rg_to_idx.find(rg_name);
  if (rgi != rg_to_idx.end()) idx_rg2 = rgi->second;
  return idx_rg2;
}

void parse_cigar(string cigar, list <char> & opts, list <int> & cnts){
  string cnt;
  int cnt_int;
  opts.clear();
  cnts.clear();
  for (size_t i = 0; i < cigar.size(); i++) {
    if ( !isdigit(cigar[i]) ) {
      opts.push_back(cigar[i]);
      if (!cnt.empty()) {
	seqan::lexicalCast2(cnt_int, cnt);
	cnts.push_back(cnt_int);
      }
      cnt = "";
    } else {
      cnt += cigar[i];
    }
  }
  seqan::lexicalCast2(cnt_int, cnt);
  cnts.push_back(cnt_int);  
}

string get_cigar(seqan::BamAlignmentRecord &record) {
  stringstream ss;
  for (unsigned li = 0; li < length(record.cigar); li++) 
    ss << record.cigar[li].operation << record.cigar[li].count;
  return ss.str();
}

void debug_print_read(seqan::BamAlignmentRecord &record, ostream & os) {
  os << record.qName << " " << record.beginPos << " " << record.beginPos + getAlignmentLengthInRef(record)  << " " ;
  os << get_cigar(record) << " " << record.pNext << " " << record.tLen << endl;
}

int numOfBestHits(seqan::BamAlignmentRecord &record){
  seqan::BamTagsDict tagsDict(record.tags);
  unsigned myIdx = 0, valInt = 1;  
  if ( findTagKey(myIdx, tagsDict, "X0")) extractTagValue(valInt, tagsDict, myIdx);
  return max((int)valInt, 1);
}

bool find_read(string &bam_input, string &bai_input, string &chrn, string &this_qName, int this_pos, seqan::BamAlignmentRecord &that_record, int flank_region) { // flank_region = 0, print this read; otherwise, print its pair
  int that_begin = this_pos - max(flank_region, 10); 
  int that_end = this_pos + max(flank_region, 10);
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;  
  seqan::Stream<seqan::Bgzf> inStream;  
  open(inStream, bam_input.c_str(), "r");
  seqan::BamIndex<seqan::Bai> baiIndex;
  assert( !read(baiIndex, bai_input.c_str())) ;
  int rID;
  assert ( !readRecord(header, context, inStream, seqan::Bam()) ); 
  assert ( getIdByName(nameStore, chrn, rID, nameStoreCache) ); // change context ??
  bool hasAlignments = false;
  jumpToRegion(inStream, hasAlignments, context, rID, that_begin, that_end, baiIndex);
  if (!hasAlignments) return false;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != rID || record.beginPos >= that_end) break;
    if (record.beginPos < that_begin) continue;
    if (record.qName == this_qName) {
      if ((flank_region > 0 and record.beginPos != this_pos) or (flank_region == 0) ) {
	that_record = record;
	return true;
      }
    }
  }
  return false;
}

RepMaskPos::RepMaskPos(string file_rmsk, int join_len){
  vector<string>::iterator ci;
  vector<int>::iterator pi;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  chrns.push_back("chrX");
  chrns.push_back("chrY");

  string line, chrn, chrn_pre="chr0";
  stringstream ss;
  int beginPos, endPos, beginPos_pre, endPos_pre;
  ifstream fin(file_rmsk.c_str());
  assert(fin);
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> chrn >> beginPos >> endPos;
    if ( find(chrns.begin(), chrns.end(), chrn) == chrns.end() ) continue;
    if (chrn != chrn_pre) {
      chrn_pre = chrn;
      beginPos_pre = beginPos;
      endPos_pre = endPos;      
    } else {
      if (beginPos - endPos_pre > join_len) {  // close this block, create new
	beginP[chrn].push_back(beginPos_pre);
	endP[chrn].push_back(endPos_pre);	
	beginPos_pre = beginPos;
      }
      endPos_pre = endPos;      
    }
  }
  fin.close();    
//  for ( vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++ ) {
//    assert(beginP[*ci].size() == endP[*ci].size() );
//    cerr << *ci << " " <<  beginP[*ci].size() << endl;
//  }
}

AluRefPosRead::AluRefPosRead(string file_alupos, int minLen) {
  ifstream fin( file_alupos.c_str());
  if (!fin) 
    try {
      throw 1;    
    } catch(int e) {
      cerr << "#### ERROR #### file: "<< file_alupos << " not exists!!!" << endl;
    }     
  string _tmp1, _tmp2, alu_type;
  char chain;
  int bp, ep;
  while (fin >> _tmp1 >> bp >> ep >> alu_type >> _tmp2 >> chain){
    if (ep - bp < minLen) continue;
    if (beginP.empty()) minP = bp;
    beginP.push(bp);
    endP.push(ep);
    strandP.push(chain);
    aluType.push(alu_type);
  }
  maxP = ep;
  cerr << "queue from " << file_alupos << " with " << beginP.size() << " loci, " << minP << " to " << maxP << endl;
  fin.close();
}

int AluRefPosRead::updatePos(int &beginPos, int &endPos){
  if (!beginP.size()) return 0;
  beginPos = beginP.front();
  endPos = endP.front();
  beginP.pop();
  endP.pop();
  return 1;
}  

int AluRefPosRead::updatePos(int &beginPos, int &endPos, char &chain, string & alu_type){
  if (!beginP.size()) return 0;
  beginPos = beginP.front();
  endPos = endP.front();
  chain = strandP.front();
  alu_type = aluType.front();
  beginP.pop();
  endP.pop();
  strandP.pop();
  aluType.pop();
  return 1;
}  

AluRefPos::AluRefPos(string file_alupos) {
  ifstream fin( file_alupos.c_str());
  assert(fin);
  string _tmp1, _tmp2, alu_type, _tmp3;
  int bp, ep;
  while (fin >> _tmp1 >> bp >> ep >> alu_type >> _tmp2 >> _tmp3){
    beginV.push_back(bp);
    endV.push_back(ep);
    typeV.push_back(alu_type);
  }
  fin.close();  
}



bool AluRefPos::insideAlu(int beginPos, int endPos, int alu_min_overlap, int &len_overlap){
  vector<int>::iterator bi, ei;
  for (bi = beginV.begin(), ei = endV.begin(); bi != beginV.end(); bi++, ei++){
    if ( *bi > endPos ) return false;
    if ( (len_overlap = min(*ei,endPos)-max(*bi, beginPos)) >= alu_min_overlap ) return true;
  }
  return false;
}

AluRefPos::~AluRefPos(void){
  beginV.clear();
  endV.clear();
  typeV.clear();
}

string phred_scaled(float p0, float p1, float p2) {
  float pmax = max(p0, max(p1, p2));
  string s0 =  ( p0 == pmax ) ?  "0" : phred_log(p0/pmax);
  string s1 =  ( p1 == pmax ) ?  "0" : phred_log(p1/pmax);
  string s2 =  ( p2 == pmax ) ?  "0" : phred_log(p2/pmax);
  return s0  + "," + s1 + "," + s2; 
}

void read_pdf_pn( string prefix, string pn, string pdf_param,  map <int, EmpiricalPdf *> & pdf_rg) {
  ifstream fin( get_name_rg(prefix, pn).c_str() );
  int idx = 0; 
  string rg;
  while (fin >> rg) 
    pdf_rg[idx++] = new EmpiricalPdf( get_name_rg_pdf(prefix, pn, rg, pdf_param));
  fin.close();
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
  for ( i = 0; i < 3; i++) {
    idx_logp[i] = log_gp[i];
    assert (log_gp[i] <= 0);
  }
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
    else gp[i] = pow(10, (idx_logp[i] - min_logp)) / ratio_sum;
  }    
}

void read_first2col(string fn, vector < pair<int, int> > & insert_pos, bool has_header, int maxi ) {
  insert_pos.clear();
  ifstream fin(fn.c_str());
  assert(fin);
  stringstream ss;
  string line;
  if (has_header)  getline(fin, line); 
  int beginPos, endPos;
  int ni = 0;
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );  
    ss >> beginPos >> endPos;
    ni ++ ;
    if ( maxi > 0 and ni == maxi)
      break;
    insert_pos.push_back( make_pair(beginPos, endPos) );
  }
  fin.close();
}
