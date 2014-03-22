#include "utils.h"

string read_config(string config_file, string key){
  string _key, _value;
  ifstream fin( config_file.c_str());
  if (!fin) 
    try {
      throw 1;    
    } catch(int e) {
      cerr << "#### ERROR #### file: "<< config_file << " not exists!!!" << endl;
    }   
  while (fin >> _key >> _value) {
    if (_key == key) 
      return _value;
  }
  fin.close();  
  try {
    throw 1;    
  } catch(int e) {
    cerr << "#### ERROR #### key: " << key <<" doesn't exist in file:" << config_file << endl;
    exit(0);
  }     
}

void print_cigar(seqan::BamAlignmentRecord &record) {
  for (unsigned li = 0; li < length(record.cigar); li++) cerr << record.cigar[li].operation << record.cigar[li].count; 
}

void print_cigar(seqan::BamAlignmentRecord &record, ofstream &fout) {
  for (unsigned li = 0; li < length(record.cigar); li++) fout << record.cigar[li].operation << record.cigar[li].count; 
}

void print_read(seqan::BamAlignmentRecord &record) {
  cerr << record.qName << " " << record.beginPos << " " << record.beginPos + getAlignmentLengthInRef(record)  << " " << record.pNext << " " << record.beginPos + record.tLen << " ";
  print_cigar(record);
  cerr << endl;
}

void print_read(seqan::BamAlignmentRecord &record, ofstream &fout) {
  fout << record.qName << " " << record.beginPos << " " << record.beginPos + getAlignmentLengthInRef(record)  << " " << record.pNext << " " << record.beginPos + record.tLen << " ";
  print_cigar(record, fout);
  fout << endl;
}

void get_rID_chrn(string & bam_input, vector<string> &chrns, map<int, seqan::CharString> &rID_chrn){
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );
  int rID;
  for ( vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    assert(getIdByName(nameStore, *ci, rID, nameStoreCache));
    rID_chrn[rID] = *ci;
  }  
  seqan::close(inStream);
}

char mappingType(seqan::BamAlignmentRecord &record){
  seqan::BamTagsDict tagsDict(record.tags);
  unsigned myIdx = 0;
  char valChar = 0;  
  if ( findTagKey(myIdx, tagsDict, "X0") and extractTagValue(valChar, tagsDict, myIdx)) return valChar;
  return 0;
}

int numOfBestHits(seqan::BamAlignmentRecord &record){
  seqan::BamTagsDict tagsDict(record.tags);
  unsigned myIdx = 0, valInt = 1;  
  if ( findTagKey(myIdx, tagsDict, "X0")) extractTagValue(valInt, tagsDict, myIdx);
  return max((int)valInt, 1);
}

seqan::CharString fasta_seq(string fa_input, string chrn, int beginPos, int endPos, bool upper){
  seqan::FaiIndex faiIndex;
  assert (!read(faiIndex, fa_input.c_str()) );
  unsigned idx = 0;
  assert (getIdByName(faiIndex, chrn, idx));
  seqan::CharString seq;
  assert (!readRegion(seq, faiIndex, idx, beginPos, endPos));
  if (upper) seqan::toUpper(seq);
  return seq;
//    seqan::ModifiedString< seqan::CharString, seqan::ModView<MyUpper> > SEQ(seq);
//    return SEQ;
}

seqan::CharString fasta_seq(seqan::FaiIndex &faiIndex, unsigned idx, int beginPos, int endPos, bool upper){
  seqan::CharString seq;
  assert (!readRegion(seq, faiIndex, idx, beginPos, endPos));
  if (upper) seqan::toUpper(seq);
  return seq;
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
  assert ( !readRecord(header, context, inStream, seqan::Bam()) ); // do i need this ??
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
	//cerr << flank_region << " " <<  record.beginPos << " " << this_pos << endl;
	that_record = record;
	return true;
      }
    }
  }
  return false;
}

bool find_read(seqan::Stream<seqan::Bgzf> &inStream, TBamIOContext &context, int rID, int this_pos, string this_qName, seqan::BamAlignmentRecord &record){
  //bool find_read(seqan::Stream<seqan::Bgzf> &inStream, int rID, int this_pos, string &this_qName ){
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != rID || record.beginPos >= this_pos + 5) break;
    if (record.beginPos < this_pos-5 ) continue;
    if (record.qName == this_qName and record.beginPos == this_pos) return true;
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

void RepMaskPos::print_begin(int ni){
  for ( vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++ ) {
    cerr << *ci << " " << beginP[*ci].size() << endl;
    vector<int>::iterator pi = beginP[*ci].begin(); 
    vector<int>::iterator ei = endP[*ci].begin(); 
    for (int i = 0; i < ni and pi != beginP[*ci].end(); i++) 
      cerr << *pi++ << " " << *ei++ << endl;
  }
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

int AluRefPosRead::updatePos(int &beginPos, int &endPos, char &chain){
  if (!beginP.size()) return 0;
  beginPos = beginP.front();
  endPos = endP.front();
  chain = strandP.front();
  beginP.pop();
  endP.pop();
  strandP.pop();
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
}

void insertLen_of_nonUniq_mapping(vector<int> &starts_ends, vector<int> &reads_insert_len) {
  // starts_ends = [start_0, end_0, start_1, end_1, start_2, end_2 ....]
  vector <int> pre_ends;
  vector <int>::iterator vi, pi;
  int cur_begin, cur_end, insert_len;
  vi = starts_ends.begin();
  pre_ends.push_back(*(++vi));
  vi++;
  while ( vi != starts_ends.end() ) {
    cur_begin = *vi++;
    cur_end = *vi++; 
    for (pi = pre_ends.begin(); pi != pre_ends.end(); pi++) {
      insert_len = cur_begin - *pi;
      if (insert_len > 0) reads_insert_len.push_back(insert_len);
    }
    pre_ends.push_back(cur_end);
  }
  pre_ends.clear();    
}

