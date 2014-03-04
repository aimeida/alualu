#include "utils.h"

string get_pn(string pn_file, int idx_pn){
  ifstream fin( pn_file.c_str());
  string pn;
  int i = 0;
  while (fin >> pn) {
    if (i == idx_pn ) {
      fin.close();  
      return pn;
    }
    i++;
  }
  return "na";
}

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

void print_cigar(seqan::String<seqan::CigarElement<> > &cigar, int len_cigar) {
  for (int li = 0; li < len_cigar; li++) cerr << cigar[li].operation << cigar[li].count; 
}

void fasta_seq(string fa_input, string chrx, int beginPos, int endPos, seqan::CharString &seq){
  seqan::FaiIndex faiIndex;
  assert (!read(faiIndex, fa_input.c_str()) );
  unsigned idx = 0;
  assert (getIdByName(faiIndex, chrx, idx));
  assert (!readRegion(seq, faiIndex, idx, beginPos, endPos));
}

bool find_read(string &bam_input, string &bai_input, string &chrx, string &this_qName, int this_pos, seqan::BamAlignmentRecord &that_record, int flank_region) { // flank_region = 0, print this read; otherwise, print its pair
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
  assert ( getIdByName(nameStore, chrx, rID, nameStoreCache) );
  bool hasAlignments = false;
  jumpToRegion(inStream, hasAlignments, context, rID, that_begin, that_end, baiIndex);
  if (!hasAlignments) return false;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != rID || record.beginPos >= that_end) break;
    if (record.beginPos < that_begin) continue;
    if (record.qName == this_qName) {
      if ((flank_region > 0 and record.beginPos != this_pos) or (flank_region == 0) ) {
	cerr << flank_region << " " <<  record.beginPos << " " << this_pos << endl;
	that_record = record;
	return true;
      }
    }
  }
  return false;
}

///// 1-based to 0-based.
AluRefPos::AluRefPos(string file_alupos) {
  ifstream fin( file_alupos.c_str());
  if (!fin) 
    try {
      throw 1;    
    } catch(int e) {
      cerr << "#### ERROR #### file: "<< file_alupos << " not exists!!!" << endl;
    }     
  string _tmp1, _tmp2, alu_type, chain;
  int bp, ep;
  while (fin >> _tmp1 >> bp >> ep >> alu_type >> _tmp2 >> chain){
    beginP.push(bp);
    endP.push(ep);
  }
  fin.close();
  cerr << "reading from " << file_alupos << " with " << beginP.size() << " loci\n" ;
}

int AluRefPos::updatePos(int &beginPos, int &endPos){
  if (!beginP.size()) return 0;
  beginPos = beginP.front();
  endPos = endP.front();
  beginP.pop();
  endP.pop();
  return 1;
}  

AluRefPos::~AluRefPos(void){
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

