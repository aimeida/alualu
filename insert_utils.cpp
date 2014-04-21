// utils function for alu_insert 
#include "insert_utils.h"

bool alu_mate_flag( string bam_input, map<int, seqan::CharString> const &rID_chrn, MapFO &fileMap){
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );

  int rID_pre = -1, ia = 0;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( !QC_insert_read(record) ) continue;
    if ( rID_chrn.find(record.rID) == rID_chrn.end() ) continue;
    if ( record.rID != rID_pre) {
      if (rID_pre > -1) {
	cerr << "done with " << rID_pre << endl; 
      }
      rID_pre = record.rID;	
    }       
    if (ia++ % 1000000 == 0) cerr << ia << " " << endl;
    int len_seq = length(record.seq);
    if ( getAlignmentLengthInRef(record) < (size_t)(len_seq - 10)) continue; // want very good quality reads 
    if (record.rID < record.rNextId and rID_chrn.find(record.rNextId) != rID_chrn.end() ) { 
      *(fileMap[record.rNextId]) << "dif_chr " << record.qName << " " << record.rNextId << " " << record.pNext << " " 
				 << record.rID << " " << record.beginPos << " " << len_seq << endl;      
      *(fileMap[record.rID]) << "dif_chr " << record.qName << " " << record.rID << " " << record.beginPos << " " 
				 << record.rNextId << " " << record.pNext << " " << len_seq << endl;            
      continue;
    } 

    if (record.rID == record.rNextId and record.beginPos < record.pNext) { 
      if (abs(record.tLen) > DISCORDANT_LEN) {
	*(fileMap[record.rID]) << "ilong " << record.qName << " " << record.rID << " " << record.pNext << " " 
			       << record.rID << " " << record.beginPos << " " << len_seq << endl;      
	*(fileMap[record.rID]) << "ilong " << record.qName << " " << record.rID << " " << record.beginPos << " " 
			       << record.rID << " " << record.pNext << " " << len_seq << endl;      
      }    
    }
  }
  seqan::close(inStream);     
  return 0;
}

// only support 4 types for now 
string parse_alu_type(string alu_name){
  assert ( !alu_name.empty() );
  if ( alu_name.substr(0,4) == "AluY") return "AluY";
  if ( alu_name.substr(0,4) == "AluS") return "AluSx";
  if ( alu_name.substr(0,4) == "AluJ") {
    if ( alu_name.substr(0,5) == "AluJo") return "AluJo";
    if ( alu_name.substr(0,5) == "AluJb") return "AluJb";
    return "AluJo";
  }
  return "AluY";
}

void get_alu_type(string fn, map<int, string> &pos_aluType){
  stringstream ss;
  int region_begin, region_end;
  string line, tmp1, tmp2, alu_type;   
  pos_aluType.clear();
  ifstream fin(fn.c_str());
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> region_begin >> region_end >> tmp1 >> tmp2 >> alu_type;
    pos_aluType[region_begin] = alu_type;
  }
  fin.close();
}

/** it considers ilong_RC, dif_chr and mul_map. too many cases and therefore might be slow
    eg: count for chr1, 1320803 dif_chr, 127630 ilong_RC, 18 mul_map */
bool alu_mate_flag_depreciate( string bam_input, map<int, seqan::CharString> const &rID_chrn, MapFO &fileMap){
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );

  int rID_pre = -1, ia = 0;
  map<seqan::CharString, pair<int,int> > same_chr_left; // value: (pos, count)
  map<seqan::CharString, pair<int,int> >::iterator sc;

  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( !QC_insert_read(record) ) continue;
    if ( rID_chrn.find(record.rID) == rID_chrn.end() ) continue;
    if ( record.rID != rID_pre) {
      if (same_chr_left.size() ) {  // do not print out left reads 
	for (map<seqan::CharString, pair<int,int> >::iterator mc = same_chr_left.begin(); mc!= same_chr_left.end(); mc++)
	  cerr << "left " << mc->first << " " << rID_pre << " " << (mc->second).first << " -1 -1 " << (mc->second).second << " " << length(record.seq) << endl;
	same_chr_left.clear();
      }
      if (rID_pre > -1) cerr << "done with " << rID_pre << endl; //  execute 'rID_chrn[rID_pre]' will add key to the dict !!
      rID_pre = record.rID;	
    }       
    if (ia++ % 1000000 == 0) cerr << ia << " " << same_chr_left.size() << endl;
 
    if (record.rID != record.rNextId) { // NB: redundant! one pair read twice, thus are written twice as well
      if (rID_chrn.find(record.rNextId) != rID_chrn.end())  // both chrn exist
	*(fileMap[record.rNextId]) << "dif_chr " << record.qName << " " << record.rNextId << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << numOfBestHits(record) << " " << length(record.seq) << endl;      
      continue;
    }
    int hits = numOfBestHits(record);      
    if (record.beginPos < record.pNext) {	// left read
      if ( hits > 1 and not_all_match(record) ) { 
	same_chr_left[record.qName] = make_pair(record.beginPos, hits);
      } else if ( (abs(record.tLen) > DISCORDANT_LEN ) or (hasFlagNextRC(record) == hasFlagRC(record))) {  // here to define iLong_RC
	same_chr_left[record.qName] = make_pair(record.beginPos, FLAG_LONG); // non-redundant
      }	
    } else if (record.beginPos > record.pNext) { // right read
      sc = same_chr_left.find(record.qName);  
      
      if (sc == same_chr_left.end()) {   // left read is unique
	if (hits > 1 and not_all_match(record) )
	  *(fileMap[record.rID]) << "mul_map " << record.qName << " " << record.rID <<  " " << record.beginPos << " " << record.rID << " " << record.pNext << " " << hits << " " << length(record.seq) << endl;
      } else {                   // left read not multi_mapped
	if ( (sc->second).second == FLAG_LONG) {  // NB, redundant! written twice when right reads are read
	  *(fileMap[record.rID]) << "ilong_RC " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " 1" << hits << " " << length(record.seq) << endl;
	  *(fileMap[record.rID]) << "ilong_RC " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rID << " " << record.pNext << " " << hits << " " << length(record.seq) << endl;
	} else if (hits == 1) {  // at least one pair is ok
	  *(fileMap[record.rID]) << "mul_map " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << (sc->second).second << " " << length(record.seq)  << endl;
	}
	same_chr_left.erase(sc); // rm left read
      }
     
    }    
  } 
  seqan::close(inStream);     
  return 0;
}

