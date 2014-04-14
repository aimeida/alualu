// utils function for alu_insert 
#include "insert_utils.h"

//write down alu_mate and flag, 
bool alu_mate_flag_slow( string bam_input, map<int, seqan::CharString> const &rID_chrn, MapFO &fileMap){
  /* it considers ilong_RC, dif_chr and mul_map. too many cases and therefore might be slow
     eg: count for chr1, 1320803 dif_chr, 127630 ilong_RC, 18 mul_map
     to improve: ignore mul_map !!
  */
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

bool readbam_clip(seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, size_t ref_begin, int ref_end, seqan::BamStream &bamStreamOut, vector<seqan::CharString> & qnames, size_t offset){
  bool hasAlignments = false;
  if (!jumpToRegion(inStream, hasAlignments, context, rID, ref_begin, ref_end, baiIndex)) return false;
  if (!hasAlignments) return false;
  seqan::BamAlignmentRecord record;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( record.rID != rID || record.beginPos >= ref_end) break;
    if ( record.beginPos + getAlignmentLengthInRef(record) <= ref_begin) continue; // has overlap with this region           
    if ( hasFlagQCNoPass(record) or hasFlagDuplicate(record) or hasFlagUnmapped(record) or hasFlagNextUnmapped(record) or (not hasFlagMultiple(record)) ) continue;
    if ( has_soft_last(record, offset) or has_soft_first(record, offset) ) {
      writeRecord(bamStreamOut, record);
      qnames.push_back(record.qName);
    }
  }
  return true;
}

void write_from_bam(string path_fout, string fin_read_bam, string fin_read_map, map <string, set< pair<int, int> > > & chr_clipReads, string format) {
  int pre_pos = 0;
  int region_begin, region_end;
  string chrn, line; 
  stringstream ss;
  ofstream fout;

  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, fin_read_bam.c_str(), "r"));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );
  ifstream fin(fin_read_map.c_str());
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    getline(fin, line);
    ss.clear(); ss.str(line);
    ss >> chrn >> region_begin >> region_end; // ignore rest of the row 
    if ( (!pre_pos) or (region_begin != pre_pos ) ) { 
      if (pre_pos)
	fout.close();	
      string fout_reads = path_fout + chrn + "_"+ int_to_string(region_begin) + "_"+ int_to_string(region_end) + "." + format;
      if ( chr_clipReads[chrn].find( make_pair(region_begin, region_end) ) == chr_clipReads[chrn].end()) {
	fout.open( fout_reads.c_str() );
	chr_clipReads[chrn].insert( make_pair(region_begin, region_end) );
      } else {
	fout.open( fout_reads.c_str(), ios::app); 	
      }
    }
    pre_pos = region_begin;
    if (format == "fastq")
      writeRecord(fout, record.qName, record.seq, record.qual, seqan::Fastq());
    else if  (format == "fa")
      writeRecord(fout, record.qName, record.seq, seqan::Fasta());
    else 
      cerr <<  "unknown option. Only support fastq or fa format\n";
  }
  fin.close();
  fout.close();
  seqan::close(inStream);
}


// called by seqcons. eg:  611,634[id=ReadId,fragId=PairId,repeatId=0]
// also need: simulated_reads.fastaF file. eg: >0[libId=1] mean, sd(libSize)
// output record.tLen can be wrong !! set simulate data example 
void write_fastq_seqcons(string path_fout, string fin_read_bam, string fin_read_map, map <string, set< pair<int, int> > > & chr_clipReads, string format) {
  // NB THINK, what the input should be !
  /*
   >15,50[id=0,fragId=0,repeatId=0]
TCACCAGCGCATAACACCTCTCACAACTCCCATTG
    >1,39[id=1,fragId=1,repeatId=0]
AAAGTCTGTGGACATCACCAGCGCATAACACCTCTCAC
    >3,38[id=2,fragId=2,repeatId=0]
AGTCTGTGGACATCACCAGCGCATAACACCTCTCA
    >8,43[id=3,fragId=3,repeatId=0]
  */  
}

void write_fastq_seqcons_nopair(string path_fout, string fin_read_bam, string fin_read_map, map <string, set< pair<int, int> > > & chr_clipReads, string format) {
  // NB THINK, what the input should be !
  /*
   >15,50[id=0,fragId=0,repeatId=0]
TCACCAGCGCATAACACCTCTCACAACTCCCATTG
    >1,39[id=1,fragId=1,repeatId=0]
AAAGTCTGTGGACATCACCAGCGCATAACACCTCTCAC
    >3,38[id=2,fragId=2,repeatId=0]
AGTCTGTGGACATCACCAGCGCATAACACCTCTCA
    >8,43[id=3,fragId=3,repeatId=0]
  */  
}
