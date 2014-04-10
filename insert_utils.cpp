// utils function for alu_insert 
#include "insert_utils.h"

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
    ss >> chrn >> region_begin >> region_end; // >> idx_pn_first;
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
