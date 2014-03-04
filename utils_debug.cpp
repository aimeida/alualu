// 1. check bam flag, 2. print details of a location 
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "utils.h"

void read_pn_chr(string &bamInput_file, string &chrx) {
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;  
  seqan::Stream<seqan::Bgzf> inStream;
  open(inStream, bamInput_file.c_str(), "r");
  readRecord(header, context, inStream, seqan::Bam());
  seqan::close(inStream);
  int rID = 0;
  if(!getIdByName(nameStore, chrx, rID, nameStoreCache)) 
    cerr << "ERROR: Reference sequence named "<< chrx << " not known.\n";
  seqan::BamStream bamStreamIn(bamInput_file.c_str());
  seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
  bamStreamOut.header = bamStreamIn.header; // cannot skip !

  size_t i = 0;
  size_t left_counts = 0, right_counts = 0;
  while (!atEnd(bamStreamIn)) {
    readRecord(record, bamStreamIn);
    if (record.rID != rID) {
      if (!i) continue;
      else break;
    }
    if (i++ > 500000 ) break;    
    if ((not hasFlagQCNoPass(record) ) and hasFlagAllProper(record) and (not hasFlagDuplicate(record)) and hasFlagMultiple(record) ) {
      if ( length(record.cigar) < 1 ) 
	cerr << record.beginPos << " " <<  length(record.seq) << "\n";
      if (record.beginPos < record.pNext) {
	left_counts++;
      } else {
	right_counts++;
	continue;	
      }
      size_t align_len = getAlignmentLengthInRef(record);
      if (align_len < 90) {   // majority is soft klipped 	
	writeRecord(bamStreamOut, record);
	// cerr << record.beginPos << " " << align_len << " " << length(record.seq) << " ";
	// cigar_line(record.cigar, length(record.cigar) ); 
	// cerr << endl;	
      }
    }    
  }  
  seqan::close(bamStreamIn);
  cerr << "## " <<  left_counts << ",  " << right_counts << endl;
}

int main( int argc, char* argv[] )
{
  int idx_pn;  // start from 0
  seqan::lexicalCast2(idx_pn, argv[1]);
  string chrx = argv[2];  // chr1

  string pn_file = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/PN_noP";
  string pn = get_pn(pn_file, idx_pn);
  cerr << "pn " << pn << endl;
  string bam_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/bamFile/" + pn + ".bam";
  //read_pn_chr(bamInput_file, chrx);
  string bai_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/bamFile/" + pn + ".bam.bai"; 
  seqan::BamAlignmentRecord record;  

  string fa_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/chromFa/" + chrx + ".fa";
  seqan::CharString fa_seq;
  string qname = "5:89:13046:5131";
  int p = 531848 - 1; // 531848 if from sam file

  if (find_read(bam_input, bai_input, chrx, qname, p, record, 0)) {
    cerr << "find! " << record.qName << " " << record.beginPos << " " << record.pNext << endl;
    cerr << length(record.seq) << " " << hasFlagRC(record) << endl;
    cerr << record.seq << endl;

    fasta_seq(fa_input, "chr1", record.beginPos, record.beginPos+20, fa_seq);
    cerr << fa_seq << endl;
  }

  qname = "7:17:12645:11069";
  p = 705847 - 1;
  if (find_read(bam_input, bai_input, chrx, qname, p, record, 1000)) {
    cerr << "find! " << record.qName << " " << record.beginPos << " " << record.pNext << endl;
    cerr << length(record.seq) << " " << hasFlagRC(record) << endl;
    cerr << record.seq << endl;

    fasta_seq(fa_input, "chr1", record.beginPos, record.beginPos+20, fa_seq);
    cerr << fa_seq << endl;    
  }

//// find read_pair
//  if (find_read(bam_input, bai_input, chrx, qname, 531847, record, 1000)) 
//    cerr << "find! " << record.qName << " " << record.beginPos << " " << record.pNext << endl;
  

  return 0;
}

