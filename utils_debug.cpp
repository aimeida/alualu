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
    cout << "ERROR: Reference sequence named "<< chrx << " not known.\n";
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
	cout << record.beginPos << " " <<  length(record.seq) << "\n";
      if (record.beginPos < record.pNext) {
	left_counts++;
      } else {
	right_counts++;
	continue;	
      }
      size_t align_len = getAlignmentLengthInRef(record);
      if (align_len < 90) {   // majority is soft klipped 	
	writeRecord(bamStreamOut, record);
	// cout << record.beginPos << " " << align_len << " " << length(record.seq) << " ";
	// cigar_line(record.cigar, length(record.cigar) ); 
	// cout << endl;	
      }
    }    
  }  
  seqan::close(bamStreamIn);
  cout << "## " <<  left_counts << ",  " << right_counts << endl;
}

int main( int argc, char* argv[] )
{
  int idx_pn;  // start from 0
  seqan::lexicalCast2(idx_pn, argv[1]);
  string chrx = argv[2];  // chr1

  string pn_file = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/PN_all";
  string pn = get_pn(pn_file, idx_pn);
  cout << "pn " << pn << endl;
  string bam_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/bamFile/" + pn + ".bam";
  //read_pn_chr(bamInput_file, chrx);
  string bai_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/bamFile/" + pn + ".bam.bai"; 
  seqan::BamAlignmentRecord record;  

  string fa_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/chromFa/" + chrx + ".fa";
  seqan::CharString fa_seq;
  string qname;
  int p, p_ref_a, p_ref_b, s, p2;
  int ref_fa = 120;
  string suba;

  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  assert(!readRecord(header, context, inStream, seqan::Bam()) );

  seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
  bamStreamOut.header = header;
  
  
  TAlign align;
  int score;
  resize(rows(align), 2);

  qname = "5:116:6455:19955";
  p = 1079720;
  p_ref_a = 1079720;
  p_ref_b = 1079749;
  find_read(bam_input, bai_input, chrx, qname, p, record, 0);
  writeRecord(bamStreamOut, record);
  cout << "\n##########################\n";
  cout << "hasRC " << hasFlagRC(record) << endl; // 1
  suba = toCString(record.seq);
  fa_seq = fasta_seq(fa_input, "chr1", p_ref_a, p_ref_b, true);
  cerr << fa_seq << endl;
  assignSource(row(align,0),suba.substr(0, 27));
  assignSource(row(align,1),fa_seq);
  // score: match, mismatch, gap[, gapOpen]
  seqan::Score<int> scoringScheme(1, -3, -2, -5);
  //score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>());
  // force to align in the beginning of the sequences 
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, false, false>());  
  cout << "Score: " << score << endl;
  cout << align << endl;
  cout << "##########################\n";

  // cross alignment 
  
  p_ref_a = 1079945;
  p_ref_b = 1080245;
  fa_seq = fasta_seq(fa_input, "chr1", p_ref_a, p_ref_b, true);
  cout << suba.substr(27) << endl;
  cout << fa_seq << endl;

  /*
  assignSource(row(align,0),suba.substr(27, 100));
  assignSource(row(align,1),fa_seq);
  //score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>());
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, false, false>());  
  cout << "Score: " << score << endl;
  cout << align << endl;
  cout << "##########################\n";
  */
  
  qname = "5:101:15310:5772";
  p = 1688516;
  find_read(bam_input, bai_input, chrx, qname, p, record, 0);
  writeRecord(bamStreamOut, record);
  cout << "\n##########################\n";
  cout << "hasRC " << hasFlagRC(record) << endl;  // 0
  suba = toCString(record.seq);
  cout << suba.substr(0,40) << endl;
  fa_seq = fasta_seq(fa_input, "chr1", 1688218 - 40, 1688218);
  cout << fa_seq << endl;
  cout << "##########################\n";

  exit(0);
  // check this cross mapping: chr1 1913466 1913765 1913764 1913810 S74M46 5:56:18835:12168
  qname = "5:56:18835:12168";
  p = 1913764;
  s = 74;
  p2 = 1913466;
  find_read(bam_input, bai_input, chrx, qname, p, record, 0);
  print_cigar(record);
  cout << " " << record.beginPos << " " << hasFlagRC(record) << " " << getAlignmentLengthInRef(record) << endl;
  suba = toCString(record.seq);
  cout << suba.substr(s) << endl;
  cout << suba.substr(0, s) << endl;

  fa_seq = fasta_seq(fa_input, "chr1", p, p+ref_fa);
  cout << fa_seq << endl;
  fa_seq = fasta_seq(fa_input, "chr1", p2 - 200 , p2 + 200); // call local alignment externally
  cout << fa_seq << endl;
  
  //chr1 1688218 1688515 1688516 1688576 S41M60 5:101:15310:5772
  qname = "5:101:15310:5772";
  p = 1688516;
  s = 41;
  p2 = 1688218;
  find_read(bam_input, bai_input, chrx, qname, p, record, 0);
  print_cigar(record);
  cout << " " << record.beginPos << " " << hasFlagRC(record) << " " << getAlignmentLengthInRef(record) << endl;
  suba = toCString(record.seq);
  cout << suba.substr(s) << endl;
  cout << suba.substr(0, s) << endl;

  // soft clip
  qname = "5:7:12281:15766";
  p = 227391;
  s = 14;
  find_read(bam_input, bai_input, chrx, qname, p, record, 0);
  print_cigar(record);
  cout << " " << record.beginPos << " " << hasFlagRC(record) << " " << getAlignmentLengthInRef(record) << endl;
  cout << record.seq << endl;
  suba = toCString(record.seq);
  cout << suba.substr(s) << endl;

  return 0;
}

