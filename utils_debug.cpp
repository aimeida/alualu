// 1. check bam flag, 2. print details of a location 
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "utils.h"

void read_pn_chr(string &bamInput_file, string &chrn) {
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
  if(!getIdByName(nameStore, chrn, rID, nameStoreCache)) 
    cout << "ERROR: Reference sequence named "<< chrn << " not known.\n";
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

bool check_region(seqan::Stream<seqan::Bgzf> &inStream, seqan::BamStream &bamStreamOut, seqan::BamIndex<seqan::Bai> &baiIndex,TBamIOContext &context, int rID, int beginPos, int endPos ){
  bool hasAlignments = false;
  if (!jumpToRegion(inStream, hasAlignments, context, rID, beginPos, endPos, baiIndex)) return 0;
  if (!hasAlignments) return 0;
  seqan::BamAlignmentRecord record;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != rID || record.beginPos >= endPos) break;
    if (record.beginPos < beginPos) continue;            
    //if ( hasFlagQCNoPass(record) or hasFlagDuplicate(record) or hasFlagUnmapped(record) or hasFlagNextUnmapped(record) or (not hasFlagMultiple(record))) continue;
    if ( !(QC_insert_read) )  continue;
    ///if ( ! hasFlagAllProper(record) )  // print out all records for now
    writeRecord(bamStreamOut, record);
  }
}

int main( int argc, char* argv[] )
{
  int idx_pn;  // start from 0
  seqan::lexicalCast2(idx_pn, argv[1]);
  string chrn = argv[2];  // chr1

  string pn_file = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/PN_all";
  string pn = get_pn(pn_file, idx_pn);
  cerr << "pn " << pn << endl;
  string bam_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/bamFile/" + pn + ".bam";
  //read_pn_chr(bamInput_file, chrn);
  string bai_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/bamFile/" + pn + ".bam.bai"; 
  seqan::BamAlignmentRecord record;  
  string fa_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/chromFa/" + chrn + ".fa";

  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  seqan::BamIndex<seqan::Bai> baiIndex;
  assert (!read(baiIndex, bai_input.c_str()));
  
  assert(!readRecord(header, context, inStream, seqan::Bam()) );
  seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
  bamStreamOut.header = header;

  int rID = 0;
  assert (getIdByName(nameStore, chrn, rID, nameStoreCache));

  string qname;
  int p, p_ref_a, p_ref_b, s, p2;
  int ref_fa = 120;
  

  // if read gap > 100 bp, might be the place of insertion
  // possible locations: 100766800  or 100766810 ! 
  check_region(inStream, bamStreamOut, baiIndex, context, rID, 183203676 - 600,  183203987 + 600);
  exit(0);

  seqan::CharString fa_seq;
  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////  

//   qname = "8:93:11692:3291";
//   p = 100766799;
//   if (find_read(bam_input, bai_input, chrn, qname, p, record, 0))
//     writeRecord(bamStreamOut, record);
//   
//   qname = "8:64:15011:5948";
//   if (find_read(bam_input, bai_input, chrn, qname, p, record, 0))
//     writeRecord(bamStreamOut, record);
  
  fa_seq = fasta_seq(fa_input, "chr1", 35351533, 35351633, true);
  cout << "fa_seq: \n";
  cout << fa_seq << endl;  

  fa_seq = fasta_seq(fa_input, "chr1", 211386483, 211386583, true);
  cout << "fa_seq: \n";
  cout << fa_seq << endl;  
  
  exit(0);

  cout << "## " << not_all_match(record) << endl;
  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
  TAlign align;
  int score;
  string suba;  
  resize(rows(align), 2);
  qname = "5:116:6455:19955";
  p = 1079720;
  p_ref_a = 1079720;
  p_ref_b = 1079749;
  find_read(bam_input, bai_input, chrn, qname, p, record, 0);
  writeRecord(bamStreamOut, record);
  assignSource(row(align,0),suba.substr(0, 27));
  assignSource(row(align,1),fa_seq);
  seqan::Score<int> scoringScheme(1, -3, -2, -5);
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, false, false>());  
  cout << "Score: " << score << endl;
  cout << align << endl;
  cout << "##########################\n";
  p_ref_a = 1079945;
  p_ref_b = 1080245;
  fa_seq = fasta_seq(fa_input, "chr1", p_ref_a, p_ref_b, true);
  cout << suba.substr(27) << endl;
  cout << fa_seq << endl;  

  return 0;
}

