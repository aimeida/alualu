#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "read_files.h"

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  string config_file = argv[1];
  string bam_input = argv[2];
  string bam_output = argv[3];
  string bai_input = bam_input + ".bai";  

  string chrx = read_config(config_file, "chr");
  int coverage_max, alu_flank;
  seqan::lexicalCast2(coverage_max, (read_config(config_file, "coverage_max")));
  seqan::lexicalCast2(alu_flank, (read_config(config_file, "alu_flank")));
  
  string file_fa = read_config(config_file, "file_fa_prefix") + chrx + ".fa";
  string file_alupos = read_config(config_file, "file_alupos_prefix") + chrx;
  
  // Open BGZF Stream for reading.
  seqan::Stream<seqan::Bgzf> inStream;
  if (!open(inStream, bam_input.c_str(), "r")) {
    std::cerr << "ERROR: Could not open " << bam_input << " for reading.\n";
    return 1;
  }  
  // Read BAI index.
  seqan::BamIndex<seqan::Bai> baiIndex;
  if (read(baiIndex, bai_input.c_str()) != 0){
    cerr << "ERROR: Could not read BAI index file " << bai_input << endl;
    return 1;
  }
  // Setup name store, cache, and BAM I/O context.
  typedef seqan::StringSet<seqan::CharString> TNameStore;
  typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
  typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;


  if (readRecord(header, context, inStream, seqan::Bam()) != 0) {
    cerr << "ERROR: Could not read header from BAM file " << bam_input << "\n";
    return 1;
  }  
  // Translate from reference name to rID.
  int rID = 0;
  if (!getIdByName(nameStore, chrx, rID, nameStoreCache)) {
    cerr << "ERROR: Reference sequence named "<< chrx << " not known.\n";
    return 1;
  }
  
  AluRefPos *alurefpos = new AluRefPos(file_alupos, alu_flank);
  int beginPos = 0, endPos = 0;  

  ////seqan::BamStream bamStreamOut(bam_output.c_str(), seqan::BamStream::WRITE);
  ////writeRecord(bamStreamOut, record);
  ofstream bamStreamOut(bam_output.c_str());
  for (int count_loci = 0; count_loci < 2; count_loci++) {
    if (!alurefpos->updatePos(beginPos, endPos)) continue;
    cerr << "read pos " << beginPos << " " << endPos << endl;

    bool hasAlignments = false;
    if (!jumpToRegion(inStream, hasAlignments, context, rID, beginPos, endPos, baiIndex)) {
      std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
      return 1;
    }
    if (!hasAlignments) continue;    
    int reads_num = 0;  
    while (!atEnd(inStream)) {
      if (readRecord(record, context, inStream, seqan::Bam()) != 0) {
	std::cerr << "ERROR: Could not read record from BAM file.\n";
	exit(0);
      }     
      // If we are on the next reference or at the end already then we stop.
      if (record.rID != rID || record.beginPos >= endPos) break;
      // If we are left of the selected posi2tion then we skip this record.
      if (record.beginPos < beginPos) continue;      
      if ( (not hasFlagDuplicate(record)) and (not hasFlagUnmapped(record)) and (not hasFlagQCNoPass(record)) ) {	
	write2(bamStreamOut, record, context, seqan::Sam());
	//cout << record.tLen << " ==? "<< length(record.seq) << endl;
	reads_num++ ;
      }
    }    
    float mean_coverage = 120.* reads_num / (endPos - beginPos);
    if (mean_coverage > coverage_max) continue; // skip high coverage region
    cout << "reads num: " << reads_num << ", mean coverage: " << mean_coverage << "/" << coverage_max << endl;
  }
  
  delete alurefpos;
  bamStreamOut.close();
  seqan::close(inStream);

  cout << "file_alupos:done  " << file_alupos << endl;
  ProbByQuantile *prob_quantile = new ProbByQuantile(read_config(config_file, "file_dist"));
  cout << "prob " << prob_quantile->prob_inputlen(100) << endl;
  delete prob_quantile;

  return 0;

}
