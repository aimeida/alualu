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

  int chrn = read_config<int>(config_file, "chrn");
  string file_fa = read_config<string>(config_file, "file_fa_prefix") + int_to_string(chrn) + ".fa";
  string file_alupos = read_config<string>(config_file, "file_alupos_prefix") + int_to_string(chrn);

  cout << int_to_string(chrn) << endl;
  cout << "file_alupos " << file_alupos << endl;

//  // Open BGZF Stream for reading.
//  seqan::Stream<seqan::Bgzf> inStream;
//  if (!open(inStream, bam_input.c_str(), "r")) {
//    std::cerr << "ERROR: Could not open " << bam_input << " for reading.\n";
//    return 1;
//  }
//
//  // Read BAI index.
//  seqan::BamIndex<seqan::Bai> baiIndex;
//  if (read(baiIndex, bai_input.c_str()) != 0){
//    cerr << "ERROR: Could not read BAI index file " << bai_input << endl;
//    return 1;
//  }
  
  // Open input stream, BamStream can read SAM and BAM files.  
  seqan::BamStream bamStreamIn(bam_input.c_str());
  seqan::BamStream bamStreamOut(bam_output.c_str(), seqan::BamStream::WRITE);
  bamStreamOut.header = bamStreamIn.header;
  seqan::BamAlignmentRecord record;

  int count_i = 0;
  while (!atEnd(bamStreamIn)) {
    if (count_i++ > 1000) break;
    readRecord(record, bamStreamIn);
    if ( (not hasFlagDuplicate(record)) and (not hasFlagUnmapped(record)) and (not hasFlagQCNoPass(record)) ) {
      writeRecord(bamStreamOut, record);
      //cout << "## " << count_i << " " << record.seq << endl;
    }
  }
  seqan::flush(bamStreamOut);
  seqan::close(bamStreamOut);

  cout << "file_alupos:done  " << file_alupos << endl;
  ProbByQuantile *prob_quantile = new ProbByQuantile(read_config<string>(config_file, "file_dist"));

  //////// fix me 
  ////////../libraries/seqan/basic/basic_exception.h:236 FAILED!  (Uncaught exception of type St11logic_error: basic_string::_S_construct NULL not valid)
  
  cout << "prob " << prob_quantile->prob_inputlen(100) << endl;
  delete prob_quantile;

  return 0;

}
