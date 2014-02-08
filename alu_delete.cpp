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
  cout << int_to_string(chrn) << endl;
  string file_fa = read_config<string>(config_file, "file_fa_prefix") + int_to_string(chrn) + ".fa";
  string file_alupos = read_config<string>(config_file, "file_alupos_prefix") + int_to_string(chrn);

  cout << "file_alupos " << file_alupos << endl;
  cout << bam_input << endl;
  cout << bam_output << endl;

  // Open input stream, BamStream can read SAM and BAM files.  
  seqan::BamStream bamStreamIn(bam_input.c_str());
  seqan::BamStream bamStreamOut(bam_output.c_str(), seqan::BamStream::WRITE);
  bamStreamOut.header = bamStreamIn.header;
  seqan::BamAlignmentRecord record;

  int count_i = 0;
  // while (!atEnd(bamStreamIn)) {
  while (count_i++ < 1000) {
    readRecord(record, bamStreamIn);
    if ( (not hasFlagDuplicate(record)) and (not hasFlagUnmapped(record)) and (not hasFlagQCNoPass(record)) ) {
      writeRecord(bamStreamOut, record);
      //cout << "## " << count_i << " " << record.seq << endl;
    }
  }
  seqan::flush(bamStreamOut);
  seqan::close(bamStreamOut);
  //exit(0);  

  // what happened in the following ???

//  ProbByQuantile *prob_quantile = new ProbByQuantile(read_config<string>(config_file, "file.dist"));
//  cout << "prob " << prob_quantile->prob_inputlen(100) << endl;
//  delete prob_quantile;

  return 0;

}
