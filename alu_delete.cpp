#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include "configfile.h"
#include "diststat.h"

int main( int argc, char* argv[] )
{
  string config_file = argv[1];
  // Open input stream, BamStream can read SAM and BAM files.
  seqan::BamStream bamStreamIn(argv[2]);
  seqan::BamStream bamStreamOut(argv[3], seqan::BamStream::WRITE);
  cout << "chrn " << read_config<int>(config_file, "chrn") << endl ; 

  ProbByQuantile *prob_quantile = new ProbByQuantile(read_config<string>(config_file, "file.dist"));
  cout << "prob " << prob_quantile->prob_inputlen(100) << endl;
  delete prob_quantile;

  return 0;

}
