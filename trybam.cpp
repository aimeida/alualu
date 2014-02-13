#define SEQAN_HAS_ZLIB 1
#include <iostream>
#include <map>
#include <string>
#include <iostream>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

int main( int argc, char* argv[] )
{
  std::cout << "test compile" << std::endl; 
  // Open input stream, BamStream can read SAM and BAM files.
  seqan::BamStream bamStreamIn(argv[1]);
  // Open output stream, "-" means stdin on if reading, else stdout.
  seqan::BamStream bamStreamOut(argv[2], seqan::BamStream::WRITE);
  bamStreamOut.header = bamStreamIn.header;

  seqan::BamAlignmentRecord record;
  //  while (!atEnd(bamStreamIn)) {
  int i = 0;
  while ( i++ < 10 ){
    readRecord(record, bamStreamIn);
    if( (not hasFlagDuplicate( record ) ) and (hasFlagUnmapped( record ) or hasFlagNextUnmapped( record )) and (not hasFlagSecondary( record )) and (not hasFlagQCNoPass(record)) )
      continue;
  }
  
  return 0;  
}
