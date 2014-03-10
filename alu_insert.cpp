#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

#define ALU_MINLEN 10
#define MAPPED_ALU_LEN 100

//search alu_mate by flag, to bam file and coordiate file
bool alu_mate_flag( string & bam_input, string &bai_input, string &file_alupos_prefix, vector<string> &chrns, ofstream &fout1) {
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  seqan::BamIndex<seqan::Bai> baiIndex;
  assert (!read(baiIndex, bai_input.c_str()));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );
  
  // for debug etc
  seqan::BamStream bamStreamOut("tmp.sam", seqan::BamStream::WRITE);
  bamStreamOut.header = header; 

  map<int, seqan::CharString> rID_chrx;
  map<int, AluRefPos *> rID_alurefpos;
  int rID, i = 0;
  ////for (i = 0; i < length( header.sequenceInfos); i++) rID_chrx[i] = header.sequenceInfos[i].i1;  
  for ( vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    assert(getIdByName(nameStore, *ci, rID, nameStoreCache));
    rID_chrx[rID] = *ci;
    rID_alurefpos[rID] = new AluRefPos(file_alupos_prefix + *ci, true);    
  }
  
  map<int, AluRefPos *>::iterator ra;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( hasFlagQCNoPass(record) or hasFlagDuplicate(record) or hasFlagUnmapped(record) or hasFlagNextUnmapped(record)) continue;
    /////if (hasFlagSecondary(record)) ==> ALWAYS FALSE

    /*
    if (record.rID < record.rNextId) { // don't check twice
      if ( ( (ra = rID_alurefpos.find(record.rNextId)) != rID_alurefpos.end() ) and (ra->second)->insideAlu(record.pNext, record.pNext + MAPPED_ALU_LEN, ALU_MINLEN) ) {
	fout1 << "dif_chr " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rNextId << " " << record.pNext << endl;	
	continue;
      } else if ( ( (ra = rID_alurefpos.find(record.rID)) != rID_alurefpos.end() ) and (ra->second)->insideAlu(record.beginPos, record.beginPos + getAlignmentLengthInRef(record), ALU_MINLEN) ) {
	fout1 << "dif_chr " << record.qName << " " << record.rNextId << " " << record.pNext  << " " << record.rID << " " << record.beginPos << endl;
	continue;
      }
      // if not in alu, maybe due to non-unique mapping ???
    } 
    */
    
    // Next step: if qName appears twice, ignore it !
    if (hasFlagMultiple(record))
      writeRecord(bamStreamOut, record);

      //fout1 << "mul_map " << record.qName << " " << record.rNextId << " " << record.pNext  << " " << record.rID << " " << record.beginPos << " " << hasFlagSecondary(record) << endl;        

  }
  
  seqan::close(inStream);     
  seqan::close(bamStreamOut);
  for (ra = rID_alurefpos.begin(); ra != rID_alurefpos.end(); ra++)
    delete ra->second;
  rID_alurefpos.clear();
  return 0;
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);

  int opt, idx_pn;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];
  seqan::lexicalCast2(idx_pn, argv[3]);

//  unsigned coverage_max;
//  seqan::lexicalCast2(coverage_max, (read_config(config_file, "coverage_max")));

  string pn = get_pn(read_config(config_file, "file_pn"), idx_pn);
  cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
  string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
  string bai_input = bam_input + ".bai";  
  string file_alu_mate_prefix = read_config(config_file, "file_alu_mate_prefix");  
  
  if (opt == 1) {
    string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
    vector<string> chrns;
    for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
    chrns.push_back("chrX");
    chrns.push_back("chrY");

    ofstream fout1( (file_alu_mate_prefix + pn + ".alumate_flag").c_str() );
    fout1 << "flag qname proper_chr proper_pos bad_chr bad_pos \n";
    alu_mate_flag(bam_input, bai_input, file_alupos_prefix, chrns, fout1);
    fout1.close();
  }  
  return 0;  
}
