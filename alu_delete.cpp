#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"

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

  EmpiricalPdf *empiricalpdf = new EmpiricalPdf(read_config(config_file, "file_dist"));  
  AluRefPos *alurefpos = new AluRefPos(file_alupos, alu_flank);
  ofstream bamStreamOut(bam_output.c_str());
  for (int count_loci = 0; count_loci < 2; count_loci++) {
    int beginPos = 0, endPos = 0;      
    if (!alurefpos->updatePos(beginPos, endPos)) continue;
    bool hasAlignments = false;
    if (!jumpToRegion(inStream, hasAlignments, context, rID, beginPos, endPos, baiIndex)) {
      std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
      return 1;
    }
    if (!hasAlignments) continue;    
    int reads_num = 0;
    vector<int> tmp_vec;
    ReadsPosStore reads_pos; 
    ReadsPosStore::iterator ri;    
    //cerr << "read pos " << beginPos << " " << endPos << endl;    
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
	//write2(bamStreamOut, record, context, seqan::Sam());
	reads_num ++;
	string record_qName = toCString(record.qName);
	if ( (ri = reads_pos.find(record_qName)) == reads_pos.end() ) { 
	  tmp_vec.push_back(record.beginPos);
	  tmp_vec.push_back(record.beginPos + length(record.seq));
	  reads_pos[record_qName].swap(tmp_vec);	  
	} else {
	  (ri->second).push_back(record.beginPos);
	  (ri->second).push_back(record.beginPos + length(record.seq));	
	}
      }
    }    
    float mean_coverage = 120.* reads_num / (endPos - beginPos);
    if (mean_coverage > coverage_max) continue; // skip high coverage region    
    cerr << "mean coverage: " << mean_coverage << "/" << coverage_max << endl;
    vector<int> reads_insert_len;
    float log_p[] = { 0, 0, 0 };
    genotype_prob(reads_insert_len, reads_pos, endPos - beginPos, log_p);
  }
  
  delete alurefpos;
  bamStreamOut.close();
  seqan::close(inStream);



  //cout << empiricalpdf->pdf_obs(100) << endl;
  //cout << empiricalpdf->pdf_obs(50) << endl;
  delete empiricalpdf;
  cout << "file_alupos:done  " << file_alupos << endl;
  
  return 0;

}
