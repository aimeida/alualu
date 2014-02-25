#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

int pn_delete_search( string & bam_input, string &bai_input, vector<string> &chrns, string &path_output, string &pn, string &file_dist_prefix, string &pdf_param, string &file_alupos_prefix, int coverage_max, int alu_flank){  
  ofstream fout((path_output + pn).c_str()); 
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

  map <string, EmpiricalPdf *> empiricalpdf_rg;
  string rg;
  ifstream fin;
  fin.open(( file_dist_prefix + "RG." + pn) .c_str());
  assert(fin);
  while (fin >> rg) 
    empiricalpdf_rg[rg] = new EmpiricalPdf( file_dist_prefix + pn + ".count." + rg + "." + pdf_param );
  fin.close();

  map <string, vector<int> > insertlen_rg;
  float *log_p = new float[3];
  unsigned idx_rg;
  int aluBegin, aluEnd;      

  boost::timer clocki;    
  for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++) {
    string chrx = *ci;
    string file_alupos = file_alupos_prefix + chrx;
    AluRefPos *alurefpos = new AluRefPos(file_alupos);
    int rID = 0;
    if (!getIdByName(nameStore, chrx, rID, nameStoreCache)) {
      cerr << "ERROR: Reference sequence named "<< chrx << " not known.\n";
      return 1;
    }

    for (int count_loci = 0; ; count_loci++) {
      bool hasAlignments = false;
      if (!alurefpos->updatePos(aluBegin, aluEnd)) break;
      if (aluBegin < 0 ) continue;
      if (!jumpToRegion(inStream, hasAlignments, context, rID, aluBegin - alu_flank, aluEnd + alu_flank, baiIndex)) {
	std::cerr << "ERROR: Could not jump to " << aluBegin << ":" << aluEnd << "\n";
	continue;
      }
      if (!hasAlignments) {
	fout << chrx << " " << alu_flank << " " << aluBegin << " " << aluEnd << " 2 2 2 0\n";
	continue;
      }
      
      int reads_cov = 0;    

      insertlen_rg.clear(); 
      clocki.restart();
      while (!atEnd(inStream)) {
	assert (!readRecord(record, context, inStream, seqan::Bam())); 
	if (record.rID != rID || record.beginPos >= aluEnd + alu_flank) break;
	if (record.beginPos < aluBegin - alu_flank) continue;            
	if ((not hasFlagQCNoPass(record) ) and hasFlagAllProper(record) and (not hasFlagDuplicate(record)) and hasFlagMultiple(record) ) {
	  if ( record.beginPos >= record.pNext ) continue;  // ==> hasFlagFirst(record) 	
	  reads_cov ++;
	  if ( (record.pNext <= aluBegin) or (record.beginPos >= aluEnd)) continue;  // ignore broken reads	  
	  seqan::BamTagsDict tags(record.tags);
	  if (!findTagKey(idx_rg, tags, "RG")) continue;
	  rg = toCString(getTagValue(tags, idx_rg));	
	  map <string, vector<int> >::iterator itr;
	  if ( (itr = insertlen_rg.find(rg)) == insertlen_rg.end()) {
	    insertlen_rg[rg].push_back(abs(record.tLen));
	  } else {
	    (itr->second).push_back(abs(record.tLen));
	  }	  
 	}        
      }
      
      float mean_coverage = length(record.seq) * reads_cov * 2. / (aluEnd - aluBegin + alu_flank);
      if (mean_coverage > coverage_max) { 
	fout << chrx << " " << alu_flank << " " << aluBegin << " " << aluEnd << " 3 3 3 " << mean_coverage << endl;
	continue;
      }
      //cerr << endl << count_loci << " time used " << aluBegin - alu_flank << " "<< clocki.elapsed() << " " << reads_cov <<  endl;      
      genotype_prob(insertlen_rg, empiricalpdf_rg, aluEnd - aluBegin, log_p);
      fout << chrx << " " << alu_flank << " " << aluBegin << " " << aluEnd << " " ;
      for (int i = 0; i < 3; i++)  fout << log_p[i] << " " ;
      fout << mean_coverage << endl;
    }
    delete alurefpos;
    cerr << "file_alupos:done  " << file_alupos << endl;  
  }
  delete log_p;
  for (map <string, EmpiricalPdf *>::iterator ri = empiricalpdf_rg.begin(); ri != empiricalpdf_rg.end(); ri++) 
    delete ri->second;
  //bamStreamOut.close();
  fout.close();
  seqan::close(inStream);     
  return 0;
}


int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  string config_file = argv[1];
  int idx_pn;  
  seqan::lexicalCast2(idx_pn, argv[2]);
  string chrn = argv[3];
  string path_output = argv[4];

  int coverage_max, alu_flank;
  string pn, bam_input, bai_input;

  seqan::lexicalCast2(coverage_max, (read_config(config_file, "coverage_max")));
  seqan::lexicalCast2(alu_flank, (read_config(config_file, "alu_flank")));
  
  ifstream fin(read_config(config_file, "file_pn").c_str());
  int i = 0;
  while (fin >> pn) {
    if (i++ == idx_pn ) {
      cerr << "reading pn: " << i << " " << pn << "..................\n";
      /// old bam files in /nfs/gpfs/data/Results/GWS/
      //bam_input = read_config(config_file, "old_bam_prefix") + pn + "/" + pn + ".bam";
      bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
      bai_input = bam_input + ".bai";  
      break;
    }
  }
  fin.close();  

  vector<string> chrns;
  if (chrn != "chr0") {
    chrns.push_back(chrn);
  } else {
    for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
    //for (int i = 21; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  }
  
  string file_dist_prefix = read_config(config_file, "file_dist_prefix");
  string pdf_param = read_config(config_file, "pdf_param");
  string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 

  pn_delete_search(bam_input, bai_input, chrns, path_output, pn, file_dist_prefix, pdf_param, file_alupos_prefix, coverage_max, alu_flank);

  return 0;  
}
