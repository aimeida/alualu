#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

#define CROSS_BP 5
#define PNEXT_BP 10 // pNext and breakpoint
#define CLIP_BP 5  // min length to be called soft clip 
#define BAD_READ -1
#define UNKNOW_READ 0
#define CROSS_READ 1  // read without deletion
#define CLIP_READ 2 

int cross_pos(seqan::BamAlignmentRecord & record, int align_len, int begin_del, int end_del){
  if (( record.beginPos <= begin_del - CROSS_BP) and ( record.beginPos + align_len >= begin_del + CROSS_BP )) return CROSS_READ;
  if (( record.beginPos <= end_del - CROSS_BP) and ( record.beginPos + align_len >= end_del + CROSS_BP )) return CROSS_READ;
  if (( record.pNext <= begin_del - CROSS_BP) and ( record.pNext + PNEXT_BP > begin_del)) return CROSS_READ;
  if (( record.pNext <= end_del - CROSS_BP) and ( record.pNext + PNEXT_BP > end_del)) return CROSS_READ;
  return UNKNOW_READ; // check ?
}

int clip_pos(int beginPos, int endPos, int begin_del, int end_del, seqan::BamAlignmentRecord & record){
  if ( abs(beginPos-begin_del) < CROSS_BP ) return CROSS_READ;
  if ( abs(endPos-end_del) < CROSS_BP ) return CROSS_READ; 
  //chr1 6018867 6019167 6018782 6018867 M85S16 8:110:3182:21081
  if ( (abs(endPos-begin_del) < CROSS_BP) and has_soft_last(record, CLIP_BP) ) return CLIP_READ;
  if ( (abs(beginPos-end_del) < CROSS_BP) and has_soft_first(record, CLIP_BP) ) return CLIP_READ;
  return UNKNOW_READ;  // check why so maorny -1 ?
}

int check_region(string &chrx, seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, int aluBegin, int aluEnd, int alu_flank, map <string, vector<int> > &insertlen_rg, ofstream &fout, int coverage_max) {
  bool hasAlignments = false;
  if (aluBegin <= alu_flank  or !jumpToRegion(inStream, hasAlignments, context, rID, aluBegin-alu_flank, aluEnd-alu_flank, baiIndex)) return 0;
  if (!hasAlignments) return 0;
  int reads_cov = 0;
  set <seqan::CharString> med_del;  // reads with deletion, 0.999 sure
  set <seqan::CharString> uden_del; // reads without deletion, 0.999 sure
  seqan::BamAlignmentRecord record;
  unsigned idx_rg;
  map <string, vector<int> >::iterator itr;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != rID || record.beginPos >= aluEnd + alu_flank) break;
    if (record.beginPos < aluBegin - alu_flank) continue;            
    if (!QC_read(record)) continue;    
    if (left_read(record)) reads_cov ++;  

    if ( (med_del.find(record.qName)!=med_del.end()) or (uden_del.find(record.qName)!=uden_del.end()) ) continue;    
    int read_len = length(record.seq);
    int align_len = getAlignmentLengthInRef(record);
    if ( cross_pos( record, align_len, aluBegin, aluEnd) == CROSS_READ) {
      uden_del.insert(record.qName);
      continue;
    }
    int status_del = clip_pos( record.beginPos, record.beginPos + align_len, aluBegin, aluEnd, record);
    if (status_del == CROSS_READ) {
      uden_del.insert(record.qName);
      continue;      
    }

    //    chr1 7960397 7960634 7960635 7960716 S20M81 7:55:17103:13566 7960461
    if (status_del == CLIP_READ) {
      med_del.insert(record.qName);
      cerr << "clip " << " " <<  chrx << " " << aluBegin << " " << aluEnd << " " << record.beginPos << " " << record.beginPos + getAlignmentLengthInRef(record) << " ";
      print_cigar(record);
      cerr << " " << record.qName << " " << record.pNext << "\n";                    
      continue;      
    }
  
    if (status_del == BAD_READ) {
      cerr << "clip " << " " <<  chrx << " " << aluBegin << " " << aluEnd << " " << record.beginPos << " " << record.beginPos + getAlignmentLengthInRef(record) << " ";
      print_cigar(record);
      cerr << " " << record.qName << " " << record.pNext << "\n";                    
    } 
    
  }
  // first check if it cross 2 regions.
       
//    if (length(record.cigar)==0) 
//      cerr << chrx << " " << aluBegin << " " << aluEnd << " " << record.qName << " " << record.beginPos << " " << record.pNext << " " << length(record.cigar) << "\n";    
//if (record.cigar[0].operation == 'S' or record.cigar[len_cigar-1].operation =='S') {
  
  //cerr << "cigar " << record.cigar[0].operation << endl;

  /*
    if ( (record.pNext <= aluBegin) or (record.beginPos >= aluEnd)) continue;  // ignore broken reads	  
    seqan::BamTagsDict tags(record.tags);
    assert (findTagKey(idx_rg, tags, "RG"));
    string rg = toCString(getTagValue(tags, idx_rg));	
    if ( (itr = insertlen_rg.find(rg)) == insertlen_rg.end()) {
      insertlen_rg[rg].push_back(abs(record.tLen));
    } else {
	(itr->second).push_back(abs(record.tLen));
    }	  
    */  
  
  //if (reads_cov > 2) cerr << " ## " << reads_cov << " " << reads_mid << endl; 
  //float mean_coverage = length(record.seq) * reads_cov * 2. / (aluEnd - aluBegin + alu_flank);
  //if ( (reads_cov <= 2) or (mean_coverage > coverage_max) ) return 0;
  return 0;
}
  
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

  float *log_p = new float[3];
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
    map <string, vector<int> > insertlen_rg;
    for (int count_loci = 0; ; count_loci++) {
      if (!alurefpos->updatePos(aluBegin, aluEnd)) break;
      clocki.restart();
      insertlen_rg.clear();
      int int_err = check_region(chrx, inStream, baiIndex, context, rID, aluBegin, aluEnd, alu_flank, insertlen_rg, fout, coverage_max);
      if (int_err < 0) cerr << chrx << " " <<  aluBegin << " " <<  aluEnd << "\n";


      /*
      genotype_prob(insertlen_rg, empiricalpdf_rg, aluEnd - aluBegin, log_p);
      fout << chrx << " " << alu_flank << " " << aluBegin << " " << aluEnd << " " ;
      for (int i = 0; i < 3; i++)  fout << log_p[i] << " " ;
      fout << mean_coverage << endl;
      */
      //cerr << endl << count_loci << " time used " << aluBegin - alu_flank << " "<< clocki.elapsed() << " " << reads_cov <<  endl;      
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
  seqan::lexicalCast2(coverage_max, (read_config(config_file, "coverage_max")));
  seqan::lexicalCast2(alu_flank, (read_config(config_file, "alu_flank")));
  
  string pn = get_pn(read_config(config_file, "file_pn"), idx_pn);
  cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
  /// old bam files in /nfs/gpfs/data/Results/GWS/
  //bam_input = read_config(config_file, "old_bam_prefix") + pn + "/" + pn + ".bam";
  string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
  string bai_input = bam_input + ".bai";  

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
