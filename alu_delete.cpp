#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

#define CROSS_BP 5
#define CLIP_BP 5  // min length to be called soft clip 
#define UNKNOW_READ 0
#define MID_READ 1  // read without deletion
#define CLIP_READ 2 

// evidence of deletion
int clip_pos(int beginPos, int endPos, int begin_del, int end_del, seqan::BamAlignmentRecord & record){
  //chr1 6018867 6019167 6018782 6018867 M85S16 8:110:3182:21081
  if ( (abs(endPos-begin_del) < CROSS_BP) and has_soft_last(record, CLIP_BP) ) return CLIP_READ;
  if ( (abs(beginPos-end_del) < CROSS_BP) and has_soft_first(record, CLIP_BP) ) return CLIP_READ;
  return UNKNOW_READ;  // check why so maorny -1 ?
}

// one of the read(in the pair ) covers the delete region
int cross_pos(seqan::BamAlignmentRecord & record, int align_len, int begin_del, int end_del){
  int beginPos = record.beginPos;
  int endPos = record.beginPos + align_len;
  if ( (endPos < begin_del - CROSS_BP) or (beginPos > end_del + CROSS_BP) ) return UNKNOW_READ; 
  if ( (endPos > begin_del + CROSS_BP) and (endPos < end_del - CROSS_BP) ) return MID_READ;
  int endPos_pNext = beginPos + record.tLen;
  if ( (endPos_pNext > begin_del + CROSS_BP) and (endPos_pNext < end_del- CROSS_BP) ) return MID_READ;
  if ( (abs(beginPos-begin_del) < CROSS_BP ) or (abs(endPos-end_del) < CROSS_BP) ) return MID_READ;
  return clip_pos(beginPos, endPos, begin_del, end_del, record);
}

void print_summary(map <seqan::CharString, int> &read_type){
  map<int, int> type_counts;
  for (map <seqan::CharString, int>::iterator rt = read_type.begin(); rt != read_type.end(); rt++) {
    addKey(type_counts, rt->second, 1);
  }
  print_map(type_counts);
}

int chr_alupos(string &chrx, seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, int aluBegin, int aluEnd, int alu_flank, ofstream &f_tmp2, unsigned coverage_max, map <seqan::CharString, int> &read_type, map <seqan::CharString, pair<string, int> > &rg_len) {
  read_type.clear();
  rg_len.clear();
  bool hasAlignments = false;
  if (aluBegin <= alu_flank  or !jumpToRegion(inStream, hasAlignments, context, rID, aluBegin-alu_flank, aluEnd-alu_flank, baiIndex)) return 0;
  if (!hasAlignments) return 0;
  int reads_cov = 0;
  map <seqan::CharString, int>::iterator rt; 
  seqan::BamAlignmentRecord record;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != rID || record.beginPos >= aluEnd + alu_flank) break;
    if (record.beginPos < aluBegin - alu_flank) continue;            
    if (!QC_read(record)) continue;    
    reads_cov ++;  

    int align_len = getAlignmentLengthInRef(record);    
    if ( ( max(record.beginPos + align_len, record.beginPos + record.tLen) < aluBegin - CROSS_BP) or 
	 ( min(record.beginPos, record.pNext) > aluEnd + CROSS_BP) ) {   // useless reads           
      //cerr << "useless reads: " << chrx << " " << aluBegin  << " " << aluEnd << endl;
      continue;        
    }
  
    if ( ( (rt = read_type.find(record.qName)) != read_type.end()) and (rt->second != UNKNOW_READ) ) continue;
    int rt_val = cross_pos( record, align_len, aluBegin, aluEnd );
    read_type[record.qName] = rt_val;

    if ( rt_val == UNKNOW_READ ) {
      seqan::BamTagsDict tags(record.tags);
      unsigned idx_rg;
      assert (findTagKey(idx_rg, tags, "RG"));
      string rg = toCString(getTagValue(tags, idx_rg));	
      rg_len[record.qName] = make_pair(rg, abs(record.tLen));
    }     
    
    if ( rt_val == CLIP_READ ) {
      f_tmp2 << "CLIP_READ " << chrx << " " << aluBegin  << " " << aluEnd << " ";
      print_read(record, f_tmp2);
    }
    
  }  
  
  if (length(record.seq) * reads_cov / (aluEnd - aluBegin + 2 * alu_flank) > coverage_max) return 0;
  print_summary(read_type);   
  return 0;
}
  
int delete_search( string & bam_input, string &bai_input, vector<string> &chrns, string &path_output, string &pn, string &file_dist_prefix, string &pdf_param, string &file_alupos_prefix, int coverage_max, int alu_flank) {  

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
  EmpiricalPdf *empiricalpdf;
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
  map <seqan::CharString, int> read_type; 
  map <seqan::CharString, int>::iterator rt;
  map <seqan::CharString, pair<string, int> > rg_len; 

  ofstream f_tmp1((path_output + pn + ".tmp1").c_str()); 
  f_tmp1 << "chr alu_flank aluBegin aluEnd p0 p1 p2 mean_coverage\n";
  ofstream f_tmp2((path_output + pn + ".tmp2").c_str()); 
  f_tmp2 << "type chr aluBegin aluEnd qName beginPos endPos pNext pNextEnd\n";
  
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
      if (!alurefpos->updatePos(aluBegin, aluEnd)) break;
      int int_err = chr_alupos(chrx, inStream, baiIndex, context, rID, aluBegin, aluEnd, alu_flank, f_tmp2, coverage_max, read_type, rg_len);
      if (int_err < 0) continue;
      int mid_read_count = 0;
      for (int i = 0; i < 3; i++) log_p[i] = 0;
      for ( rt = read_type.begin(); rt != read_type.end(); rt++) {
	if (rt->second == MID_READ) { 
	  mid_read_count++;
	} else if (rt->second == UNKNOW_READ) {
	  empiricalpdf = empiricalpdf_rg[rg_len[rt->first].first];
	  int insert_len = rg_len[rt->first].second;
	  float p_y = empiricalpdf->pdf_obs(insert_len);
	  float p_z = empiricalpdf->pdf_obs(insert_len + aluEnd - aluBegin );
	  log_p[0] += log(p_y);
	  log_p[1] += log(0.67 * p_y + 0.33 * p_z);
	  log_p[2] += log(p_z); 
	}	
      }
      
      f_tmp1 << chrx << " " << alu_flank << " " << aluBegin << " " << aluEnd << " " ;
      for (int i = 0; i < 3; i++)  f_tmp1 << log_p[i] << " " ;
      f_tmp1 << mean_coverage << endl;
      //cerr << endl << count_loci << " time used " << aluBegin - alu_flank << " "<< clocki.elapsed() << " " << reads_cov <<  endl;      
    }
    delete alurefpos;
    cerr << "file_alupos:done  " << file_alupos << endl;  
  }
  delete log_p;
  for (map <string, EmpiricalPdf *>::iterator ri = empiricalpdf_rg.begin(); ri != empiricalpdf_rg.end(); ri++) 
    delete ri->second;
  //bamStreamOut.close();
  f_tmp1.close();
  f_tmp2.close();
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

  unsigned coverage_max, alu_flank;
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

  delete_search(bam_input, bai_input, chrns, path_output, pn, file_dist_prefix, pdf_param, file_alupos_prefix, coverage_max, alu_flank);
  
  
  return 0;  
}
