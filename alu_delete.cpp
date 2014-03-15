#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

#define LOG_RATIO_UB 7 // max log(Odds Ratio) for a single read, estimated from quantile(0.95)/quantile(0.05), 1e3 to 1e4
#define CROSS_BP 5
#define CLIP_BP 10  // min length to be called soft clip 
#define BAD_READ -1  // not use! eg. one read is mid-read, the other is clip read
#define UNKNOW_READ 0
#define MID_READ 1  // read without deletion
#define CLIP_READ 2 
typedef map<int, int> MapII;
typedef map<int, int>::iterator MapIIt;


inline string get_name1(string &path, string &fn){
  return path + fn + ".tmp1";
}

inline string get_name2(string &path, string &fn){
  return path + fn + ".tmp2";
}

// evidence of deletion
int clip_pos(int beginPos, int endPos, int begin_del, int end_del, seqan::BamAlignmentRecord & record){
  //chr1 6018867 6019167 6018782 6018867 M85S16 8:110:3182:21081
  if ( (abs(endPos-begin_del) < CROSS_BP) and has_soft_last(record, CLIP_BP) ) return CLIP_READ;
  if ( (abs(beginPos-end_del) < CROSS_BP) and has_soft_first(record, CLIP_BP) ) return CLIP_READ;
  return UNKNOW_READ;  // check why so many unknown ?
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

void count_reads(map <seqan::CharString, int> &special_read, int &mid_read_count, int &clip_read_count){
  for (map <seqan::CharString, int>::iterator rt = special_read.begin(); rt != special_read.end(); rt++) 
    if (rt->second == MID_READ) mid_read_count++ ;
    else clip_read_count++;
}
 
int chr_alupos(string &chrx, seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, int aluBegin, int aluEnd, int alu_flank, ofstream &f_tmp2, unsigned coverage_max, unsigned &coverage_mean, map <seqan::CharString, int> &special_read, map <seqan::CharString, pair<string, int> > &rg_len){
  bool hasAlignments = false;
  if (aluBegin <= alu_flank  or !jumpToRegion(inStream, hasAlignments, context, rID, aluBegin-alu_flank, aluEnd-alu_flank, baiIndex)) return 0;
  if (!hasAlignments) return 0;
  int reads_cov = 0;
  seqan::BamAlignmentRecord record;
  rg_len.clear();
  special_read.clear();
  map <seqan::CharString, int>::iterator sr;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != rID || record.beginPos >= aluEnd + alu_flank) break;
    if (record.beginPos < aluBegin - alu_flank) continue;            
    if ( (!QC_read(record)) or (abs(record.tLen) > 2000)) continue;    // ignore extreme large tLen 
    reads_cov ++;  // only left read counts
    int align_len = getAlignmentLengthInRef(record);    
    if ( ( max(record.beginPos + align_len, record.beginPos + record.tLen) < aluBegin - CROSS_BP) or 
	 ( min(record.beginPos, record.pNext) > aluEnd + CROSS_BP) ) {   // useless reads           
      continue;        
    }

    if ( special_read.find(record.qName) != special_read.end()) continue;    
    int rt_val = cross_pos( record, align_len, aluBegin, aluEnd );
    if (rt_val != UNKNOW_READ) {
      special_read[record.qName] = rt_val;
      if (rt_val == CLIP_READ) {
	f_tmp2 << "CLIP_READ " << chrx << " " << aluBegin  << " " << aluEnd << " ";
	print_read(record, f_tmp2);	
      }
    } else if (left_read(record)) {
      seqan::BamTagsDict tags(record.tags);
      unsigned idx_rg;
      assert (findTagKey(idx_rg, tags, "RG"));
      string rg = toCString(getTagValue(tags, idx_rg));	
      rg_len[record.qName] = make_pair(rg, abs(record.tLen));      
    }    
  }   
  coverage_mean = length(record.seq) * reads_cov / (aluEnd - aluBegin + 2 * alu_flank) ;
  if ( coverage_mean > coverage_max) return 0;
  return 1;
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
  map <seqan::CharString, int> special_read;
  map <seqan::CharString, pair<string, int> > rg_len; 
  map <seqan::CharString, pair<string, int> >::iterator rl; 
  
  ofstream f_tmp1( get_name1(path_output, pn).c_str()); 
  f_tmp1 << "chr alu_flank aluBegin aluEnd p0 p1 p2 mean_coverage mid_read_count clip_read_count unknow_read_count\n";
  ofstream f_tmp2( get_name2(path_output, pn).c_str()); 
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
      unsigned coverage_mean = 0;
      if (!alurefpos->updatePos(aluBegin, aluEnd)) break;
      if (!chr_alupos(chrx, inStream, baiIndex, context, rID, aluBegin, aluEnd, alu_flank, f_tmp2, coverage_max, coverage_mean, special_read, rg_len)) continue;
      int mid_read_count = 0, clip_read_count = 0, unknow_read_count = 0;
      count_reads(special_read, mid_read_count, clip_read_count);      
      
      // cerr << "count " << mid_read_count << " " <<  clip_read_count << endl;
      // clip_read: evidence for deletion, mid_read: evidence for insertion

      log_p[0] = clip_read_count * (-LOG_RATIO_UB) ; 
      log_p[1] = (mid_read_count + clip_read_count) * log(0.5) ;
      log_p[2] = mid_read_count * (-LOG_RATIO_UB) ; 
      
      for ( rl = rg_len.begin(); rl != rg_len.end(); rl++) {
	if (special_read.find(rl->first) != special_read.end()) continue;
	unknow_read_count ++;
	empiricalpdf = empiricalpdf_rg[(rl->second).first];
	int insert_len = (rl->second).second;
	float p_y = empiricalpdf->pdf_obs(insert_len);
	float p_z = empiricalpdf->pdf_obs(insert_len + aluEnd - aluBegin );
	log_p[0] += log(p_y);
	log_p[1] += log(0.67 * p_y + 0.33 * p_z);
	log_p[2] += log(p_z); 
      }	
      if ( (!unknow_read_count) and (!clip_read_count)) continue;
      if ( !normalize_prob(log_p) ) continue;
      f_tmp1 << chrx << " " << alu_flank << " " << aluBegin << " " << aluEnd << " " ;
      for (int i = 0; i < 3; i++)  f_tmp1 << log_p[i] << " " ;
      f_tmp1 << coverage_mean << " "<< mid_read_count << " " << clip_read_count << " " <<  unknow_read_count << endl;
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

bool combine_pns(vector <string> &pns, string path1, string f_out){
  ofstream fout( f_out.c_str() );
  map <string, MapII > del_count;
  ifstream fin;
  stringstream ss;
  string line, chrn, tmp1, tmp2;
  int aluBegin, aluEnd;  
  float p0, p1, p2;
  for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
    fin.open( get_name1(path1, *pi).c_str() );
    assert(fin);
    while ( getline(fin, line) ) {
      ss.clear(); ss.str( line );
      ss >> chrn >> tmp1 >> aluBegin >> aluEnd >> p0 >> p1 >> p2 ;
      if (p1 > p0 or p2 > p0) addKey(del_count[chrn], aluBegin, p1 > p2 ? 1 : 2); 
    }
    fin.close();
  }
  float pn_chr = 2. * pns.size();
  map<string, MapII >::iterator dc;
  for (  dc = del_count.begin(); dc != del_count.end(); dc++ )
    for (MapIIt dc2 = (dc->second).begin(); dc2 != (dc->second).end(); dc2++)
      if (dc2->second > 1) fout << dc->first << " " << dc2->first << " " << dc2->second << endl;
  fout.close();
  return true;
}

void check_clip(string bam_input, string bai_input, string path1, string pn){
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::Stream<seqan::Bgzf> inStream;  
  open(inStream, bam_input.c_str(), "r");
  seqan::BamIndex<seqan::Bai> baiIndex;
  assert ( !read(baiIndex, bai_input.c_str())) ;
  assert ( !readRecord(header, context, inStream, seqan::Bam()) ); // do i need this ??

  seqan::BamAlignmentRecord record;    
  string tmp1, chrn, qName, cigar;
  int rID, aluBegin, aluEnd, beginPos, endPos, pNext, pNextEnd;
  string chrn_pre = "chr-1";
  int ir = 0;
  ifstream fin( get_name2(path1, pn).c_str() );
  assert(fin);
  getline(fin, tmp1 );
  while ( fin >> tmp1 >> chrn >> aluBegin >> aluEnd >> qName >> beginPos >> endPos >> pNext >> pNextEnd >> cigar ) {
    if (chrn != chrn_pre) {
      assert ( getIdByName(nameStore, chrn, rID, nameStoreCache) );
      chrn_pre = chrn;
    }
    bool hasAlignments = false;      
    jumpToRegion(inStream, hasAlignments, context, rID, beginPos - 5, beginPos + 5, baiIndex);
    if (!hasAlignments) {
      continue;
    }
    if (find_read(inStream, context, rID, beginPos, qName, record)) {
      cerr << qName << " " << record.seq << endl;
      if (ir++ > 2) exit(0);
    }
  }
  fin.close();
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);

  int opt, idx_pn;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];
  string path1 = read_config(config_file, "file_delete_pn_prefix");
  string file_pn = read_config(config_file, "file_pn");
  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  if ( opt == 1 ) {  // write potential candidate
    seqan::lexicalCast2(idx_pn, argv[3]);
    string chrn = argv[4];
    string pn = get_pn(file_pn, idx_pn);
    unsigned coverage_max, alu_flank;
    seqan::lexicalCast2(coverage_max, (read_config(config_file, "coverage_max")));
    seqan::lexicalCast2(alu_flank, (read_config(config_file, "alu_flank")));
    cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";  

    if (chrn != "chr0") {  
      chrns.clear();
      chrns.push_back(chrn);
    } 
    string file_dist_prefix = read_config(config_file, "file_dist_prefix");
    string pdf_param = read_config(config_file, "pdf_param");
    string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
    delete_search(bam_input, bai_input, chrns, path1, pn, file_dist_prefix, pdf_param, file_alupos_prefix, coverage_max, alu_flank);
  } else if (opt == 2) {   // filter candidate
    vector <string> pns;
    ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
    while (fin >> idx_pn) pns.push_back(get_pn(file_pn, idx_pn));
    fin.close();
    //print_vec(pns);
    string f_out = path1+"test_pn";
    combine_pns(pns, path1, f_out); 
    cerr << "output to " << f_out << endl;
  }  else if (opt == 3) {   // check clip_read
    seqan::lexicalCast2(idx_pn, argv[3]);
    string pn = get_pn(file_pn, idx_pn);
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";  
    check_clip(bam_input, bai_input, path1, pn);
  }
  return 0;  
}
