#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/vcf_io.h>

#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

#define LOG_RATIO_UB 7 // max log(Odds Ratio) for a single read, estimated from quantile(0.95)/quantile(0.05), 1e3 to 1e4
#define MIN_CLIP_ALIGN_SCORE 5
#define CROSS_BP 5
#define CLIP_BP 10  // min length to be called soft clip 
#define CLIP_ALIGN 10 // extend bp in both end for split mapping 
#define BAD_READ -1  // not use! eg. one read is mid-read, the other is clip read
#define UNKNOW_READ 0  // we use them for length estimation
#define MID_READ 1  // read without deletion
#define CLIP_READ 2 

typedef map<int, int> MapII;
typedef map<int, int>::iterator MapIIt;

inline string get_name1(string &path, string &fn){ return path + fn + ".tmp1";}
inline string get_name2(string &path, string &fn){ return path + fn + ".tmp2";}
inline string get_name3(string &path, string &fn){ return path + fn + ".tmp3";}


bool get_align_pos(int aluBegin, int aluEnd, int beginPos, int endPos, int &ref_a, int &ref_b, int &read_a, int &read_b, seqan::BamAlignmentRecord &record){
  int len_match = 0;
  unsigned i = length(record.cigar) - 1; // ignore complicated alignment between S and M
  int len_read = length(record.seq);
  if ( abs(beginPos - aluEnd) < CROSS_BP and record.cigar[i].operation == 'M') {  //eg. S41M60, has_soft_first 
    len_match = record.cigar[i].count;
    ref_a = aluBegin - (len_read - len_match);
    ref_b = aluBegin;
    read_a = 0;
    read_b = len_read - len_match;
  } else if (abs(endPos - aluBegin) < CROSS_BP and record.cigar[0].operation == 'M' ) { // eg. M27S93, has_soft_last
    len_match = record.cigar[0].count;
    ref_a = aluEnd;
    ref_b = aluEnd + (len_read - len_match);
    read_a = len_match;
    read_b = len_read;
  }
  return len_match ? true : false;
}

bool split_global_align(seqan::CharString &fa_seq, seqan::BamAlignmentRecord &record, int read_a, int read_b, int &score, int &align_len){
  // end quality if bad!!! use external mapping tool ?
  // no need to use seed or chain, since don't know exact position to start mapping
  read_b = read_a + (read_b - read_a) * 0.75 + 2; // end might have bad quality
  seqan::Score<int> scoringScheme(1, -3, -2, -5);
  TAlign align;
  int align_len_min = min(CLIP_ALIGN, read_b - read_a );
  resize(rows(align), 2);
  assignSource(row(align,0), infix(record.seq, read_a, read_b) );  // 2,3
  assignSource(row(align,1), fa_seq);   // 1,4
  // AlignConfig[1,2,3,4], 1=true: targetsite duplication ? 2=true:free end gaps for record.seq
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<false, true, true, false>()); 
  align_len =  toViewPosition(row(align, 0), read_b - read_a) - max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  if (score >= MIN_CLIP_ALIGN_SCORE and align_len >= align_len_min) {
    //cerr << "score1 is " << score << " " << align_len << "\n" << align << endl;
    return true;
  }
  score = globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, false>()); 
  align_len = toViewPosition(row(align, 0), read_b - read_a) - max(toViewPosition(row(align, 0), 0), toViewPosition(row(align, 1), 0));
  return (score >= MIN_CLIP_ALIGN_SCORE and align_len >= align_len_min); 
  ////TPosition prefixLength = toSourcePosition(row(align, 1), toViewPosition(row(align, 0), 0));  
  ////cerr << clippedBeginPosition(row(align, 0)) << " " << clippedBeginPosition(row(align, 1)) << " " <<  endl;
}

// evidence of deletion
int clip_pos(int beginPos, int endPos, int aluBegin, int aluEnd, seqan::BamAlignmentRecord & record, seqan::FaiIndex &faiIndex, unsigned fa_idx){
  //chr1 6018867 6019167 6018782 6018867 M85S16 8:110:3182:21081
  int ref_a, ref_b, read_a, read_b, score, align_len;
  seqan::CharString fa_seq;
  if ( (abs(endPos-aluBegin) < CROSS_BP and has_soft_last(record, CLIP_BP) ) or   // maybe clip_read
       (abs(beginPos-aluEnd) < CROSS_BP and has_soft_first(record, CLIP_BP)) ) {   
    if ( get_align_pos(aluBegin, aluEnd, beginPos, endPos, ref_a, ref_b, read_a, read_b, record)) {
      fa_seq = fasta_seq(faiIndex, fa_idx, ref_a - CLIP_ALIGN, ref_b + CLIP_ALIGN, true);      
      if (split_global_align(fa_seq, record, read_a, read_b, score, align_len)) return CLIP_READ;
    }
    if ( max(record.beginPos, record.pNext) < aluEnd - length(record.seq) or min(record.beginPos, record.pNext) > aluBegin) 
      return BAD_READ;
  }
  return UNKNOW_READ;  // check why so many unknown ?
}

// one of the read(in the pair ) covers the delete region
int cross_pos(seqan::BamAlignmentRecord & record, int align_len, int aluBegin, int aluEnd, seqan::FaiIndex &faiIndex, unsigned fa_idx){
  int beginPos = record.beginPos;
  int endPos = record.beginPos + align_len;
  if ( (endPos < aluBegin - CROSS_BP) or (beginPos > aluEnd + CROSS_BP) ) return UNKNOW_READ; 
  if ( (endPos > aluBegin + CROSS_BP) and (endPos < aluEnd - CROSS_BP) ) return MID_READ;
  int endPos_pNext = beginPos + record.tLen;
  if ( (endPos_pNext > aluBegin + CROSS_BP) and (endPos_pNext < aluEnd- CROSS_BP) ) return MID_READ;
  if ( (abs(beginPos-aluBegin) < CROSS_BP ) or (abs(endPos-aluEnd) < CROSS_BP) ) return MID_READ;
  return clip_pos(beginPos, endPos, aluBegin, aluEnd, record, faiIndex, fa_idx);
}

void count_reads(map <seqan::CharString, int> &special_read, int &mid_read_count, int &clip_read_count){
  for (map <seqan::CharString, int>::iterator rt = special_read.begin(); rt != special_read.end(); rt++) 
    if (rt->second == MID_READ) mid_read_count++ ;
    else if (rt->second == CLIP_READ) clip_read_count++;
}
 
int chr_alupos(seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, seqan::FaiIndex &faiIndex, unsigned fa_idx, TBamIOContext &context, int rID, int aluBegin, int aluEnd, int alu_flank, unsigned coverage_max, unsigned &coverage_mean, map <seqan::CharString, int> &special_read, map <seqan::CharString, pair<string, int> > &rg_len){
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
    int rt_val = cross_pos( record, align_len, aluBegin, aluEnd, faiIndex, fa_idx);
    if (rt_val != UNKNOW_READ) {
      if (rt_val != BAD_READ) special_read[record.qName] = rt_val;
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
 
int delete_search( string & bam_input, string &bai_input, string file_fa_prefix, vector<string> &chrns, string &path_output, string &pn, string &file_dist_prefix, string &pdf_param, string &file_alupos_prefix, int coverage_max, int alu_flank) {  

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
  assert(!readRecord(header, context, inStream, seqan::Bam()));
  seqan::FaiIndex faiIndex;
  unsigned fa_idx = 0;

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
  //ofstream f_tmp2( get_name2(path_output, pn).c_str()); 
  //f_tmp2 << "type chr aluBegin aluEnd qName beginPos endPos pNext pNextEnd\n";
  
  for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++) {
    string chrn = *ci;
    string file_alupos = file_alupos_prefix + chrn;
    AluRefPosRead *alurefpos = new AluRefPosRead(file_alupos, 200);
    int rID = 0;
    if (!getIdByName(nameStore, chrn, rID, nameStoreCache)) {
      cerr << "ERROR: Reference sequence named "<< chrn << " not known.\n";
      return 1;
    }
    
    assert (!read(faiIndex, (file_fa_prefix + chrn + ".fa").c_str()) );      
    assert (getIdByName(faiIndex, chrn, fa_idx));

    for (int count_loci = 0; ; count_loci++) {
      unsigned coverage_mean = 0;
      if (!alurefpos->updatePos(aluBegin, aluEnd)) break;
      if (!chr_alupos(inStream, baiIndex, faiIndex, fa_idx, context, rID, aluBegin, aluEnd, alu_flank, coverage_max, coverage_mean, special_read, rg_len)) continue;
      int mid_read_count = 0, clip_read_count = 0, unknow_read_count = 0;
      count_reads(special_read, mid_read_count, clip_read_count);      
      ///  other reads in rg_len

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
      f_tmp1 << chrn << " " << alu_flank << " " << aluBegin << " " << aluEnd << " " ;
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
  seqan::close(inStream);     
  return 0;
}

bool combine_pns_count(vector <string> &pns, string path1, string f_out){
  map <string, MapII > del_count;
  map <string, MapII >::iterator dc;
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
  //float pn_chr = 2. * pns.size();
  ofstream fout( f_out.c_str() );
  for (  dc = del_count.begin(); dc != del_count.end(); dc++ )
    for (MapIIt dc2 = (dc->second).begin(); dc2 != (dc->second).end(); dc2++)
      if (dc2->second > 1) fout << dc->first << " " << dc2->first << " " << dc2->second << endl;
  fout.close();
  return true;
}

string phred_log (float p) {
  if (p) return int_to_string ( -(int)(log10 (p) * 10));
  else return "255";
}

string phred_scaled(float p0, float p1, float p2){
  return phred_log(p0)  + "," + phred_log(p1) + "," + phred_log(p2);
}

void combine_pns_vcf(vector <string> &pns, string path1, string f_out, vector <string> &chrns){
  ifstream fin;
  stringstream ss;
  string line, chrn, tmp1, tmp2;
  int aluBegin, aluEnd, cii, flag;

  map < pair<int, int>, map<string, string> > pos_pn_prob;
  float p0, p1, p2;      
  map < pair<int, int>, map<string, string> >::iterator ppp;
  map<string, string>::iterator pp;
  vector <string>::iterator ci, pi;
  //seqan::VcfStream out("-", seqan::VcfStream::WRITE);
  seqan::VcfStream vcfout(seqan::toCString(f_out), seqan::VcfStream::WRITE);
  for ( ci = chrns.begin(); ci != chrns.end(); ci++)
    appendValue(vcfout.header.sequenceNames, *ci);
  for ( pi = pns.begin(); pi != pns.end(); pi++) 
    appendValue(vcfout.header.sampleNames, *pi);
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("fileformat", "VCFv4.1"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("fileDate", "201402"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("reference", "hg18"));
  appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FILTER", "<ID=PL,Number=3,Type=Integer, Description=\"Phred-scaled likelihoods for genotypes\">"));
  seqan::VcfRecord record;    
  record.ref = ".";
  record.alt = "1";
  record.qual = 0;
  record.filter = ".";
  record.info = ".";
  record.format = ".";
    
  for ( cii = 0, ci = chrns.begin(); ci != chrns.end(); ci++, cii++) {
    pos_pn_prob.clear();
    for ( pi = pns.begin(); pi != pns.end(); pi++) {
      fin.open( get_name1(path1, *pi).c_str() );
      assert(fin);
      flag = 1;
      while ( getline(fin, line) ) {
	ss.clear(); ss.str( line );
	ss >> chrn >> tmp1 >> aluBegin >> aluEnd >> p0 >> p1 >> p2 ;
	if (chrn != *ci) {
	  if (flag) continue;
	  else break;
	}
	flag = 0;
	if (p1 > p0 or p2 > p0) pos_pn_prob[ make_pair(aluBegin, aluEnd) ][*pi] = phred_scaled(p0, p1, p2);
      }
      fin.close();
    }

    // Write out the records.
    for (ppp = pos_pn_prob.begin(); ppp != pos_pn_prob.end(); ppp++) {
      record.beginPos = (ppp->first).first;
      record.rID = cii;
      record.id = int_to_string((ppp->first).second);
      for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
	if ( (pp = ppp->second.find(*pi)) != ppp->second.end() )  appendValue(record.genotypeInfos, pp->second);
	else  appendValue(record.genotypeInfos, "0,255,255"); 	  
      }
      writeRecord(vcfout, record);
      clear(record.genotypeInfos);
    }        
  }  // chrn finished
  clear(record);
  seqan::close(vcfout);
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
    string file_fa_prefix = read_config(config_file, "file_fa_prefix");

    if (chrn != "chr0") {  
      chrns.clear();
      chrns.push_back(chrn);
    } 
    string file_dist_prefix = read_config(config_file, "file_dist_prefix");
    string pdf_param = read_config(config_file, "pdf_param");
    string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
    delete_search(bam_input, bai_input, file_fa_prefix, chrns, path1, pn, file_dist_prefix, pdf_param, file_alupos_prefix, coverage_max, alu_flank);
  } else if (opt == 2) {   // filter candidate
    vector <string> pns;
    ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
    while (fin >> idx_pn) pns.push_back(get_pn(file_pn, idx_pn));
    fin.close();
    //print_vec(pns);
    string f_out = path1+"test_pn";
    //combine_pns_count(pns, path1, f_out); 
    f_out = path1+"test_vcf";
    combine_pns_vcf(pns, path1, f_out, chrns);
    cerr << "output to " << f_out << endl;
  } 


return 0;  
}
