#define SEQAN_HAS_ZLIB 1
#include <seqan/vcf_io.h>
#include "common.h"
#include "delete_utils.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

inline string get_name_tmp(string path, string fn, string suffix){ return path + fn + suffix;}
inline string get_name_rg(string prefix, string pn){ return prefix + "RG." + pn;}
inline string get_name_rg_pdf(string prefix, string pn, string rg, string pdf_param){ return prefix + pn + ".count." + rg + "." + pdf_param; }

void count_reads(map <seqan::CharString, T_READ> &special_read, map < T_READ, int > &readCnt) {
  readCnt.clear();
  for (map <seqan::CharString, T_READ>::iterator rt = special_read.begin(); rt != special_read.end(); rt++) 
    addKey(readCnt, rt->second); 
}
 
int check_chr_alupos(seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, seqan::FaiIndex &faiIndex, unsigned fa_idx, map <string, int> &rg_to_idx, TBamIOContext &context, int rID, int aluBegin, int aluEnd, unsigned coverage_max, float &coverage_mean, map <seqan::CharString, T_READ> &special_read,  string &rg_str){
  
  special_read.clear();
  stringstream rg_ss;

  bool hasAlignments = false;
  if (aluBegin <= ALU_FLANK  or !jumpToRegion(inStream, hasAlignments, context, rID, aluBegin-ALU_FLANK, aluEnd-ALU_FLANK, baiIndex)) return 0;
  if (!hasAlignments) return 0;
  int reads_cnt = 0;
  T_READ rt_val;
  seqan::BamAlignmentRecord record;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    //  flank_bound -- alu_begin -- alu_end --(begin)-- flank_bound
    if (record.rID != rID || record.beginPos >= aluEnd + ALU_FLANK) break;
    if (record.beginPos + DEFAULT_READ_LEN < aluBegin - ALU_FLANK) continue;
    if ( !QC_delete_read(record) )  continue;    // ignore extreme large tLen 
    
    reads_cnt ++;  
    int align_len = getAlignmentLengthInRef(record);    
    bool read_is_left = left_read(record);
    if ( (!read_is_left) and special_read.find(record.qName) != special_read.end() ) continue;    	
    rt_val = classify_read( record, align_len, aluBegin, aluEnd, faiIndex, fa_idx, read_is_left);
    
    if ( rt_val != useless_read)  
      special_read[record.qName] = rt_val;
    if ( rt_val == unknow_read ) { // (unknow = jump across alu region)
	seqan::BamTagsDict tags(record.tags);
	unsigned idx_rg;
	assert (findTagKey(idx_rg, tags, "RG"));
	rg_ss << rg_to_idx[seqan::toCString(getTagValue(tags, idx_rg))] << ":" << abs(record.tLen) << " ";
    }	    
    /*    if ( rt_val == clip_read ) 
	  if ( rt_val == useless_read and (!read_is_left) )  
	  cout << "useless " << record.beginPos << " " << get_cigar(record) << " " << record.pNext << " " << record.tLen << endl;
    */
  }
  coverage_mean = length(record.seq) * (float) reads_cnt / (aluEnd - aluBegin + 2 * ALU_FLANK) ;
  if ( coverage_mean > coverage_max) return COVERAGE_HIGH;
  rg_str = rg_ss.str();
  return 1;
}

int delete_search( string & bam_input, string &bai_input, string file_fa_prefix, vector<string> &chrns, string &f_out, string &f_log, string &file_alupos_prefix, int coverage_max, map<string, int> &rg_to_idx) {    
  // Open BGZF Stream for reading.
  seqan::Stream<seqan::Bgzf> inStream;
  if (!open(inStream, bam_input.c_str(), "r")) {
    std::cerr << "ERROR: Could not open " << bam_input << " for reading.\n";
    return 1;
  }  
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
  int aluBegin, aluEnd;      
  map < seqan::CharString, T_READ> special_read;  
  ofstream f_tmp1( f_out.c_str()); 
  f_tmp1 << "chr aluBegin aluEnd mean_coverage midCnt clipCnt unknowCnt unknowStr\n";
  ofstream f_log1( f_log.c_str());  // print out info for clip reads 

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
      float coverage_mean = 0;
      if (! alurefpos->updatePos(aluBegin, aluEnd)) break;      
      string rg_str;
      int check_alu = check_chr_alupos(inStream, baiIndex, faiIndex, fa_idx, rg_to_idx, context, rID, aluBegin, aluEnd, coverage_max, coverage_mean, special_read, rg_str);
      if ( check_alu == 0 ) continue;
      if ( check_alu == COVERAGE_HIGH) {
	f_log1 << "COVERAGE_HIGH " << chrn << " " << aluBegin << " " << aluEnd << " " << setprecision(2) << coverage_mean << endl;
	continue;
      }
      map < T_READ, int > readCnt;
      count_reads(special_read, readCnt); 
      if ( readCnt[clip_read] or readCnt[unknow_read] ) // not interesting if all are mid_reads
	f_tmp1 << chrn << " " << aluBegin << " " << aluEnd << " " << setprecision(2) << coverage_mean << " " <<  readCnt[mid_read]
	       << " " << readCnt[clip_read] << " " << readCnt[unknow_read] << " " << rg_str << endl;      
    }
    delete alurefpos;
    cerr << "file_alupos:done  " << file_alupos << endl;  
  }
  f_tmp1.close();
  f_log1.close();
  seqan::close(inStream);     
  return 0;
}

bool genoProb_per_line(string &line, string & output_line, map <int, EmpiricalPdf *> & empiricalpdf_rg, float *log10_pl, float *log10_pm, float *gp) {
  stringstream ss, ss_out;
  string phred_str_l, phred_str_m,chrn, meanCov;
  int aluBegin, aluEnd, midCnt, clipCnt, unknowCnt;
  ss.clear(); ss.str(line); 
  ss >> chrn >> aluBegin >> aluEnd >> meanCov >> midCnt >> clipCnt >> unknowCnt ;
  // based on special reads
  log10_pl[0] = clipCnt * (-LOG10_RATIO_UB) ; 
  log10_pl[1] = (midCnt + clipCnt) * log10 (0.5) ;
  log10_pl[2] = midCnt * (-LOG10_RATIO_UB) ; 
  
  if ( unknowCnt == 0 )  {
    if ( p00_is_dominant(log10_pl, - LOG10_GENO_PROB) ) return false;
    ss_out << chrn << " " << aluBegin << " " << aluEnd << " " << unknowCnt ;
    //cout << "### " << log10_pl[0] << " " << log10_pl[1] << " " << log10_pl[2] << " " << endl;
    log10P_to_P(log10_pl, gp, LOG10_GENO_PROB);
    genoProb_print(gp, ss_out, 6);
    genoProb_print(gp, ss_out, 6);
  }

  if (unknowCnt) { 
    int insert_len, idx;
    string token;
    for (int i = 0; i < 3; i++)  log10_pm[i] = log10_pl[i];
    for (int i = 0; i < unknowCnt; i++) {
      getline(ss, token, ':');
      seqan::lexicalCast2(idx, token);
      getline(ss, token, ' ');
      seqan::lexicalCast2(insert_len, token);      
      float p_y = empiricalpdf_rg[idx]->pdf_obs(insert_len);
      float p_z = empiricalpdf_rg[idx]->pdf_obs(insert_len - aluEnd + aluBegin);

      //cout << "## insert_len " << insert_len << "  alu_len "  << aluEnd - aluBegin << endl;
      //cout << "### " << p_y << " " << p_z << endl;

      log10_pm[0] += log10 (p_y);
      log10_pm[1] += log10 (0.67 * p_y + 0.33 * p_z);
      log10_pm[2] += log10 (p_z); 
    }      
    if ( p00_is_dominant(log10_pl, - LOG10_GENO_PROB) and p00_is_dominant(log10_pm, - LOG10_GENO_PROB) ) 
      return false;
    ss_out << chrn << " " << aluBegin << " " << aluEnd << " " << unknowCnt ;
    log10P_to_P(log10_pl, gp, LOG10_GENO_PROB);
    genoProb_print(gp, ss_out, 6);
    log10P_to_P(log10_pm, gp, LOG10_GENO_PROB);
    genoProb_print(gp, ss_out, 6);
  }
  
  output_line = ss_out.str();
  return true;
}

void calculate_genoProb(string fn_tmp1, string fn_tmp2, map <int, EmpiricalPdf *> & empiricalpdf_rg){
  // clip_read: evidence for deletion, mid_read: evidence for insertion
  ofstream fout(fn_tmp2.c_str());
  fout << "chr aluBegin aluEnd unknowCnt l0 l1 l2 m0 m1 m2\n";  
  string line, output_line;
  ifstream fin(fn_tmp1.c_str());
  assert(fin);
  getline(fin, line); // read header
  float *log10_pl = new float[3];
  float *log10_pm = new float[3];
  float *gp = new float[3];
  while (getline(fin, line))
    if (genoProb_per_line(line, output_line, empiricalpdf_rg, log10_pl, log10_pm, gp)) 
      fout << output_line << endl;
  fin.close();
  fout.close();
  delete log10_pl;
  delete log10_pm;
  delete gp;
}
 
void filter_by_llh(string path1, string f_in_suffix, string f_out, vector <string> &pns, vector <string> &chrns, int col_00) {
  stringstream ss;
  string line, chrn, tmpfield;
  int aluBegin, flag;
  float p0, p1, p2;
  int alleleCnt = pns.size() * 2;
  map < int, int > pos_altCnt; // count of alternative alleles 
  map < int, map<string, GENO_PROB * > > pos_pnProb;  
  map < int, map<string, GENO_PROB * > >::iterator pp;
  map < string, GENO_PROB * >::iterator ppi;

  ofstream fout(f_out.c_str());  
  for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
      ifstream fin( get_name_tmp(path1, *pi, f_in_suffix).c_str() );
      getline(fin, line); 
      assert(fin);
      flag = 1;
      while ( getline(fin, line) ) {
	ss.clear(); ss.str( line );
	ss >> chrn >>  aluBegin;
	for (int ti = 0; ti < col_00 - 3; ti++) ss >> tmpfield;
	ss >> p0 >> p1 >> p2 ;
	if (chrn != *ci) {
	  if (flag) continue;
	  else break;
	}
	flag = 0;
	pos_pnProb[ aluBegin ][*pi] = new GENO_PROB( p0, p1, p2);	  
	if (p1 > p0 or p2 > p0) 
	  addKey(pos_altCnt, aluBegin, p1 > p2 ? 1 : 2);	
      }
      fin.close();
    }
    for ( map<int, int>::iterator pa = pos_altCnt.begin(); pa != pos_altCnt.end(); pa++) {
      float altFreq = (pa->second) / (float)alleleCnt;
      float freq0 = (1 - altFreq) * (1 - altFreq);
      float freq1 = 2 * altFreq * (1 - altFreq);
      float freq2 = altFreq * altFreq;
      float llh_alt = 0, llh_00 = 0;
      for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
	if ( (ppi = pos_pnProb[pa->first].find(*pi)) != pos_pnProb[pa->first].end() ) {
	  llh_00 +=  ( (ppi->second)->g0 > 0 ? log( (ppi->second)->g0 ) : ( - LOG10_GENO_PROB ) ) ;
	  llh_alt += log (freq0 * (ppi->second)->g0 + freq1 * (ppi->second)->g1 + freq2 * (ppi->second)->g2);  
	}
      }
      if (llh_alt > llh_00) 
	fout << *ci << " " << pa->first << " " << altFreq << " " << llh_alt - llh_00 << endl;
    }
    pos_altCnt.clear();
    for (pp = pos_pnProb.begin(); pp != pos_pnProb.end(); pp++) {
      for ( ppi = (pp->second).begin(); ppi != (pp->second).end(); ppi++) 
	delete ppi->second;
    }
    pos_pnProb.clear();
  } 
  fout.close();
}

void combine_pns_vcf(string path1, string f_in_suffix, string f_out, vector <string> &pns, vector <string> &chrns, int col_00) {
  ifstream fin;
  stringstream ss;
  string line, chrn, tmp1, tmp2;
  int aluBegin, aluEnd, cii, flag;

  map < pair<int, int>, map<string, string> > pos_pnProb;
  map < pair<int, int>, map<string, string> > pos_pnProb_all;
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
  float p0, p1, p2;      
  string tmpfield;
  for ( cii = 0, ci = chrns.begin(); ci != chrns.end(); ci++, cii++) {
    pos_pnProb.clear();
    pos_pnProb_all.clear();
    for ( pi = pns.begin(); pi != pns.end(); pi++) {
      fin.open( get_name_tmp(path1, *pi, f_in_suffix).c_str() );
      getline(fin, line); 
      assert(fin);
      flag = 1;
      while ( getline(fin, line) ) {
	ss.clear(); ss.str( line );
	ss >> chrn >>  aluBegin >> aluEnd;
	for (int ti = 0; ti < col_00 - 4; ti++) ss >> tmpfield;
	ss >> p0 >> p1 >> p2 ;
	if (chrn != *ci) {
	  if (flag) continue;
	  else break;
	  }
	flag = 0;
	string phred_scaled_str = phred_scaled(p0, p1, p2);
	pos_pnProb_all[ make_pair(aluBegin, aluEnd) ][*pi] = phred_scaled_str;
	if (p1 > p0 or p2 > p0) pos_pnProb[ make_pair(aluBegin, aluEnd) ][*pi] = phred_scaled_str;
      }

      fin.close();
    }
    // Write out the records.
    for (ppp = pos_pnProb.begin(); ppp != pos_pnProb.end(); ppp++) {
      record.beginPos = (ppp->first).first;
      record.rID = cii;
      record.id = int_to_string((ppp->first).second);
      for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
	if ( (pp = ppp->second.find(*pi)) != ppp->second.end() )  
	  appendValue(record.genotypeInfos, pp->second);
	else if ((pp = pos_pnProb_all[ppp->first].find(*pi)) !=  pos_pnProb_all[ppp->first].end() )
	  appendValue(record.genotypeInfos,  pos_pnProb_all[ppp->first][*pi]); 	  
	else 
	  appendValue(record.genotypeInfos, "0,255,255"); // otherwise vcf consider it as missing
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
  if (argc < 3) exit(1);

  int opt, idx_pn;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];
  string path1 = read_config(config_file, "file_alu_delete0");
  if ( access( path1.c_str(), 0 ) != 0 ) system( ("mkdir " + path1).c_str() );
  
  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  map<int, string> ID_pn;
  get_pn(read_config(config_file, "file_pn"), ID_pn);
  string file_fa_prefix = read_config(config_file, "file_fa_prefix");
  string file_dist_prefix = read_config(config_file, "file_dist_prefix");
  boost::timer clocki;    
  clocki.restart();

  if ( opt == 1 ) {  // write down counts 
    seqan::lexicalCast2(idx_pn, argv[3]);
    string chrn = argv[4];
    if (argc != 5) {
      cerr << "not enough parameters \n";
      exit(0);
    }
    string pn = ID_pn[idx_pn];
    cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
    if (chrn != "chr0") { chrns.clear();  chrns.push_back(chrn); }
    map<string, int> rg_to_idx;
    int idx = 0; 
    ifstream fin;
    string rg;
    fin.open( get_name_rg(file_dist_prefix, pn).c_str());
    assert(fin);
    while (fin >> rg) rg_to_idx[rg] = idx++;
    fin.close();    
    unsigned coverage_max;
    seqan::lexicalCast2(coverage_max, (read_config(config_file, "coverage_max")));
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";  
    string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
    string fn_tmp1 = get_name_tmp(path1, pn, ".tmp1");
    string fn_log1 = get_name_tmp(path1, pn, ".log1");
    // delete_search(bam_input, bai_input, file_fa_prefix, chrns, fn_tmp1, fn_log1, file_alupos_prefix, coverage_max, rg_to_idx);
    if ( chrn != "chr0") 
      return 0;
    // step 2: calculate prob
    string fn_tmp2 = get_name_tmp(path1, pn, ".tmp2");
    map <int, EmpiricalPdf *> empiricalpdf_rg;    
    string pdf_param = read_config(config_file, "pdf_param"); // 100_1000_5  
    fin.open( get_name_rg(file_dist_prefix, pn).c_str());
    idx = 0; 
    while (fin >> rg) 
      empiricalpdf_rg[idx++] = new EmpiricalPdf( get_name_rg_pdf(file_dist_prefix, pn, rg, pdf_param));
    fin.close();
    calculate_genoProb(fn_tmp1, fn_tmp2, empiricalpdf_rg); 
    for (map <int, EmpiricalPdf *>::iterator ri = empiricalpdf_rg.begin(); ri != empiricalpdf_rg.end(); ri++) 
      delete ri->second;
    
  } else if (opt == 2) {   // write vcf for all pn

    vector <string> pns;
    ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
    while (fin >> idx_pn) pns.push_back( ID_pn[idx_pn]);
    fin.close();
    filter_by_llh(path1, ".tmp2", path1+"filter_mr.pos", pns, chrns, 8);
    combine_pns_vcf(path1, ".tmp2", path1+"test_mr.vcf", pns, chrns, 8);
    combine_pns_vcf(path1, ".tmp2", path1+"test_lr.vcf", pns, chrns, 5);

  } else if (opt == 0) { // debugging and manually check some regions 
    string pn, chrn, bam_input, bai_input, fa_input;
    // check pdf 
    /*
    pn = "AAOOCST";
    map <int, EmpiricalPdf *> empiricalpdf_rg;    
    string rg;
    int idx = 0;
    string pdf_param = read_config(config_file, "pdf_param"); // 100_1000_5  
    ifstream fin( get_name_rg(file_dist_prefix, pn).c_str());
    while (fin >> rg) 
      empiricalpdf_rg[idx++] = new EmpiricalPdf( get_name_rg_pdf(file_dist_prefix, pn, rg, pdf_param));
    fin.close();
    string line = "chr16 19137943 19138261 11 1 0 3 1:700 1:603 1:675";
    string output_line;
    float *log10_pl = new float[3];
    float *log10_pm = new float[3];
    float *gp = new float[3];
    genoProb_per_line(line, output_line, empiricalpdf_rg, log10_pl, log10_pm, gp);
    cout << output_line << endl;
    delete log10_pl;
    delete log10_pm;
    delete gp;
    return 0;
    */

    // check some regions 
    // genotype call /nfs_mount/bioinfo/users/yuq/work/Alu/outputs/jon_chr0/pn1.check
    pn = "AAFDUTV";
    chrn = "chr3";
    int pa = 160911683;
    int pb = 160911872;
    bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    bai_input = bam_input + ".bai";  
    fa_input = file_fa_prefix + chrn + ".fa";    
    check_delete_region(bam_input, bai_input, fa_input, chrn, pa - 600,  pb + 600);
    
  }

  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
