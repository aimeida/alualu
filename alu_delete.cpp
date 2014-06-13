#define SEQAN_HAS_ZLIB 1
#include <seqan/vcf_io.h>
#include "delete_utils.h"

inline string get_name_tmp(string path, string fn, string suffix){ return path + fn + suffix;}

void count_reads(map <seqan::CharString, T_READ> &qName_info, map < T_READ, int > &readCnt) {
  readCnt.clear();
  for (map <seqan::CharString, T_READ>::iterator rt = qName_info.begin(); rt != qName_info.end(); rt++) 
    addKey(readCnt, rt->second); 
}
 
int check_chr_alupos(BamFileHandler* bam_fh, FastaFileHandler *fasta_fh, map <string, int> &rg_to_idx, string chrn, int aluBegin, int aluEnd, unsigned coverage_max, float &coverage_mean, map <seqan::CharString, T_READ> &qName_info,  map<seqan::CharString, string> & rg_str){  
  qName_info.clear();
  int reads_cnt = 0;
  seqan::BamAlignmentRecord record;
  
  while ( true ) {
    string read_status = bam_fh->fetch_a_read(chrn, aluBegin - ALU_FLANK, aluEnd + ALU_FLANK, record);
    if (read_status == "stop" ) break;
    if (read_status == "skip" or !QC_delete_read(record)) continue;
    reads_cnt ++;  
    int align_len = getAlignmentLengthInRef(record);    
    map <seqan::CharString, T_READ >::iterator qItr = qName_info.find(record.qName);
    if ( qItr != qName_info.end() and qItr->second != unknow_read) continue;
    T_READ iread = classify_read( record, align_len, aluBegin, aluEnd, fasta_fh);
    if ( iread == useless_read) continue;
    qName_info[record.qName] = iread;
    if ( qItr == qName_info.end() and iread == unknow_read) {
      int rgIdx = get_rgIdx(rg_to_idx, record);
      stringstream rg_ss;
      rg_ss << rgIdx << ":" << abs(record.tLen);
      rg_str[record.qName] = rg_ss.str();
    }    
  }
  coverage_mean = length(record.seq) * (float) reads_cnt / (aluEnd - aluBegin + 2 * ALU_FLANK) ;
  if ( coverage_mean > coverage_max) return COVERAGE_HIGH;
  return 1;
}

int delete_search(int minLen_alu_del, string & bam_input, string &bai_input, string file_fa_prefix, vector<string> &chrns, string &f_out, string &f_log, string &file_alupos_prefix, int coverage_max, map<string, int> &rg_to_idx) {    
  map < seqan::CharString, T_READ> qName_info;  
  ofstream f_tmp1( f_out.c_str()); 
  f_tmp1 << "chr aluBegin aluEnd mean_coverage midCnt clipCnt unknowCnt unknowStr\n";
  ofstream f_log1( f_log.c_str());  // print out info for clip reads 

  BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bai_input);
  for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++) {
    string chrn = *ci;
    string file_alupos = file_alupos_prefix + chrn;
    AluRefPosRead *alurefpos = new AluRefPosRead(file_alupos, minLen_alu_del); // default 200
    FastaFileHandler *fasta_fh = new FastaFileHandler(file_fa_prefix + chrn + ".fa", chrn);
    
    int aluBegin, aluEnd;
    for (int count_loci = 0; ; count_loci++) {
      float coverage_mean = 0;
      if (! alurefpos->updatePos(aluBegin, aluEnd)) 
	break;      
      if (aluBegin <= ALU_FLANK or !bam_fh->jump_to_region(chrn, aluBegin-ALU_FLANK, aluEnd + ALU_FLANK))
	continue;
      map<seqan::CharString, string> rg_str;
      int check_alu = check_chr_alupos(bam_fh, fasta_fh, rg_to_idx, chrn, aluBegin, aluEnd, coverage_max, coverage_mean, qName_info, rg_str);
      if ( check_alu == 0 ) continue;
      if ( check_alu == COVERAGE_HIGH) {
	f_log1 << "COVERAGE_HIGH " << chrn << " " << aluBegin << " " << aluEnd << " " << setprecision(2) << coverage_mean << endl;
	continue;
      }
      map < T_READ, int > readCnt;
      count_reads(qName_info, readCnt); 
      if ( readCnt[clip_read] or readCnt[unknow_read] ) {// not interesting if all are mid_reads
	f_tmp1 << chrn << " " << aluBegin << " " << aluEnd << " " << setprecision(2) << coverage_mean << " " <<  readCnt[mid_read]
	       << " " << readCnt[clip_read] << " " << readCnt[unknow_read];
	for (map<seqan::CharString, string>::iterator ri = rg_str.begin(); ri != rg_str.end(); ri++) 
	  if (qName_info[ri->first] == unknow_read)
	    f_tmp1 << " " << ri->second ;
	f_tmp1 << endl;
      }
    }
    delete alurefpos;
    delete fasta_fh;
    cerr << "file_alupos:done  " << file_alupos << endl;  
  }
  f_tmp1.close();
  f_log1.close();
  delete bam_fh;
  return 0;
}

bool parseline_del_tmp1(string &line, string & output_line, map <int, EmpiricalPdf *> & pdf_rg){
  float *log10_gp = new float[3];
  stringstream ss, ss_out;
  string chrn, meanCov;
  int aluBegin, aluEnd, midCnt, clipCnt, unknowCnt;
  ss.clear(); ss.str(line); 
  ss >> chrn >> aluBegin >> aluEnd >> meanCov >> midCnt >> clipCnt >> unknowCnt ;
  float prob_ub = pow(10, -LOG10_RATIO_UB);
  float prob_known = (midCnt+clipCnt)/(float)(midCnt + clipCnt + unknowCnt);

  for (int i = 0; i < 3; i++) log10_gp[i] = 0;
  if (prob_known) {
    log10_gp[0] = clipCnt * log10 ( prob_known * prob_ub ) + midCnt * log10 ( prob_known * (1 - prob_ub) );
    log10_gp[1] = (midCnt + clipCnt) * log10 (prob_known * 0.5) ; 
    log10_gp[2] = midCnt * log10 ( prob_known * prob_ub ) + clipCnt * log10 ( prob_known * (1 - prob_ub) );
  }
  if (unknowCnt) { 
    int insert_len, idx;
    string token;
    for (int i = 0; i < unknowCnt; i++) {
      getline(ss, token, ':');
      seqan::lexicalCast2(idx, token);
      getline(ss, token, ' ');
      seqan::lexicalCast2(insert_len, token);      
      float p_y = pdf_rg[idx]->pdf_obs(insert_len);
      float p_z = pdf_rg[idx]->pdf_obs(insert_len - aluEnd + aluBegin);
      //float freq0 = 0.67;  // high FP ratio      
      float freq0 = ( midCnt + 1 )/(float)(midCnt + clipCnt + 2); // 1 and 2 are psudo count
      log10_gp[0] += log10 (p_y * (1 - prob_known));
      log10_gp[1] += log10 ((freq0 * p_y + (1 - freq0) * p_z) * (1 - prob_known) ) ;
      log10_gp[2] += log10 (p_z * (1 - prob_known));
    }
  }
  bool use_this_line = false;
  if ( !p00_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    ss_out << chrn << " " << aluBegin << " " << aluEnd << " " << midCnt << " " << clipCnt << " " << unknowCnt ;
    float *gp = new float[3];
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
    ss_out << " " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2]; 
    delete gp;    
    output_line = ss_out.str();
    use_this_line = true;
  }
  delete log10_gp;
  return use_this_line;
}

void calculate_genoProb(string fn_tmp1, string fn_tmp2, map <int, EmpiricalPdf *> & pdf_rg){
  // clip_read: evidence for deletion, mid_read: evidence for insertion
  ofstream fout(fn_tmp2.c_str());
  fout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";  
  string line, output_line;
  ifstream fin(fn_tmp1.c_str());
  assert(fin);
  getline(fin, line); // read header
  while (getline(fin, line))
    if (parseline_del_tmp1(line, output_line, pdf_rg))
      fout << output_line << endl;
  fin.close();
  fout.close();
}

void filter_by_llh_noPrivate(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> &chrns, int col_00) {
  stringstream ss;
  string line, chrn, tmpfield;
  int aluBegin, flag;
  float p0, p1, p2;

  map < int, int > pos_altCnt; // count of alternative alleles 
  map < int, int > pos_pnCnt;  // count of pns having polymorphism at this loci
  map < int, map<string, GENO_PROB * > > pos_pnProb;  
  map < int, map<string, GENO_PROB * > >::iterator pp;
  map < string, GENO_PROB * >::iterator ppi;

  ofstream fout(f_out.c_str());  
  fout <<  "chrn aluBegin pnCnt alleleFreq Log_LikelihoodRatio\n"; 
  int alleleCnt = pns.size() * 2; 
  for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
      ifstream fin( get_name_tmp(path0, *pi, f_in_suffix).c_str() );
      //cout << "reading " << get_name_tmp(path0, *pi, f_in_suffix) << endl;
      assert(fin);
      getline(fin, line); 
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
	if (p1 > p0 or p2 > p0) {
	  addKey(pos_altCnt, aluBegin, p1 > p2 ? 1 : 2);	
	  addKey(pos_pnCnt, aluBegin, 1);
	}
      }
      fin.close();
    }
    
    map<int, int>::iterator pan = pos_pnCnt.begin();
    for ( map<int, int>::iterator pa = pos_altCnt.begin(); pa != pos_altCnt.end(); pa++, pan++) {
      if ( pan->second <= 1) continue; // no private 
      float altFreq = (pa->second) / (float)alleleCnt;
      float freq0 = (1 - altFreq) * (1 - altFreq);
      float freq1 = 2 * altFreq * (1 - altFreq);
      float freq2 = altFreq * altFreq;
      float llh_alt = 0, llh_00 = 0;
      for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
	if ( (ppi = pos_pnProb[pa->first].find(*pi)) != pos_pnProb[pa->first].end() ) {
	  llh_00 +=  ( (ppi->second)->g0 > 0 ? log( (ppi->second)->g0 ) : ( - LOG10_GENO_PROB ) ) ;
	  llh_alt += log (freq0 * (ppi->second)->g0 + freq1 * (ppi->second)->g1 + freq2 * (ppi->second)->g2);  
	} else
	  llh_alt += log (freq0);  // g0 = 1	
      }
      fout << *ci << " " << pa->first << " " << pan->second << " " << altFreq << " " << llh_alt - llh_00 << endl;
    }
    cout << "done with " << *ci << endl;

    pos_altCnt.clear(); 
    pos_pnCnt.clear();
    for (pp = pos_pnProb.begin(); pp != pos_pnProb.end(); pp++) {
      for ( ppi = (pp->second).begin(); ppi != (pp->second).end(); ppi++) 
	delete ppi->second;
    }
    pos_pnProb.clear();
  } 
  fout.close();
}

bool combine_pns_vcf_noPrivate(string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> & chrns, int col_00) {
  ifstream fin;
  stringstream ss;
  string line, chrn, tmp1, tmp2;
  int aluBegin, aluEnd, cii, flag;

  map < pair<int, int>, map<string, string> > pos_pnProb;
  map < pair<int, int>, map<string, string> > pos_pnProb_all;
  map < pair<int, int>, map<string, string> >::iterator ppp;
  map<string, string>::iterator pp;
  vector <string>::iterator ci, pi;
  //seqan::VcfStream vcfout("-", seqan::VcfStream::WRITE);
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
      fin.open( get_name_tmp(path0, *pi, f_in_suffix).c_str() );
      assert(fin);
      getline(fin, line); 
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
      //cout << get_name_tmp(path0, *pi, f_in_suffix) << endl;
      fin.close();
    }
    cout << *ci << " size " << pos_pnProb.size() << endl;
    // Write out the records.
    for (ppp = pos_pnProb.begin(); ppp != pos_pnProb.end(); ppp++) {
      record.beginPos = (ppp->first).first;
      record.rID = cii;
      record.id = int_to_string((ppp->first).second);
      int n_pn = 0;
      for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
	if ( (pp = ppp->second.find(*pi)) != ppp->second.end() ) { 
	  appendValue(record.genotypeInfos, pp->second);
	  n_pn++ ;
	} else if ((pp = pos_pnProb_all[ppp->first].find(*pi)) !=  pos_pnProb_all[ppp->first].end() ) {
	  appendValue(record.genotypeInfos,  pos_pnProb_all[ppp->first][*pi]); 	  
	} else {
	  appendValue(record.genotypeInfos, "0,255,255"); // otherwise vcf consider it as missing
	}
      }
      if (n_pn > 1) writeRecord(vcfout, record);  // don't output private loci
      clear(record.genotypeInfos);
    }        
    cout << "done with " << *ci << endl;
  }  // chrn finished
  clear(record);
  seqan::close(vcfout);
  return true;
}

void read_highCov_region(string f_input, map < string, set<int> > & chrn_aluBegin) {
  ifstream fin(f_input.c_str());
  string tmp1, tmp2, chrn;
  int aluBegin, aluEnd;
  while ( fin >> tmp1>> chrn >> aluBegin >> aluEnd >> tmp2 ) chrn_aluBegin[chrn].insert(aluBegin);
  fin.close();
}

void remove_highCov_region(string f_input, string f_output, int offset, map <string, set<int> > &chrn_aluBegin) {
  string line, chrn;
  int aluBegin;
  stringstream ss;
  ifstream fin(f_input.c_str());
  assert(fin);
  ofstream fout(f_output.c_str());
  getline(fin, line);
  fout << line << endl;
  while (getline(fin, line)) {
    if (line[0]=='#') {
      fout << line << endl;
      continue;
    }
    ss.clear(); ss.str(line); 
    ss >> chrn >> aluBegin;
    aluBegin = aluBegin + offset;
    if (chrn_aluBegin[chrn].find(aluBegin) == chrn_aluBegin[chrn].end())
      fout << line << endl;
  }
  fin.close();
  fout.close();
}

int main( int argc, char* argv[] )
{
  string config_file = argv[1];
  string opt = argv[2];
  if (argc < 3) exit(1);

  boost::timer clocki;    
  clocki.restart();  
  ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
  map<int, string> ID_pn;
  get_pn(cf_fh.get_conf("file_pn"), ID_pn);

  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  string file_fa_prefix = cf_fh.get_conf("file_fa_prefix");
  string file_dist_prefix = cf_fh.get_conf("file_dist_prefix");

  if ( opt == "write_tmp1" ) { 
    int idx_pn = seqan::lexicalCast<int> (argv[3]);
    assert(argc == 4);
    string pn = ID_pn[idx_pn];
    cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";

    string path0 = cf_fh.get_conf("file_alu_delete0");
    check_folder_exists(path0);
    map<string, int> rg_to_idx;
    parse_reading_group( get_name_rg(file_dist_prefix, pn), rg_to_idx );
    
    unsigned coverage_max = seqan::lexicalCast<unsigned> (cf_fh.get_conf("coverage_max"));
    string bam_input = cf_fh.get_conf("file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";  
    string file_alupos_prefix = cf_fh.get_conf("file_alupos_prefix"); 
    string fn_tmp1 = get_name_tmp(path0, pn, ".tmp1");
    string fn_log1 = get_name_tmp(path0, pn, ".log1");
    int minLen_alu_del = seqan::lexicalCast <int> (cf_fh.get_conf("minLen_alu_del"));
    delete_search(minLen_alu_del, bam_input, bai_input, file_fa_prefix, chrns, fn_tmp1, fn_log1, file_alupos_prefix, coverage_max, rg_to_idx);
    string path_move = path0 + "log1s/";
    check_folder_exists(path_move);
    system(("mv " + path0 + pn + ".log1 " + path_move).c_str());
    path_move = path0 + "tmp1s/";
    check_folder_exists(path_move);
    system(("mv " + path0 + pn + ".tmp1 " + path_move).c_str());

  } else if ( opt == "write_tmp2" ) {
    int idx_pn = seqan::lexicalCast<int> (argv[3]);
    assert(argc != 4);
    string pn = ID_pn[idx_pn];
    string path0 = cf_fh.get_conf("file_alu_delete0");    
    string fn_tmp1 = get_name_tmp(path0 + "tmp1s/", pn, ".tmp1");
    string fn_tmp2 = get_name_tmp(path0, pn, ".tmp2");
    map <int, EmpiricalPdf *> pdf_rg;    
    string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
    read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
    calculate_genoProb(fn_tmp1, fn_tmp2, pdf_rg); 
    EmpiricalPdf::delete_map(pdf_rg);

    string path_move = path0 + "tmp2s/";
    check_folder_exists(path_move);
    system(("mv " + path0 + pn + ".tmp2 " + path_move).c_str());
    
  } else if (opt == "write_vcf") {   // write vcf for all pn
    string path0 = cf_fh.get_conf("file_alu_delete0");
    string path1 = cf_fh.get_conf("file_alu_delete1");
    check_folder_exists(path1);
    string pn;
    vector <string> pns;
    int i = 0, ni = 3000;
    ifstream fin(cf_fh.get_conf("file_pn").c_str());
    while ( i++ < ni and  fin >> pn) pns.push_back( pn );
    fin.close();    
    string path_input = path0 + "tmp2s/";
    string fn_pos, fn_vcf;
    fn_pos = path1 + int_to_string( pns.size()) + ".pos";
    //filter_by_llh_noPrivate(path_input, ".tmp2", fn_pos + ".tmp", pns, chrns, 7);
    fn_vcf = path1 + int_to_string( pns.size()) + ".vcf";  
    //combine_pns_vcf_noPrivate(path_input, ".tmp2", fn_vcf + ".tmp", pns, chrns, 7);  //  10 mins
    
    map <string, set<int> > chrn_aluBegin;
    for ( vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) {
      string fn_log1 = get_name_tmp(path0 + "log1s/", *pi, ".log1");  
      read_highCov_region(fn_log1, chrn_aluBegin);
    }
    remove_highCov_region(fn_pos+".tmp", fn_pos, 0, chrn_aluBegin);
    remove_highCov_region(fn_vcf+".tmp", fn_vcf, -1, chrn_aluBegin);

  } else if (opt == "debug") { // debugging and manually check some regions 
    string pn = ID_pn[0];

    /* print header info
    string bam_input = cf_fh.get_conf("file_bam_prefix") + pn + ".bam";
    seqan::BamStream bamIO(bam_input.c_str());
    for (unsigned i = 0; i < length(bamIO.header.records); ++i) {
      cout << i << " " << bamIO.header.records[i].tags[0].i1 
	   << " " << bamIO.header.records[i].tags[0].i2 
	   << " " << bamIO.header.records[i].tags[1].i1 
	   << " " << bamIO.header.records[i].tags[1].i2 << endl;
    }
    */

    float *log10_gp = new float[3];
    float *gp = new float[3];
    float offset = 0.4;
    log10_gp[0] = -4.168 + offset;
    log10_gp[1] = -0.60814 + offset;
    log10_gp[2] = -23.9364 + offset;
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);
    cout << "test1 done " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2]; 
     
    string line, output_line;
    map <int, EmpiricalPdf *> pdf_rg;    
    string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
    
    cout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";
    read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
    line = "chr1 1455790 1456098 6 2 1 0";
    parseline_del_tmp1(line, output_line, pdf_rg);
    cout << output_line << endl;
    line = "chr1 1455790 1456098 6 2 0 2 0:733 1:748";
    parseline_del_tmp1(line, output_line, pdf_rg);
    cout << output_line << endl;
    line = "chr1 1455790 1456098 6 9 1 2 0:733 1:748";
    parseline_del_tmp1(line, output_line, pdf_rg);
    cout << output_line << endl;
    cout << "test2 done\n";
    delete log10_gp;
    delete gp;
    // genotype call /nfs_mount/bioinfo/users/yuq/work/Alu/outputs/jon_chr0/pn1.check
  } else {
    cout << "unknown options !\n";
  }
  
  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
