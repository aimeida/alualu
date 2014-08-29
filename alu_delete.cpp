#define SEQAN_HAS_ZLIB 1
#include "delete_utils.h"

void count_reads(map <seqan::CharString, T_READ> &qName_info, map < T_READ, int > &readCnt) {
  readCnt.clear();
  for (map <seqan::CharString, T_READ>::iterator rt = qName_info.begin(); rt != qName_info.end(); rt++) 
    addKey(readCnt, rt->second); 
}

int check_one_pos(BamFileHandler* bam_fh, FastaFileHandler *fasta_fh, map <string, int> &rg_to_idx, string chrn, int aluBegin, int aluEnd, unsigned coverage_max, int &totalCnt, map <seqan::CharString, T_READ> &qName_info,  map<seqan::CharString, string> & rg_str){  
  qName_info.clear();
  int reads_cnt = 0;
  seqan::BamAlignmentRecord record;  
  while ( true ) {
    string read_status = bam_fh->fetch_a_read(chrn, aluBegin - ALU_FLANK, aluEnd + ALU_FLANK, record);
    if (read_status == "stop" ) break;
    if (read_status == "skip" or !QC_delete_read(record)) continue;
    reads_cnt ++;  
    map <seqan::CharString, T_READ >::iterator qItr = qName_info.find(record.qName);
    if ( qItr != qName_info.end() and qItr->second != unknow_read) 
      continue;
    T_READ iread = classify_read( record, aluBegin, aluEnd, fasta_fh);   
    if ( iread == useless_read) continue;    
    qName_info[record.qName] = iread;
    if ( qItr == qName_info.end() and iread == unknow_read) {
      int rgIdx = get_rgIdx(rg_to_idx, record);
      stringstream rg_ss;
      rg_ss << rgIdx << ":" << abs(record.tLen);
      rg_str[record.qName] = rg_ss.str();
    }    
  }
  totalCnt = reads_cnt;
  if ( length(record.seq) * (float) reads_cnt / (aluEnd - aluBegin + 2 * ALU_FLANK) > coverage_max) return COVERAGE_HIGH;
  return 1;
}

int delete_search(int minLen_alu_del, BamFileHandler *bam_fh, string file_fa_prefix, vector<string> &chrns, string &f_out, string &f_log, string &file_alupos_chr0, int coverage_max, map<string, int> &rg_to_idx) {    
  map < seqan::CharString, T_READ> qName_info;  
  ofstream f_tmp1( f_out.c_str()); 
  f_tmp1 << "chr aluBegin aluEnd totalCnt midCnt clipCnt unknowCnt unknowStr\n";
  ofstream f_log1( f_log.c_str());  // p    rint out info for clip reads 
  for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++) {
    string chrn = *ci;
    string file_alupos = replace_str0_str(file_alupos_chr0, chrn, "chr0");
    AluRefPos *alurefpos = new AluRefPos(file_alupos, minLen_alu_del, 300);  // min distance to neighbor is 300 bp 
    FastaFileHandler *fasta_fh = new FastaFileHandler(file_fa_prefix + chrn + ".fa", chrn);    
    int aluBegin, aluEnd;
    for (int count_loci = 0; ; count_loci++) {
      int totalCnt = 0;
      if ( !alurefpos->nextdb() ) break;      
      aluBegin = alurefpos->get_beginP();
      aluEnd = alurefpos->get_endP();
      if (aluBegin <= ALU_FLANK or !bam_fh->jump_to_region(chrn, aluBegin-ALU_FLANK, aluEnd + ALU_FLANK))
	continue;
      map<seqan::CharString, string> rg_str;
      int check_alu = check_one_pos(bam_fh, fasta_fh, rg_to_idx, chrn, aluBegin, aluEnd, coverage_max, totalCnt, qName_info, rg_str);
      if ( check_alu == 0 ) continue;
      if ( check_alu == COVERAGE_HIGH) {
	f_log1 << "COVERAGE_HIGH " << chrn << " " << aluBegin << " " << aluEnd << " " << totalCnt << endl;
	continue;
      }
      map < T_READ, int > readCnt;
      count_reads(qName_info, readCnt); 
      if ( readCnt[clip_read] or readCnt[unknow_read] ) {// not interesting if all are mid_reads
	f_tmp1 << chrn << " " << aluBegin << " " << aluEnd << " " << totalCnt << " " <<  readCnt[mid_read]
	       << " " << readCnt[clip_read] << " " << readCnt[unknow_read];
	for (map<seqan::CharString, string>::iterator ri = rg_str.begin(); ri != rg_str.end(); ri++) 
	  if (qName_info[ri->first] == unknow_read)
	    f_tmp1 << " " << ri->second ;
	f_tmp1 << endl;
      } else if (totalCnt > 0)
	f_tmp1 << chrn << " " << aluBegin << " " << aluEnd << " " <<  -totalCnt << endl;  // different as missing data 
    }
    delete alurefpos;
    delete fasta_fh;
    cout << "file_alupos:done  " << file_alupos << endl;  
  }
  f_tmp1.close();
  f_log1.close();
  return 0;
}

int parseline_del(ostream & fout, string &line, map <int, EmpiricalPdf *> & pdf_rg, float logPE){
  const int read_len = 100;
  float *log10_gp = new float[3];
  stringstream ss;
  string chrn;
  int aluBegin, aluEnd, totalCnt, midCnt, clipCnt, unknowCnt;
  ss.clear(); ss.str(line); 
  ss >> chrn >> aluBegin >> aluEnd >> totalCnt;
  if (totalCnt < 0) {    // G00 > 1e-5, not missing data
    fout << chrn << " " << aluBegin << " " << aluEnd << " " << totalCnt << endl;  
    return 0;
  } 
  // G00 < 1e-5
  ss >> midCnt >> clipCnt >> unknowCnt ;
  float bias1 = (aluEnd - aluBegin + 3. * read_len) / (3. * read_len); 
  float ph0 = bias1 / (bias1 + 1);    
  logPE = - abs(logPE);
  float logPM = log10 ( 1 - pow(10, logPE) );
  
  float *gp = new float[3];
  for (int i = 0; i < 3; i++) log10_gp[i] = 0;
  if (midCnt + clipCnt > 0) {
    log10_gp[0] = clipCnt * logPE + midCnt * logPM;
    log10_gp[1] = midCnt * log10 (ph0) + clipCnt * log10 (1 - ph0) + (midCnt + clipCnt) * logPM  ; 
    log10_gp[2] = midCnt * logPE + clipCnt * logPM;
  }
  //cout << "debug1# " << log10_gp[0] << " " << log10_gp[1] << " " << log10_gp[2] << endl;   
  if (unknowCnt) { 
    int insert_len, idx;
    string token;
    for (int i = 0; i < unknowCnt; i++) {
      getline(ss, token, ':');
      seqan::lexicalCast2(idx, token);
      getline(ss, token, ' ');
      seqan::lexicalCast2(insert_len, token);      
      float p_y, p_z;
      pdf_rg[idx]->ratio_obs(insert_len, insert_len - aluEnd + aluBegin, abs(logPE) - 0.5, p_y, p_z);
      log10_gp[0] += log10 (p_y);
      log10_gp[1] += log10 (ph0 * p_y + (1 - ph0) * p_z) ;
      log10_gp[2] += log10 (p_z);
    }
  }  
  //cout << "debug1# " << log10_gp[0] << " " << log10_gp[1] << " " << log10_gp[2] << endl;   
  if ( !p00_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    fout << chrn << " " << aluBegin << " " << aluEnd << " " << midCnt << " " << clipCnt << " " << unknowCnt ;
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
    fout << " " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2] << endl;   
  } else {
    fout << chrn << " " << aluBegin << " " << aluEnd << " " << -totalCnt << endl;  
  }
  delete gp;    
  delete log10_gp;
  return 1;
}

void read_highCov_region(string f_input, map < string, set<int> > & chrn_aluBegin) {
  ifstream fin(f_input.c_str());
  string tmp1, tmp2, chrn;
  int aluBegin, aluEnd;
  while ( fin >> tmp1>> chrn >> aluBegin >> aluEnd >> tmp2 ) chrn_aluBegin[chrn].insert(aluBegin);
  fin.close();
}

void write_rm1(string f_output, map < string, std::set<int> > & chrn_aluBegin) {
  ofstream fout(f_output.c_str());
  for (map < string, std::set<int> >::iterator it = chrn_aluBegin.begin(); it != chrn_aluBegin.end(); it++) 
    for ( std::set<int>::iterator i2 = (it->second).begin(); i2 != (it->second).end(); i2++ )
      fout << it->first << " " << *i2 << " highCov\n";
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
  string s_chrns = cf_fh.get_conf("chrns");
  parse_chrns(s_chrns, chrns);

  string file_fa_prefix = cf_fh.get_conf("file_fa_prefix");
  string file_dist_prefix = cf_fh.get_conf("file_dist_prefix");

  string path0 = cf_fh.get_conf("file_alu_delete0");
  check_folder_exists(path0);
  
  float log10RatioUb = seqan::lexicalCast<float> (cf_fh.get_conf("LOG10_RATIO_UB"));
  string file_alupos_prefix = cf_fh.get_conf("file_alupos_prefix"); 
  check_folder_exists(file_alupos_prefix);

  if ( opt == "preprocess" ) {  
    // sort alu files 
    for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++)   // combine Alu if less than 10 
      sort_file_by_col<int> (file_alupos_prefix + "alu_" + *ci, 2, false);
    
  } else if ( opt == "write_tmps_pn" ) { 
    int idx_pn = seqan::lexicalCast<int> (argv[3]);
    assert(argc == 4);
    string pn = ID_pn[idx_pn];
    map<string, int> rg_to_idx;
    parse_reading_group( get_name_rg(file_dist_prefix, pn), rg_to_idx );

    string fn_tmp1 = path0 + pn + ".tmp1";
    string fn_log1 = path0 + pn + ".log1";
    string bam_input = cf_fh.get_conf("file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";      
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bai_input);
    string file_alupos = file_alupos_prefix + "alu_chr0";

    int minLen_alu_del = seqan::lexicalCast <int> (cf_fh.get_conf("minLen_alu_del"));
    unsigned coverage_max = seqan::lexicalCast<unsigned> (cf_fh.get_conf("coverage_max"));
    delete_search(minLen_alu_del, bam_fh, file_fa_prefix, chrns, fn_tmp1, fn_log1, file_alupos, coverage_max, rg_to_idx);
    delete bam_fh;
    move_files(path0+"log1s/", path0 + pn + ".log1") ;

    string fn_tmp2 = path0 + pn + ".tmp2";
    map <int, EmpiricalPdf *> pdf_rg;    
    string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
    read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
    
    ofstream fout(fn_tmp2.c_str());
    fout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";  
    string line;
    ifstream fin(fn_tmp1.c_str());
    assert(fin);
    getline(fin, line); // read header
    while (getline(fin, line))
      parseline_del(fout, line, pdf_rg, log10RatioUb);
    fin.close();
    fout.close();
    EmpiricalPdf::delete_map(pdf_rg);
    move_files(path0+"tmp1s/", path0 + pn + ".tmp1") ;
    move_files(path0+"tmp2s/", path0 + pn + ".tmp2") ;
    
  } else if (opt == "write_vcf_pns") {   // write vcf for all pn

    vector <string> pns;
    read_file_pn_used(cf_fh.get_conf("pn_del_vcf"), pns);  // select some pns for writing vcf files 
    map <string, std::set<int> > chrn_aluBegin;
    for ( vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) 
      read_highCov_region(path0 + "log1s/" + *pi + ".log1",chrn_aluBegin);

    string path_input = path0 + "tmp2s/";
    string fn_prefix = path0 + int_to_string( pns.size()) ;
    string ref_name = cf_fh.get_conf("ref_name");
    float llh_th = 1.92; //qchisq(0.95, df=1)  [1] 3.841459 
    combine_pns_vcf(path_input, ".tmp2", fn_prefix + ".vcf", pns, chrns, chrn_aluBegin, llh_th, ref_name);  

  } else if (opt == "debug1") { // manually check some regions 
    
    string pn = argv[3];
    string chrn = "chr21";
    int aluBegin = 34850160;
    int aluEnd = 34850464;
    map < seqan::CharString, T_READ> qName_info;
    FastaFileHandler *fasta_fh = new FastaFileHandler(file_fa_prefix + chrn + ".fa", chrn);
    string bam_input = cf_fh.get_conf("file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bai_input, "debug.bam");
    seqan::BamAlignmentRecord record;
    bam_fh->jump_to_region(chrn, aluBegin-ALU_FLANK, aluEnd + ALU_FLANK);
    while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, aluBegin - ALU_FLANK, aluEnd + ALU_FLANK, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_delete_read(record)) continue;
      map <seqan::CharString, T_READ >::iterator qItr = qName_info.find(record.qName);
      if ( qItr != qName_info.end() and qItr->second != unknow_read)
        continue;
      T_READ iread = classify_read( record, aluBegin, aluEnd, fasta_fh);
      // check a specific read
      if (string(seqan::toCString(record.qName)) == string("hap:26|sim:2019607")) {
	iread = classify_read( record, aluBegin, aluEnd, fasta_fh, true);
	bam_fh->write_a_read(record);  // write sam file 
	cout << tread_toStr(iread) << " ";
	debug_print_read(record);
      }

      if ( iread == useless_read) continue;
      qName_info[record.qName] = iread;
    }
    bam_fh->jump_to_region(chrn, aluBegin-ALU_FLANK, aluEnd + ALU_FLANK);
    while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, aluBegin - ALU_FLANK, aluEnd + ALU_FLANK, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_delete_read(record)) continue;
      T_READ iread = useless_read;
      get_mapVal(qName_info, record.qName, iread);
      if ( iread == unknow_read)  {
	cout << aluBegin << " " << aluEnd << " ";
	debug_print_read(record);
      }
    }
    delete bam_fh;
     
  } else if (opt == "debug2") { // debugging and manually check some regions 

    string pn = argv[3];
    string line, output_line;
    map <int, EmpiricalPdf *> pdf_rg;    
    string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
    read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
    cout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";    
    //line = "chr21 21903875 21904163 26 48 5 2 0:525 0:534";
    //line = "chr21 41931103 41931403 30 62 6 2 0:514 0:550";
    //line = "chr21 33299674 33299968 22 56 0 7 0:496 0:805 0:826 0:461 0:809 0:457 0:834";
    line = "chr21 21903875 21904163 26 10 1 0";
    parseline_del(cout, line, pdf_rg, log10RatioUb);
    line = "chr21 21903875 21904163 26 1 10 0";
    parseline_del(cout, line, pdf_rg, log10RatioUb);

    line = "chr21 21903875 21904163 26 0 6 0";
    parseline_del(cout, line, pdf_rg, log10RatioUb);

    line = "chr21 21903875 21904163 26 6 0 0";

  } else {
    cout << "unknown options !\n";
  }
  
  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
