#define SEQAN_HAS_ZLIB 1
#include "delete_utils.h"

void count_reads(map <seqan::CharString, T_READ> &qName_info, map < T_READ, int > &readCnt) {
  readCnt.clear();
  for (map <seqan::CharString, T_READ>::iterator rt = qName_info.begin(); rt != qName_info.end(); rt++) 
    addKey(readCnt, rt->second); 
}

int check_one_pos(BamFileHandler* bam_fh, FastaFileHandler *fasta_fh, map <string, int> &rg_to_idx, string chrn, int aluBegin, int aluEnd, unsigned coverage_max, float &coverage_mean, map <seqan::CharString, T_READ> &qName_info,  map<seqan::CharString, string> & rg_str){  
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
  coverage_mean = length(record.seq) * (float) reads_cnt / (aluEnd - aluBegin + 2 * ALU_FLANK) ;
  if ( coverage_mean > coverage_max) return COVERAGE_HIGH;
  return 1;
}

int delete_search(int minLen_alu_del, BamFileHandler *bam_fh, string file_fa_prefix, vector<string> &chrns, string &f_out, string &f_log, string &file_alupos_prefix, int coverage_max, map<string, int> &rg_to_idx) {    
  map < seqan::CharString, T_READ> qName_info;  
  ofstream f_tmp1( f_out.c_str()); 
  f_tmp1 << "chr aluBegin aluEnd mean_coverage midCnt clipCnt unknowCnt unknowStr\n";
  ofstream f_log1( f_log.c_str());  // print out info for clip reads 
  for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++) {
    string chrn = *ci;
    string file_alupos = file_alupos_prefix + chrn;
    AluRefPos *alurefpos = new AluRefPos(file_alupos, minLen_alu_del); // default 200
    FastaFileHandler *fasta_fh = new FastaFileHandler(file_fa_prefix + chrn + ".fa", chrn);    
    int aluBegin, aluEnd;
    for (int count_loci = 0; ; count_loci++) {
      float coverage_mean = 0;
      if ( !alurefpos->nextdb() ) break;      
      aluBegin = alurefpos->get_beginP();
      aluEnd = alurefpos->get_endP();
      if (aluBegin <= ALU_FLANK or !bam_fh->jump_to_region(chrn, aluBegin-ALU_FLANK, aluEnd + ALU_FLANK))
	continue;
      map<seqan::CharString, string> rg_str;
      int check_alu = check_one_pos(bam_fh, fasta_fh, rg_to_idx, chrn, aluBegin, aluEnd, coverage_max, coverage_mean, qName_info, rg_str);
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
    cout << "file_alupos:done  " << file_alupos << endl;  
  }
  f_tmp1.close();
  f_log1.close();
  return 0;
}

bool parseline_del_tmp1(string &line, string & output_line, map <int, EmpiricalPdf *> & pdf_rg, float Log10RatioUb){
  const int adj_cnt = 10;
  float *log10_gp = new float[3];
  stringstream ss, ss_out;
  string chrn, meanCov;
  int aluBegin, aluEnd, midCnt, clipCnt, unknowCnt;
  ss.clear(); ss.str(line); 
  ss >> chrn >> aluBegin >> aluEnd >> meanCov >> midCnt >> clipCnt >> unknowCnt ;
  
  if (midCnt > adj_cnt or clipCnt > adj_cnt) {
    float bias = (aluEnd - aluBegin + 200) / 200.; // read_len = 100
    cout << "bias adj " << bias << endl;
    bias /= 2;
    midCnt = (int) (ceil)(midCnt/bias);
  }

  float prob_ub = pow(10, -Log10RatioUb);
  float *gp = new float[3];

  for (int i = 0; i < 3; i++) log10_gp[i] = 0;
  if (midCnt+clipCnt > 0) {
    log10_gp[0] = clipCnt * log10 ( prob_ub ) + midCnt * log10 ( (1 - prob_ub) );
    log10_gp[1] = (midCnt + clipCnt) * log10 (0.5) ; 
    log10_gp[2] = midCnt * log10 ( prob_ub ) + clipCnt * log10 ( (1 - prob_ub) );
  }

  log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
//  cout << "debug1# " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2]; 
//  cout << "debug3# " << log10_gp[0] << " " << log10_gp[1] << " " << log10_gp[2] << endl;   
  //float freq0 = 0.67;  // high FP ratio      
  float freq0 = ( midCnt + 1 )/(float)(midCnt + clipCnt + 2); // 1 and 2 are psudo count
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
      log10_gp[0] += log10 (p_y);
      log10_gp[1] += log10 ((freq0 * p_y + (1 - freq0) * p_z)) ;
      log10_gp[2] += log10 (p_z);
    }
  }
  bool use_this_line = false;
  if ( !p00_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    ss_out << chrn << " " << aluBegin << " " << aluEnd << " " << midCnt << " " << clipCnt << " " << unknowCnt ;
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
    ss_out << " " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2]; 
    output_line = ss_out.str();
    use_this_line = true;
  }
  delete log10_gp;
  delete gp;    
  return use_this_line;
}

bool binomial_del_tmp1(string &line, string & output_line, map <int, EmpiricalPdf *> & pdf_rg, float Log10RatioUb){
  const int insertLen_Ratio = 5;
  float *log10_gp = new float[3];
  stringstream ss, ss_out;
  string chrn, meanCov;
  int aluBegin, aluEnd, midCnt, clipCnt, unknowCnt;
  ss.clear(); ss.str(line); 
  ss >> chrn >> aluBegin >> aluEnd >> meanCov >> midCnt >> clipCnt >> unknowCnt ;
  float prob_ub = pow(10, -Log10RatioUb);
  float *gp = new float[3];
  int n0 = 0;
  int n1 = 0;
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
      if (p_y >= p_z * insertLen_Ratio) n0 ++;
      else if (p_z >= p_y * insertLen_Ratio) n1 ++;
    }
  }
  n0 += midCnt;
  n1 += clipCnt;
  if (n0 + n1 == 0) 
    return false;

  //cout << "debug2# " << n0 << " " << n1 << endl;

  //  chr21 40968644 40968939 27 61 4 3 0:491 0:506 0:490

  log10_gp[0] = n1 * log10 (prob_ub) + n0 * log10 (1 - prob_ub);
  log10_gp[1] = (n0+n1) * log10(0.5);
  log10_gp[2] = n0 * log10 (prob_ub) + n1 * log10 (1 - prob_ub);
  cout << "debug3# " << log10_gp[0] << " " << log10_gp[1] << " " << log10_gp[2] << endl; 
  
  log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
  bool use_this_line = false;
  if ( !p00_is_dominant(log10_gp, - LOG10_GENO_PROB) ) {
    ss_out << chrn << " " << aluBegin << " " << aluEnd << " " << midCnt << " " << clipCnt << " " << unknowCnt ;
    log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
    ss_out << " " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2]; 
    output_line = ss_out.str();
    use_this_line = true;
  }
  delete gp;    
  delete log10_gp;  
  return use_this_line;
}



void calculate_genoProb(string fn_tmp1, string fn_tmp2, map <int, EmpiricalPdf *> & pdf_rg, int Log10RatioUb){
  // clip_read: evidence for deletion, mid_read: evidence for insertion
  ofstream fout(fn_tmp2.c_str());
  fout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";  
  string line, output_line;
  ifstream fin(fn_tmp1.c_str());
  assert(fin);
  getline(fin, line); // read header
  while (getline(fin, line))
    if (parseline_del_tmp1(line, output_line, pdf_rg, Log10RatioUb))
      fout << output_line << endl;
  fin.close();
  fout.close();
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
  if(!fin) {
    cerr << f_input << " does not exist!\n";
    exit(0);
  }
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

  string path0 = cf_fh.get_conf("file_alu_delete0");
  check_folder_exists(path0);

  float Log10RatioUb = seqan::lexicalCast<float> (cf_fh.get_conf("Log10_RATIO_UB"));

  if (opt == "debug2") {

    string pn = argv[3];
    string line, output_line;
    map <int, EmpiricalPdf *> pdf_rg;    
    string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5  
    read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);

    cout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";    
    //line = "chr21 21903875 21904163 26 48 5 2 0:525 0:534";
    //line = "chr21 41931103 41931403 30 62 6 2 0:514 0:550";
    line = "chr21 33299674 33299968 22 56 0 7 0:496 0:805 0:826 0:461 0:809 0:457 0:834";

    cout << line << endl;
    parseline_del_tmp1(line, output_line, pdf_rg, Log10RatioUb);
    cout << output_line << endl;

    binomial_del_tmp1(line, output_line, pdf_rg, Log10RatioUb);  // biased !!!     
    cout << output_line << endl;
    
  } else if (opt == "debug3") {  // check some region 

    string pn = argv[3];
    string chrn = "chr21";
    int aluBegin = 34791400;
    int aluEnd = 34791714;
    map < seqan::CharString, T_READ> qName_info;  
    
    FastaFileHandler *fasta_fh = new FastaFileHandler(file_fa_prefix + chrn + ".fa", chrn);    
    string bam_input = cf_fh.get_conf("file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";      
    BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bai_input);
    seqan::BamAlignmentRecord record;  

    bam_fh->jump_to_region(chrn, aluBegin-ALU_FLANK, aluEnd + ALU_FLANK);
    ///seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
    //bamStreamOut.header = bamStreamIn.header;

    while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, aluBegin - ALU_FLANK, aluEnd + ALU_FLANK, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_delete_read(record)) continue;
      map <seqan::CharString, T_READ >::iterator qItr = qName_info.find(record.qName);

//      if (seqan::toCString(record.qName) == string("hap:26|sim:1047173")) {
//	///writeRecord(bamStreamOut, record);
//	cout << get_cigar(record) << " " << record.beginPos << " " << record.pNext <<  " " << record.seq << endl;
//      }

      if ( qItr != qName_info.end() and qItr->second != unknow_read) 
	continue;
      T_READ iread = classify_read( record, aluBegin, aluEnd, fasta_fh);    
      if ( has_soft_first(record, CLIP_BP) or ( has_soft_last(record, CLIP_BP) ) ) {
	cout << "soft ";
	debug_print_read(record, cout); 
      }

      if ( iread == useless_read) continue;    
      qName_info[record.qName] = iread;
    }    
    //for (map < seqan::CharString, T_READ>::iterator qi = qName_info.begin(); qi != qName_info.end(); qi ++ ) {
    //  if ( qi->second == unknow_read) 

    bam_fh->jump_to_region(chrn, aluBegin-ALU_FLANK, aluEnd + ALU_FLANK);
    while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, aluBegin - ALU_FLANK, aluEnd + ALU_FLANK, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_delete_read(record)) continue;
      map <seqan::CharString, T_READ >::iterator qItr = qName_info.find(record.qName);
//      if ( qItr != qName_info.end() )
//	if (qItr->second == mid_read) {
//	  cout << "mid_read " ;
//	  debug_print_read(record, cout);
//	} else if (qItr->second == clip_read) {
//	  cout << "clip_read " ;
//	  debug_print_read(record, cout);
//	}
    }    
    delete bam_fh;    

  } else {
    cout << "unknown options !\n";
  }
  
  cout << "time used " << clocki.elapsed() << endl;
  return 0;  
}
