#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include "delete_utils.h"
#include "consensus.h"

typedef seqan::Dna5String TSeq;

void parse_line_aludb(string line, string &pn, int &beginPos, int & sameRC, int & read_len, string & qName){
  stringstream ss;
  string alu_type, rid;
  ss.str(line);
  ss >> pn >> alu_type >> rid >> beginPos >> sameRC >> read_len >> qName;
}

bool write_cons_input(string fn_input, map < int, string > & rid_chrn, string file_fa, string fn_output, int cut_bp) {
  ifstream fin( fn_input.c_str() );
  if ( !fin ) return false;
  ofstream fout(fn_output.c_str());;
  string line, pn, tmpv;
  int rid;
  stringstream ss;
  map <int, vector<string> > rid_aludbs;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pn >> tmpv;
    if ( tmpv.substr(0,3) != "Alu" ) {
      fout << pn << " " << tmpv << " clipreads\n";
    }  else {
      ss >> rid;
      rid_aludbs[rid].push_back(line);
    }
  }
  fin.close();
  for (map <int, vector<string> >::iterator ri = rid_aludbs.begin(); ri != rid_aludbs.end(); ri++) {
    string chrn = "";
    get_mapVal(rid_chrn, ri->first, chrn);
    if (chrn.empty()) continue;
    FastaFileHandler * fasta_fh = new FastaFileHandler(replace_str0_str(file_fa, chrn, "chr0"), chrn);
    for (vector<string>::iterator bi = (ri->second).begin(); bi != (ri->second).end(); bi++ ) {
      string pn, qName;
      int beginPos, sameRC, read_len;
      parse_line_aludb( *bi, pn, beginPos, sameRC, read_len, qName);
      seqan::CharString seq0;
      fasta_fh->fetch_fasta_upper(beginPos, beginPos + read_len , seq0);          
      if ( sameRC ) reverseComplement(seq0);
      int seqlen = length(seq0) ; 
      fout << pn << " " << infix( seq0, cut_bp, seqlen - cut_bp) << " " << qName << endl;
    }
    delete fasta_fh;
  }
  fout.close();
  return true;
}

bool read_cons_input(string fn_input,  seqan::StringSet <TSeq> & pos_Tseqs) {
  clear(pos_Tseqs);
  ifstream fin( fn_input.c_str() );
  if ( !fin ) return false;
  std::set<string> pos_seqs;
  string line, pn, tmpv;
  stringstream ss;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pn >> tmpv;
    pos_seqs.insert(tmpv);
  }
  for (std::set<string>::iterator pi = pos_seqs.begin(); pi != pos_seqs.end(); pi++ ) 
    appendValue(pos_Tseqs, (TSeq) *pi);
  return (!pos_seqs.empty());
}

bool read_cons_input2(string fn_input,  seqan::StringSet <TSeq> & pos_Tseqs) {
  clear(pos_Tseqs);
  ifstream fin( fn_input.c_str() );
  if ( !fin ) return false;
  std::set<string> pos_seqs;
  string line, pn, tmpv;
  stringstream ss;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pn >> tmpv;
    pos_seqs.insert(tmpv);
  }
  for (std::set<string>::iterator pi = pos_seqs.begin(); pi != pos_seqs.end(); pi++ ) 
    appendValue(pos_Tseqs, (TSeq) *pi);
  return (!pos_seqs.empty());
}

void write_cons_output(int clipLeft, seqan::StringSet <TSeq> & pos_Tseqs, ofstream &fout1) {
  seqan::StringSet<TSeq> consensus;
  int multi_align_fail = compute_consensus(consensus, pos_Tseqs);	
  if ( multi_align_fail )
    cout << clipLeft << " multiple read alignment failed ? \n";
  size_t max_len = length(consensus[0]);
  int max_id = 0;
  for ( size_t i = 1; i < length(consensus); ++i) 
    if (length(consensus[i]) > max_len) {
      max_len = length(consensus[i]);
      max_id = i;
    }
  fout1 << consensus[max_id] << "\n";  
  for ( size_t i = 0; i < length(consensus); ++i)  
    fout1 << ">" << clipLeft << "_" << i << "\n"<< consensus[i] << "\n";  // ok, i might want to debug, why so many contigs ?	
  
}

void get_aluType( string file_alu, map<string, string> & qname_aluType) {
  qname_aluType.clear();
  ifstream fin;
  string line, pn, alu_type, tmpv1, tmpv2, tmpv3, tmpv4, qname;
  stringstream ss;
  fin.open( file_alu.c_str() );
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pn >> alu_type; 
    if ( alu_type.substr(0,3) != "Alu" ) continue;
    ss >> tmpv1 >> tmpv2 >> tmpv3 >> tmpv4 >> qname;
    qname_aluType[qname] = alu_type;
    //qname_aluType[qname] = alu_type + "_" +tmpv3;
  }
  fin.close();
}

bool get_cons_seq (string &cons_seq, string file_input, string seq_name, int minLen_alu_ins ) {
  seqan::FaiIndex faiIndex;
  build(faiIndex, file_input.c_str() ); 
  unsigned idx = 0;
  getIdByName(faiIndex, seq_name, idx);
  seqan::CharString seq;
  readSequence(seq, faiIndex, idx);
  cons_seq = toCString(seq); 
  return  (int)length(seq) >= minLen_alu_ins;
}

int get_cons_insertLen(string file_clip,  string file_alu, const string & cons_seq, int minLen_alu_ins, float consensus_freq, map<string, int> & pn_newAluCnt, int & cnt0, int & cnt1) {
  map < pair<string, string>, int >  seqPn_clipLen;
  ifstream fin;
  string line, pn, seq, qname;
  int clipLen;
  stringstream ss;
  // pns with alu reads
  fin.open( file_alu.c_str() );
  cnt0 = 0; cnt1 = 0;
  map<string, string> qname_aluType;
  get_aluType(file_clip, qname_aluType);
  //map<string, int> aluType_cnt;
  while ( getline(fin, line) ) {
    cnt0++;
    ss.clear(); ss.str( line );
    ss >> pn >> seq >> qname;
    string alu_type = qname_aluType[qname];
    if ( align_alu_to_consRef(seq, cons_seq, 0.15, alu_type))  {
      addKey(pn_newAluCnt, pn, 1); 
      //addKey(aluType_cnt, alu_type, 1);
      cnt1++;
    } else 
      align_alu_to_consRef(seq, cons_seq, 0.15, "print");
  }
  fin.close();
  cout << "align alu " << cnt0 << ", with success of " << cnt1 << endl;

  // pns with clip reads
  fin.open( file_clip.c_str() );
  if (!fin) return 0;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pn >> seq;
    if ( seq.substr(0,3) == "Alu" ) continue;
    ss >> clipLen;
    seqPn_clipLen[make_pair(seq, pn)] = clipLen;
  }
  fin.close();

  if (seqPn_clipLen.empty()) return 0;
  map <int, int> refBegin_cnt;
  map <int, int> refEnd_cnt;
  vector < pair<string, int> > pn_refPos;
  int refBegin, refEnd;
  for (map < pair<string, string>, int >::iterator si = seqPn_clipLen.begin(); si != seqPn_clipLen.end(); si++) {
    refBegin = -200; refEnd = -200;
    align_clip_to_LongConsRef( (si->first).first, cons_seq, refBegin, refEnd, si->second );    

    cout <<  "align " << (si->first).second << " " << refBegin << " " << refEnd << endl;

    if (refBegin > -200) {
      pn_refPos.push_back( make_pair(pn, refBegin) );
      addKey(refBegin_cnt, round_by_resolution(refBegin, CLIP_BP_LEFT), 1);
    } 
    else if (refEnd > -200 ) {
      pn_refPos.push_back( make_pair(pn, refEnd) );
      addKey(refEnd_cnt, round_by_resolution(refEnd, CLIP_BP_LEFT), 1);
    }
  }
  if (refBegin_cnt.empty() or refEnd_cnt.empty()) 
    return 0;

  multimap<int, int> cnt_refBegin = flip_map(refBegin_cnt);
  multimap<int, int> cnt_refEnd = flip_map(refEnd_cnt);
  refEnd = (cnt_refEnd.rbegin())->second;
  refBegin = (cnt_refBegin.rbegin())->second;
  int seqlen = refEnd - refBegin;
  if (seqlen < minLen_alu_ins)
    return 0;
  // re_check clip_reads
  int clip_pass_cnt = 0;
  for (vector < pair<string, int> >::iterator pi = pn_refPos.begin(); pi != pn_refPos.end(); pi++ ) {
    if ( abs( (*pi).second - refBegin) < 2 * CLIP_BP_LEFT or abs( (*pi).second - refEnd) < 2 * CLIP_BP_LEFT )  
      clip_pass_cnt ++;
  }
  if (clip_pass_cnt >= consensus_freq * seqPn_clipLen.size() )
    return seqlen;
  return 0;
}

void read_consensus (string fn_input, float consensus_freq, map<int, int> & pos_consLen, string pn, map<int, int > & pos_aluCnt) {
  pos_consLen.clear(); 
  pos_aluCnt.clear();
  ifstream fin( fn_input.c_str() );
  assert (fin);
  string line, tmpv, _pn, _cnt;
  int pos, consLen, n1, n2, cnt;
  stringstream ss;
  getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pos >> tmpv;
    if ( ss >> consLen) { 
      ss  >> n1 >> n2;
      if ( n2/(float)n1 < consensus_freq) continue;
      pos_consLen[pos] = consLen;
      while (ss >> tmpv) {
	split_by_sep(tmpv, _pn, _cnt, ':');
	seqan::lexicalCast2(cnt, _cnt);
	if ( _pn == pn ) {
	  pos_aluCnt[pos] = cnt;
	  break;
	}
      }
    }
  }
}

void empty_alumates_only(map <int, EmpiricalPdf *> & pdf_rg, string chrn, map <int, int> & pos_aluCnt, list< IntString> & rows_list, map<int, int> & pos_consLen) {
  for (map <int, int>::iterator bi = pos_aluCnt.begin(); bi != pos_aluCnt.end(); bi ++ ) {
    stringstream ss;
    ss << chrn << " " << bi->first << " " << bi->first << " " << pos_consLen[bi->first] << " 0 0 0";
    string line0 = ss.str();
    string output_line;
    if (parseline_del_tmp0(line0, output_line, pdf_rg, bi->second, 299 ))
      rows_list.push_back( make_pair(bi->first, output_line) );
  }	
  pos_aluCnt.clear();
  pos_consLen.clear();
}


void write_tmp2_chrn( list< IntString> & rows_list, ofstream & fout) {
  rows_list.sort(compare_IntString);
  for (list< IntString>::iterator ri = rows_list.begin(); ri!=rows_list.end(); ri++) 
    fout << (*ri).second << endl;
  rows_list.clear();
}

void write_tmp2(string pn, string f1_tmp0, string f2_tmp2, map <int, EmpiricalPdf *> & pdf_rg, string fn_alu0, float consensus_freq ) {
  ofstream fout(f2_tmp2.c_str());
  fout << "chr insertBegin insertLen midCnt clipCnt unknowCnt 00 01 11\n";  
  ifstream fin(f1_tmp0.c_str());
  assert(fin);
  stringstream ss;
  string line, output_line, chrn, fn_alu, tmpv1;
  string pre_chrn = "";
  int pos, len_alucon;
  map<int, int> pos_consLen; 
  map<int, int> pos_aluCnt; 
  list< IntString> rows_list;
  getline(fin, line); // read header
  while (getline(fin, line)) {
    ss.clear(); ss.str(line); 
    ss >> chrn >> pos >> tmpv1 >> len_alucon;
    if ( chrn != pre_chrn) {
      if ( pre_chrn != "") {
	empty_alumates_only(pdf_rg, pre_chrn, pos_aluCnt, rows_list, pos_consLen);
	write_tmp2_chrn(rows_list, fout);
      }
      pre_chrn = chrn;      
      fn_alu = replace_str0_str(fn_alu0, chrn, "chr0");
      read_consensus (fn_alu, consensus_freq, pos_consLen, pn, pos_aluCnt);
    }
    
    int cnt_alumate = 0;
    get_mapVal(pos_aluCnt, pos, cnt_alumate);
    if (cnt_alumate == 0 )
      continue;   // failed to build consensus
    
    pos_aluCnt.erase( pos );
    if (parseline_del_tmp0(line, output_line, pdf_rg, cnt_alumate, pos_consLen[pos]))
      rows_list.push_back( make_pair(pos, output_line) );
  } 
  empty_alumates_only(pdf_rg, chrn, pos_aluCnt, rows_list, pos_consLen);
  write_tmp2_chrn(rows_list, fout);
  fin.close();  
  fout.close();
  
}

int main( int argc, char* argv[] )
 {
   string config_file = argv[1];
   string opt =  argv[2];
   if (argc < 3) exit(1);
   
   boost::timer clocki;    
   clocki.restart();

   ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
   vector<string> chrns;  
   for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
   chrns.push_back("chrX");
   chrns.push_back("chrY");  

   map <int, string> ID_pn;
   get_pn(cf_fh.get_conf( "file_pn"), ID_pn);
   std::set <string> pns_used;
   read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
   string path1 = cf_fh.get_conf( "file_alu_insert1") ;    
   string pathDel0 = path1 + "fixed_delete0/";
   string pathDel1 = path1 + "cons_delete0/";
   check_folder_exists(pathDel1);
   check_folder_exists(pathDel1+"tmp2s/");
   string pathClip = path1 + "clip/";
   string pathCons = path1 + "cons/";
   for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++ )  
     check_folder_exists( pathCons + *ci + "_pos/") ;
   int minLen_alu_ins = seqan::lexicalCast <int> (cf_fh.get_conf("minLen_alu_ins"));
   float consensus_freq = seqan::lexicalCast <float> (cf_fh.get_conf("consensus_freq"));
   
   if (opt == "consReads_build" ) {  // eg: /home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/cons/chr1_pos/159408
     
     if (argc < 4 ) {
       cerr << "usage: alu_insert2 config.dk " << opt << " chr1 \n";
       cerr << "usage: alu_insert2 config.dk " << opt << " chr1 pos \n";
       return 1;
     }     
     string chrn = argv[3] ;
     int clipPos = (argc == 5) ? seqan::lexicalCast<int> (argv[4]) : 0;
     
     vector < pair<int, int> > insert_pos;
     if (clipPos) {
       insert_pos.push_back( make_pair ( clipPos, clipPos) );
     } else if (!read_first2col( pathClip + chrn + ".clip_pn", insert_pos, true) ) {
       cerr << "ERROR, file not exists: " << pathClip + chrn + ".clip_pn" << endl;
       return 1;
     }
     
     const int BUILD_CONS_MAX_INPUT = 1000;
     for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
       int clipLeft = (*pi).first;
       string file_input = pathCons + chrn + "_pos/" + int_to_string(clipLeft);
       string file_output0 = file_input + ".cons_input";
       string file_output1 = file_input + ".cons_output";
       
       map < int, string > rid_chrn;
       get_chrn(cf_fh.get_conf("bam_rid_chrn"), rid_chrn);
       string file_fa = cf_fh.get_conf("file_fa_prefix") + "chr0.fa";
       int cut_bp = seqan::lexicalCast<int> (cf_fh.get_conf("cut_bp_aluread"));
       
       if (!write_cons_input(file_input, rid_chrn, file_fa, file_output0, cut_bp))
	continue;

      if ( check_file_size(file_output1) > 0) {
	cout << "file " << file_output1 << " has been built, skipping. \n";
	cout << "remove "<< file_output1 << " to rerun. \n";
	continue;
      }
      seqan::StringSet <TSeq>  pos_Tseqs;    
      if (!read_cons_input(file_output0, pos_Tseqs))
	continue;
      cout << "read first, cnt: " << length(pos_Tseqs) << endl;

      ofstream fout1(file_output1.c_str());
      fout1 << ">" << clipLeft << "\n";
      if ( length(pos_Tseqs) < 2 ) {   // check 34776066
	fout1 << "AAAAAAA\n";   // too few reads for building consensus
      } else if ( (int)length(pos_Tseqs) > BUILD_CONS_MAX_INPUT ) {
	read_cons_input2(file_output0, pos_Tseqs); // select reads while reading 
	cout << "read again, cnt: " << length(pos_Tseqs) << endl;
	if ( (int)length(pos_Tseqs) < BUILD_CONS_MAX_INPUT ) 
	  write_cons_output(clipLeft, pos_Tseqs, fout1); 
	else 
	  fout1 << "TTTTTTT\n";   // too many reads, need to check again !!
      }	else {
	write_cons_output(clipLeft, pos_Tseqs, fout1);
      }
      fout1.close();   

      // get sequence length
      string file_output2 = file_input + ".cons_seqlen";
      ofstream fout2( file_output2.c_str() );
      string cons_seq;
      if ( !get_cons_seq ( cons_seq, file_output1, int_to_string(clipLeft), minLen_alu_ins ) ) {
	fout2 << clipLeft << " 0 " << cons_seq << endl;  // failed to build consensus
      } else {  // use clipreads to find insertion length 
	map<string, int> pn_newAluCnt;   // count of clip and alu counts for each pn 
	int cnt_alu_db, cnt_alu_aligned;
	int insertLen = get_cons_insertLen ( file_input, file_output0, cons_seq, minLen_alu_ins, consensus_freq, pn_newAluCnt, cnt_alu_db, cnt_alu_aligned );
	fout2 << clipLeft << " " << insertLen <<  " " << length(cons_seq) << " " << cnt_alu_db << " " << cnt_alu_aligned  ;
	for (map<string,int>::iterator ci = pn_newAluCnt.begin(); ci != pn_newAluCnt.end(); ci++ )
	  fout2 << " " <<  ci->first << ":" << ci->second;
	fout2 << endl;
      }
      fout2.close();
    }   
    
   } else if (opt == "delete0_pn" ) {  // only use loci that have consensus sequence

     map <int, EmpiricalPdf *> pdf_rg;    
     vector <string> pns;
     read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns);
     for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi ++ ) {
       string pn = *pi;
       read_pdf_pn(cf_fh.get_conf("file_dist_prefix"), pn, cf_fh.get_conf("pdf_param"), pdf_rg);
       string f1_tmp0 = pathDel0 + "tmp0s/" + pn + ".tmp0";
       string f2_tmp2 = pathDel1 + pn + ".tmp2";
       string fn_alu0 = pathCons + "chr0_cons_seqlen";
       write_tmp2( pn, f1_tmp0, f2_tmp2, pdf_rg, fn_alu0, consensus_freq);	    	 
       move_files(pathDel1 + "tmp2s/", f2_tmp2);
       pdf_rg.clear();
     }
     
     // write vcf file
     string path_input = pathDel1 + "tmp2s/";
     string tmp_file_pn = path_input + *(pns.begin()) + ".tmp2";
     int col_idx =  get_col_idx(tmp_file_pn, "00");
     assert (col_idx == 7 );
     string fn_pos = pathDel1 + int_to_string( pns.size()) + ".pos";
     filter_by_llh_noPrivate(path_input, ".tmp2", fn_pos, pns, chrns, col_idx);
     string fn_vcf = pathDel1 + int_to_string( pns.size()) + ".vcf";  
     combine_pns_vcf_noPrivate(path_input, ".tmp2", fn_vcf, pns, chrns, col_idx);  
     
   } else if (opt == "debug" ) {  // only use loci that have consensus sequence
     
     string chrn = "chr21";
     int clipLeft = 10163489;
     string file1 = pathCons + chrn + "_pos/" + int_to_string(clipLeft);
     string file2 = file1 + ".cons_output";
     string cons_seq;
     get_cons_seq ( cons_seq, file2, int_to_string(clipLeft), minLen_alu_ins );
     map<string, int> pn_newAluCnt;   
     string file3 = file1 + ".cons_input";
     //int cnt_alu_db, cnt_alu_aligned;
     //int cons_seqlen = get_cons_seqlen ( file1, file3, cons_seq, minLen_alu_ins, pn_newAluCnt, cnt_alu_db, cnt_alu_aligned, consensus_freq ); 
     //cout << clipLeft << " " << cons_seqlen <<  " " << length(cons_seq) << " " << cnt_alu_db << " " << cnt_alu_aligned << endl ;
     
   } else {
     
     cout << "unknown option ! \n";
     return 1;
     
   }
  
   cout << "total time used " << clocki.elapsed() << endl;
   clocki.restart();  
   return 0;
 }

