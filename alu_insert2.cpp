#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include "consensus.h"

typedef seqan::Dna5String TSeq;

void read_seq(string pn, string fn, map <string, vector<string> > & pos_seqs, string skip_mk = ""){
  ifstream fin( fn.c_str() );
  assert (fin);
  string line, pos, tmpv;
  stringstream ss;
  getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pos >> tmpv;
    if (tmpv == skip_mk) continue;
    pos_seqs[pos].push_back( replace_str0_str(line, pn, pos));
  }
  fin.close();
}

void parse_line_aludb(string line, string &pn, int &beginPos, int & sameRC, int & read_len, string & qName){
  stringstream ss;
  string alu_type, rid;
  ss.str(line);
  ss >> pn >> alu_type >> rid >> beginPos >> sameRC >> read_len >> qName;
}

void read_alumate(ifstream & fin, seqan::StringSet <TSeq> & seqs, map < int, string > & rid_chrn, string file_fa, string fn_output) {
  string line, pn, alu_type;
  int rid;
  stringstream ss;
  map <int, vector<string> > rid_aludbs;
  ofstream fout;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pn >> alu_type;
    if ( alu_type.substr(0,3) != "Alu" ) continue;
    ss >> rid;
    rid_aludbs[rid].push_back(line);
  }
  bool write_tmp_out = !fn_output.empty();
  if (write_tmp_out) fout.open(fn_output.c_str());
  for (map <int, vector<string> >::iterator ri = rid_aludbs.begin(); ri != rid_aludbs.end(); ri++) {
    string chrn = "";
    get_mapVal(rid_chrn, ri->first, chrn);
    if (chrn.empty()) continue;
    FastaFileHandler * fasta_fh = new FastaFileHandler(replace_str0_str(file_fa, chrn, "chr0"), chrn);
    for (vector<string>::iterator bi = (ri->second).begin(); bi != (ri->second).end(); bi++ ) {
      string pn, qName;
      int beginPos, sameRC, read_len;
      parse_line_aludb( *bi, pn, beginPos, sameRC, read_len, qName);
      seqan::CharString seq;
      fasta_fh->fetch_fasta_upper(beginPos, beginPos + read_len , seq);          
      if ( sameRC ) reverseComplement(seq);
      appendValue(seqs, (TSeq) seq );      
      if (write_tmp_out)
	fout << pn << " " << seq << " " << qName << endl;
    }
    delete fasta_fh;
  }
  if (write_tmp_out) fout.close();
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

int get_cons_seqlen(string file_input, const string & cons_seq, int minLen_alu_ins, map<string, int> & clipPns) {
  map < pair<string, string>, int >  seqPn_clipLen;
  ifstream fin( file_input.c_str() );
  string line, pn, seq;
  int clipLen;
  if (!fin) return 0;
  stringstream ss;
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
    align_clip_to_consRef( (si->first).first, cons_seq, refBegin, refEnd, si->second );
    if (refBegin) {
      pn_refPos.push_back( make_pair(pn, refBegin) );
      addKey(refBegin_cnt, ceil_by_resolution(refBegin, CLIP_BP_LEFT), 1);
    } 
    else if (refEnd) {
      pn_refPos.push_back( make_pair(pn, refEnd) );
      addKey(refEnd_cnt, ceil_by_resolution(refEnd, CLIP_BP_LEFT), 1);
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
  for (vector < pair<string, int> >::iterator pi = pn_refPos.begin(); pi != pn_refPos.end(); pi++ ) 
    if ( abs( (*pi).second - refBegin) < 2 * CLIP_BP_LEFT or abs( (*pi).second - refEnd) < 2 * CLIP_BP_LEFT )
      addKey(clipPns, (*pi).first, 1);
  return seqlen;
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

   std::set <string> pns_used;
   read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
   string path1 = cf_fh.get_conf( "file_alu_insert1") ;    
   string pathClip = path1 + "clip/";
   string pathCons = path1 + "cons/";
   for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++ ) 
     check_folder_exists( pathCons + *ci + "_pos/") ;
   
   if (opt == "consReads_pns" ) {  // rewrite reads to another folder 
     
     for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
       string pathCons1 = pathCons + *ci + "/" ;
       string pathCons2 = pathCons + *ci + "_pos/" ;
       system( ("rm "+ pathCons2 + "*").c_str() );      
       map <string, vector<string> > pos_seqs;
       for (std::set <string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) {
	 read_seq(*pi, pathCons1 + *pi + ".clip", pos_seqs);
	 read_seq(*pi, pathCons1 + *pi + ".aludb", pos_seqs, "na");
       }
       if (pos_seqs.empty()) continue;      
       for (map <string,  vector<string> >::iterator pi = pos_seqs.begin(); pi != pos_seqs.end() ; pi++ ) {
	 string file_cons1 = pathCons2 + pi->first;
	 ofstream fout1(file_cons1.c_str());
	 for ( vector< string > ::iterator si = (pi->second).begin(); si != (pi->second).end(); si++ )
	   fout1 << *si << endl;
	 sort_file_by_col<string> (file_cons1, 1, true); // sort by pn
	 fout1.close();		
       }
     }
     
  } else if (opt == "consReads_build" ) {  // eg: /home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/cons/chr1_pos/159408

    if (argc < 4 ) {
      cerr << "usage: alu_insert2 config.dk " << opt << " chr1 \n";
      cerr << "usage: alu_insert2 config.dk " << opt << " chr1 pos \n";
      return 1;
    }     
    string chrn = argv[3] ;
    int clipPos = (argc == 5) ? seqan::lexicalCast<int> (argv[4]) : 0;
    
    vector < pair<int, int> > insert_pos;
    if (clipPos)
      insert_pos.push_back( make_pair ( clipPos, clipPos) );
    else
      read_first2col( pathClip + chrn + ".clip_pn", insert_pos, true, 2);   

    if (insert_pos.empty() ) {
      cerr << "ERROR, file not exists: " << pathClip + chrn + ".clip_pn" << endl;
      return 1;
    }

    for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
      int clipLeft = (*pi).first;
      string file_input = pathCons + chrn + "_pos/" + int_to_string(clipLeft);
      string file_output0 = file_input + ".cons_input";
      string file_output = file_input + ".cons_output";
      if ( check_file_size(file_output) > 0) {
	cout << "file " << file_output << " has been built, skipping. \n";
	continue;
      }
      ifstream fin( file_input.c_str() );
      if ( !fin ) continue; 
      map < int, string > rid_chrn;
      get_chrn(cf_fh.get_conf("bam_rid_chrn"), rid_chrn);
      string file_fa = cf_fh.get_conf("file_fa_prefix") + "chr0.fa";
      seqan::StringSet <TSeq>  pos_seqs;    
      read_alumate(fin, pos_seqs, rid_chrn, file_fa, file_output0);
      fin.close();
      
      ofstream fout(file_output.c_str());
      if ( length(pos_seqs) < 2 ) {
	fout << ">" << clipLeft << "\nAAAAAAA\n";   // too few reads for building consensus
      } else if ( length(pos_seqs) > 999 ) {
	fout << ">" << clipLeft << "\nTTTTTTT\n";   // too many reads, need to check again !!
      }	else {
          seqan::StringSet<TSeq> consensus;
          int multi_align_fail = compute_consensus(consensus, pos_seqs);	
          cout <<  chrn << " " << clipLeft << ", multiple read alignment failed ? "<< multi_align_fail << endl;
          size_t max_len = length(consensus[0]);
          int max_id = 0;
          for ( size_t i = 1; i < length(consensus); ++i) 
	    if (length(consensus[i]) > max_len) {
	      max_len = length(consensus[i]);
	      max_id = i;
	    }
          fout << ">" << clipLeft << "\n"<< consensus[max_id] << "\n";  
          for ( size_t i = 0; i < length(consensus); ++i)  
	    fout << ">" << clipLeft << "_" << i << "\n"<< consensus[i] << "\n";  // ok, i might want to debug, why so many contigs ?	
      }
      fout.close();     
    }   
    
   } else if (opt == "consReads_chr" ) {  

    if (argc < 4 ) {
      cerr << "usage: alu_insert2 config.dk " << opt << " chr1 \n";
      return 1;
    }     

    string chrn = argv[3] ;
    int minLen_alu_ins = seqan::lexicalCast <int> (cf_fh.get_conf("minLen_alu_ins"));
    vector < pair<int, int> > insert_pos;
    read_first2col( pathClip + chrn + ".clip_pn", insert_pos, true, 2);   
    ofstream fout( (pathCons + chrn + "_cons_seqlen").c_str() );
    fout << "clipLeft len_insertion others\n";   // default 0 if failed 
    for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
      int clipLeft = (*pi).first;
      string file1 = pathCons + chrn + "_pos/" + int_to_string(clipLeft);
      string file2 = file1 + ".cons_output";
      string cons_seq;
      if ( !get_cons_seq ( cons_seq, file2, int_to_string(clipLeft), minLen_alu_ins ) ) {
	fout << clipLeft << " 0 " << cons_seq << endl;
      } else {  // use clipreads to find insertion length 
	map<string, int> clipPns;
	int cons_seqlen = get_cons_seqlen ( file1, cons_seq, minLen_alu_ins, clipPns ); 
	fout << clipLeft << " " << cons_seqlen <<  " " << length(cons_seq) ;
	for (map<string,int>::iterator ci = clipPns.begin(); ci != clipPns.end(); ci++ )
	  fout << " " <<  ci->first << ":" << ci->second;
	fout << endl;
      }
    }
    fout.close();
    cout << "output to " << pathCons + chrn + "_cons_seqlen\n" ;
    
   } else {
     
     cout << "unknown option ! \n";
     return 1;
     
   }
  
   cout << "total time used " << clocki.elapsed() << endl;
   clocki.restart();  
   return 0;
 }

