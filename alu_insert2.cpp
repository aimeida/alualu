#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include "consensus.h"

typedef seqan::Dna5String TSeq;

void read_seq(string pn, string fn, map <string, vector<string> > & pos_seqs){
  ifstream fin( fn.c_str() );
  assert (fin);
  string line, pos, seq , info1, info2;
  stringstream ss;
  getline(fin, line);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pos >> seq >> info1 >> info2;
    pos_seqs[pos].push_back( pn + " " +  seq + " " + info1 + " " + info2);
  }
  fin.close();
}

bool read_Tseq(string fn, seqan::StringSet <TSeq> & seqs, bool allow_duplicate ) {
  ifstream fin( fn.c_str() );
  string line, pn, seq;
  map <string, int> read_cnt;
  stringstream ss;
  if (!fin) return false;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pn >> seq ;
    addKey(read_cnt, seq, 1);
  }
  fin.close();  
  for (map <string, int>::iterator si = read_cnt.begin(); si != read_cnt.end(); si++) {
    seq = si->first;
    if (allow_duplicate) 
      for ( int i = 0; i < si->second; i++) appendValue(seqs, (TSeq) seq );
    else
      appendValue(seqs, (TSeq) seq );
  }
  return true;
}

bool read_Clipseq(string fn, map < string, int > & seq_clipLen) {
  ifstream fin( fn.c_str() );
  if (!fin) return false;
  string line, pn, seq, info1, info2;
  stringstream ss;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pn >> seq >> info1 >> info2;
    if ( info1 == "aluRead")
      continue;
    list <char> cigar_opts;
    list <int> cigar_cnts;
    parse_cigar(info1, cigar_opts, cigar_cnts);
    if ( *(cigar_opts.begin()) == 'S' and *(cigar_cnts.begin()) >= CLIP_BP)
      seq_clipLen[seq] = *(cigar_cnts.begin());   // right clip
    else if ( cigar_opts.back() == 'S' and cigar_cnts.back() >= CLIP_BP )
      seq_clipLen[seq] =  - cigar_cnts.back();    // left clip
  }
  fin.close();  
  return !seq_clipLen.empty();
}

string get_cons_seq (string file_input, string seq_name ) {
  const size_t MIN_CONS_LEN = 250;
  seqan::FaiIndex faiIndex;
  build(faiIndex, file_input.c_str() ); 
  unsigned idx = 0;
  getIdByName(faiIndex, seq_name, idx);
  seqan::CharString seq;
  readSequence(seq, faiIndex, idx);
  if ( length(seq) >= MIN_CONS_LEN)
    return toCString(seq);
  return "";
}

int get_cons_seqlen(string file_input, const string & cons_seq ) {
  map < string, int >  seq_clipLen;
  if (cons_seq.empty() or !read_Clipseq( file_input, seq_clipLen) )
    return 0;
  map <int, int> refBegin_cnt;
  map <int, int> refEnd_cnt;
  //  cout << "clipreads cnt " << seq_clipLen.size() << endl;  
  for (map < string, int >::iterator si = seq_clipLen.begin(); si != seq_clipLen.end(); si++) {
    int refBegin, refEnd;
    align_clip_to_consRef( si->first, cons_seq, refBegin, refEnd, si->second );
    if (refBegin) addKey(refBegin_cnt, ceil_by_resolution(refBegin, CLIP_BP_LEFT), 1);
    if (refEnd) addKey(refEnd_cnt, ceil_by_resolution(refEnd, CLIP_BP_LEFT), 1);
    //if (refBegin or refEnd) cout << "##" << refBegin << " " << refEnd << endl;      
  }
  if (refBegin_cnt.empty() or refEnd_cnt.empty()) 
    return 0;
  multimap<int, int> cnt_refBegin = flip_map(refBegin_cnt);
  multimap<int, int> cnt_refEnd = flip_map(refEnd_cnt);
  int seqlen = (cnt_refEnd.rbegin())->second  - (cnt_refBegin.rbegin())->second;
  return seqlen > 200 ? seqlen : 0;
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

#ifdef DEBUG_MODE   // only look at chr1 for now 
  cerr << "!!!! DEBUG_MODE: only chr1 tested\n";
  chrns.clear();
  chrns.push_back("chr1");
#endif       

  std::set <string> pns_used;
  read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
  string path1 = cf_fh.get_conf( "file_alu_insert1") ;    
  string pathClip = path1 + "clip/";
  string pathCons = path1 + "cons/";
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++ ) {
    check_folder_exists( pathCons + *ci + "_pos/") ;
  }


  if (opt == "cons_reads_pns" ) {  // rewrite reads to another folder 
    
    for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
      map <string, vector<string> > pos_seqs;
      for (std::set <string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) 
	read_seq(*pi, pathCons + *ci + "/" + *pi, pos_seqs);
      if (pos_seqs.empty()) continue;      
      for (map <string,  vector<string> >::iterator pi = pos_seqs.begin(); pi != pos_seqs.end() ; pi++ ) {
	string file_cons1 = pathCons + *ci + "_pos/" + pi->first;
	ofstream fout1(file_cons1.c_str());
	for ( vector< string > ::iterator si = (pi->second).begin(); si != (pi->second).end(); si++ )
	  fout1 << *si << endl;
	sort_file_by_col<string> (file_cons1, 1, true); // sort by pn
	fout1.close();		
      }
    }

  } else if (opt == "cons_reads_build" ) {  // eg: /home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/cons/chr1_pos/159408

    if (argc < 4 ) {
      cout << "usage: alu_insert2 config.dk " << opt << " chr1 \n";
      cout << "usage: alu_insert2 config.dk " << opt << " chr1 pos \n";
      return 1;
    }     
    string chrn = argv[3] ;
    int clipPos = (argc == 5) ? seqan::lexicalCast<int> (argv[4]) : 0;
    
    vector < pair<int, int> > insert_pos;
    if (clipPos)
      insert_pos.push_back( make_pair ( clipPos, clipPos) );
    else
      read_first2col( pathClip + chrn, insert_pos, true, 2);   

    bool allow_duplicate = false; 
    for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
      int clipLeft = (*pi).first;
      string file_input = pathCons + chrn + "_pos/" + int_to_string(clipLeft);
      string file_output = file_input + ".cons";
      if ( check_file_size(file_output) > 0) {
	cout << "file " << file_output << " has been built, skipping. \n";
	continue;
      }
      seqan::StringSet <TSeq>  pos_seqs;    
      if ( !read_Tseq(file_input, pos_seqs, allow_duplicate) ) {
	cout << "file " << file_input << " not exists, skipping. \n";
	continue;
      }
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
    
  } else if (opt == "delete0_pn" ) {  
    
    int min_pn = 2;    
    if (min_pn != 1)
      cout << "NB: private insertions are ignored!\n";

    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      string chrn = *ci;
      vector < pair<int, int> > insert_pos;
      read_first2col( pathClip + chrn ,insert_pos, true, min_pn);   
      ofstream fout( (pathCons + chrn + "_cons_seqlen").c_str() );
      fout << "clipLeft insert_len cons_len\n";
      for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
	int clipLeft = (*pi).first;
	string file1 = pathCons + chrn + "_pos/" + int_to_string(clipLeft);
	string file2 = file1 + ".cons";
	string cons_seq = get_cons_seq ( file2, int_to_string(clipLeft) );
	int cons_seqlen = get_cons_seqlen ( file1, cons_seq );
	fout << clipLeft << " " << cons_seqlen <<  " " << length(cons_seq) << endl;
      }
      fout.close();
    }
    
    
  } else {
    
    cout << "unknown option ! \n";
    return 1;

  }
  
  cout << "total time used " << clocki.elapsed() << endl;
  clocki.restart();  
  return 0;
 }

