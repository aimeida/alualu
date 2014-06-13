#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"

class ALUREAD_INFO {
public:
  seqan::CharString qName;
  int clipLeft, pos;
  string pn;
  bool sameRC;
  ALUREAD_INFO(seqan::CharString & qn, int cl, int p, string pn, bool t) : qName(qn), clipLeft(cl), pos(p), pn(pn), sameRC(t) {}
};

void read_file_pn_used(string fn, std::set <string> & pns_used) {
  ifstream fin(fn.c_str()) ;
  assert(fin);
  string pn;
  while (fin >> pn) pns_used.insert(pn);
  fin.close();
}

bool align_clip_to_ref(char left_right, int adj_clipPos,  int clipPos, int align_len, seqan::BamAlignmentRecord &record, FastaFileHandler *fasta_fh, ofstream &fout, string  header) {
  // modified based on pn_with_insertion() from previous commits, only keep clip read with very good mapping quality 
  int ref_plus_bp = 10; // allow small indels
  int score;
  seqan::CharString ref_fa;
  seqan::CharString read_clip_fa;
  int read_len = length(record.seq);
  if (left_right == 'R') {
    fasta_fh->fetch_fasta_upper(adj_clipPos, adj_clipPos + align_len + ref_plus_bp, ref_fa);    
    read_clip_fa = infix(record.seq, read_len - align_len, read_len);
  } else if (left_right == 'L') {
    fasta_fh->fetch_fasta_upper(adj_clipPos - align_len - ref_plus_bp, adj_clipPos, ref_fa);
    read_clip_fa = infix(record.seq, 0, align_len);
  }

//      debug_print_read(record, cout);
//      cout << "#R " << get_cigar(record) << " " << adj_clipPos << " " << adj_clipPos + align_len << endl;
//      cout << record.seq << endl;
//      fasta_fh->fetch_fasta_upper(adj_clipPos - 10, adj_clipPos + align_len + 10, ref_fa);    
//      cout << infix(ref_fa, 0, 10) << " " <<  infix(ref_fa, 10, 10 + align_len) << " " << infix(ref_fa, 10 + align_len, align_len + 20) << endl;
//      global_align_insert( hasFlagRC(record), infix(record.seq, read_len - align_len, read_len), ref_fa, score, ALIGN_END_CUT, 0.7, true);

  if (global_align_insert( hasFlagRC(record), read_clip_fa, ref_fa, score, ALIGN_END_CUT, 0.7) ) {
    //cout << record.qName << " " << score << endl;
    fout << header << left_right << " " << adj_clipPos << " " << record.qName << " " << clipPos << " " <<  get_cigar(record) << " " << record.tLen << " " << record.seq << endl;
    return true;
   }
  return false;
}

bool clipreads_at_insertPos(string pn, string chrn, BamFileHandler *bam_fh, FastaFileHandler *fasta_fh,  string fin_pos, string fout_reads) {
  ifstream fin( (fin_pos).c_str());
  if (!fin) return false;
  string line, numAluRead, _pn;
  getline(fin, line); // skip header
  ofstream fout(fout_reads.c_str());
  fout << "pn region_begin region_end left_right adj_clipPos qName clipPos cigar tLen seq\n" ;
  stringstream ss, ss_header;
  int region_begin, region_end, refBegin, refEnd;
  
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> numAluRead >> region_begin >> region_end;
    bool pn_has_alumate = false;
    while ( ss >> _pn ) 
      if ( _pn == pn) { pn_has_alumate=true; break;}
    if (!pn_has_alumate) continue;

    refBegin = region_begin - 2 *  DEFAULT_READ_LEN;
    refEnd = region_end + 2 * DEFAULT_READ_LEN;    
    if (! bam_fh->jump_to_region(chrn, refBegin, refEnd) ) continue;
    ss_header.clear(); ss_header.str("");
    ss_header << pn << " " << region_begin << " " << region_end << " " ;    

    seqan::BamAlignmentRecord record;
    while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, refBegin, refEnd, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_delete_read(record)) continue;  // search for clip reads 
      int beginPos = record.beginPos;
      int endPos = beginPos + getAlignmentLengthInRef(record);
      int adj_clipPos, align_len;
      seqan::CharString ref_fa;
      list <char> cigar_opts;
      list <int> cigar_cnts;
      
      if (has_soft_last(record, CLIP_BP) and endPos >= region_begin - FLANK_REGION_LEN and endPos < region_end + FLANK_REGION_LEN) {
	adj_clipPos = endPos;
	if (length(ref_fa) == 0) {
	  fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
	  parse_cigar(get_cigar(record), cigar_opts, cigar_cnts);
	}
	if (clipLeft_move_right( record.seq, ref_fa, cigar_cnts, refBegin, adj_clipPos, align_len) and 
	    align_clip_to_ref('L', adj_clipPos, endPos, align_len, record, fasta_fh, fout, ss_header.str()) )
	  continue;
      }
      
      if ( has_soft_first(record, CLIP_BP) and beginPos >= region_begin - FLANK_REGION_LEN and beginPos < region_end + FLANK_REGION_LEN) {
	adj_clipPos = beginPos;
	if (length(ref_fa) == 0) {
	  fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
	  parse_cigar(get_cigar(record), cigar_opts, cigar_cnts);
	}
	if (clipRight_move_left( record.seq, ref_fa, cigar_cnts, refBegin, adj_clipPos, align_len) and 
	    align_clip_to_ref('R', adj_clipPos, beginPos, align_len, record, fasta_fh, fout, ss_header.str()) )
	  continue;
      }
    }
  }
  fin.close();
  fout.close();
  return true;
}

// for each insertion region, combine reads from different pns 
void combine_clipreads_by_pos(std::set<string> &pns_used, ConfigFileHandler & cf_fh, string fn_suffix, string chrn) {
  map< pair<string, string>, vector<string> > regionPos_lines;
  for (std::set<string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) {
    string file_input_clip = cf_fh.get_conf( "file_clip_reads") + chrn + "/" + *pi + fn_suffix;
    ifstream fin(file_input_clip.c_str());
    if (!fin) {
      cout << "error " << file_input_clip <<  " missing \n";
      continue;
    }
    stringstream ss;
    string line, pn, region_begin, region_end;
    getline(fin, line);  // read header 
    while ( getline(fin, line) ) {
      ss.clear(); ss.str(line);
      ss >> pn >> region_begin >> region_end;
      regionPos_lines[make_pair(region_begin, region_end)].push_back(line); 
    }
    fin.close();
  } 
  map< pair<string, string>, vector<string> >::iterator ri;
  for (ri = regionPos_lines.begin(); ri != regionPos_lines.end(); ri ++) {
    string file_output_clip = cf_fh.get_conf( "file_clip_reads") + chrn + "_pos/" + (ri->first).first + "_" + (ri->first).second;	
    ofstream fout(file_output_clip.c_str() );
    for (vector<string>::iterator ii = (ri->second).begin(); ii != (ri->second).end(); ii++) 
      fout << *ii << endl;
    fout.close();
  }
}

void parse_fn(string fn, int & regionBegin, int & regionEnd) {
  stringstream ss;
  ss.str(fn);
  string token;
  getline(ss, token, '_');
  seqan::lexicalCast2(regionBegin, token);
  getline(ss, token, '_');
  seqan::lexicalCast2(regionEnd, token);
}

bool nonempty_files_sorted(string & path1, map < pair<int, int>, string> & regionPos_fn) {
  DIR *d;
  struct dirent *dir;
  d = opendir(path1.c_str());
  string fn;
  int regionBegin, regionEnd;
  if (d) {
    while ( (dir = readdir(d))!=NULL )  {
      fn = dir->d_name;
      if (fn[0] == '.' or !is_nonempty_file(fn)) continue;
      parse_fn(fn, regionBegin, regionEnd);
      regionPos_fn[ make_pair(regionBegin, regionEnd)] = fn;
    }
    closedir(d);
    return true;
  }
  return false;
}

void addcnt_if_pos_match( map <int, float> & clipLeft_cnt, int pos, int match_offset, float inc_val) {
  int match_key = 0;
  for ( map <int, float>::iterator mi = clipLeft_cnt.begin(); mi != clipLeft_cnt.end(); mi ++) {
    if ( abs(mi->first - pos) <= match_offset ) {
      match_offset = abs(mi->first - pos);
      match_key = mi->first;
    }
  }
  if ( match_key == 0) {
    addKey(clipLeft_cnt, pos, inc_val);
  } else {
    if (match_key > pos) { // previously add clipRight
      float cnt = clipLeft_cnt[match_key] + inc_val;
      clipLeft_cnt.erase(match_key);
      clipLeft_cnt[pos] = cnt;
      // cout << "offset " << match_key << " " << pos << endl;
    } else {
      addKey(clipLeft_cnt, match_key, inc_val);
    }
  }
}

// each pn can have max 2 pos
bool pn_vote_max2pos( vector <pair<char, int> > & pos_of_pn, int & p1, int & p2, float & f1, float &f2){
  p1 = p2 = 0;  
  f1 = f2 = 0;
  int vsize = pos_of_pn.size();
  if ( vsize < 1) return false;
  map <int, float> clipLeft_cnt;
  for (vector <pair<char, int> >::iterator pi = pos_of_pn.begin(); pi != pos_of_pn.end(); pi ++ ) {
    int match_offset =  ((*pi).first == 'L') ? CLIP_BP_LEFT : CLIP_BP_RIGHT ;
    addcnt_if_pos_match(clipLeft_cnt, (*pi).second, match_offset, 1);
  }

  //debugprint_map(clipLeft_cnt);
  multimap<float, int> cnt_clipLeft = flip_map(clipLeft_cnt);
  multimap<float, int>::reverse_iterator li = cnt_clipLeft.rbegin();
  p1 = (li->second);
  f1 = (li->first) / vsize;
  if ( cnt_clipLeft.size() == 1) return true;
  li ++; 
  if ( (f2 = (li->first) / vsize) > 0.3) { // 2nd vote
    p2 = (li->second);
    float f12 = f1 + f2;
    f1 /= f12;
    f2 /= f12;
  } else {
    f1 = 1;
    f2 = 0;
  }
  return true;
}

void region_pos_vote( string fn, int & clipLeft1, int & clipLeft2 ) {
  ifstream fin(fn.c_str());
  assert(fin);
  stringstream ss;
  string line, pn, sleft_right;
  int region_begin, region_end, clipPos;
  map <string, vector <pair<char, int> > > pn_splitori_pos;
  while ( getline(fin, line)) {  
    ss.clear(); ss.str( line );
    ss >> pn >> region_begin >> region_end >> sleft_right >> clipPos;
    pn_splitori_pos[pn].push_back( make_pair( sleft_right[0], clipPos) );
  }
  fin.close();  
  map < int, float> clipLeft_cnt;
  float vsize = 0;  
  for (map <string, vector <pair<char, int> > >::iterator pi = pn_splitori_pos.begin(); pi != pn_splitori_pos.end(); pi++) {
    int p1, p2;
    float f1, f2;
    if (!pn_vote_max2pos(pi->second, p1, p2, f1, f2)) 
      continue;
    //cout << pi->first << " " << p1 << " " << p2 << " " << f1 << " " << f2 << endl;
    assert ( abs (f1 + f2 - 1.) < 1e-5) ;
    // fixme : add weight to the count ?
    if (p1 > 0) vsize += 1;
    if (p1 > 0) addcnt_if_pos_match(clipLeft_cnt, p1, CLIP_BP_RIGHT, f1);
    if (p2 > 0) addcnt_if_pos_match(clipLeft_cnt, p2, CLIP_BP_RIGHT, f2);
  }
  multimap<float, int> cnt_clipLeft = flip_map(clipLeft_cnt);
  multimap<float, int>::reverse_iterator li = cnt_clipLeft.rbegin();
  clipLeft1 = 0, clipLeft2 = 0 ;
  if ( li->first / vsize >= MIN_VOTE_FREQ) clipLeft1 = li->second;
  if ( cnt_clipLeft.size() > 1) {
    li++;
    if ( li->first / vsize >= MIN_VOTE_FREQ) clipLeft2 = li->second;
  }
}

void regions_pos_vote( string path1, string fn_output) {
  map < pair<int,int>,  string> regionPos_fn;
  assert( nonempty_files_sorted( path1, regionPos_fn) );
  ofstream fout(fn_output.c_str());
  fout << "region_begin region_end clipLeft\n";
  for (map < pair <int, int>, string >::iterator ri = regionPos_fn.begin(); ri != regionPos_fn.end(); ri++ ) {
    int clipLeft1, clipLeft2 ;
    //if ( (ri->first).first > 26984) break;
    //cout << path1 + ri->second << endl;
    region_pos_vote( path1 + ri->second, clipLeft1, clipLeft2);
    if (clipLeft1 > 0) {
      fout << (ri->first).first << " " << (ri->first).second << " " << clipLeft1 ;
      if ( clipLeft2 > 0 )
	fout << " " << clipLeft2 ;
      fout << endl;
    }
  }  
}

void exactpos_pns(string fn_input, string path1, string fn_output){
  ifstream fin;
  fin.open(fn_input.c_str());
  assert(fin);
  map< int, vector<string> > clipLeft_fn;
  stringstream ss;
  string line, region_begin, region_end;
  int clipLeft;
  getline(fin,line);
  while (getline(fin,line)) {
    ss.clear(); ss.str( line );
    ss >> region_begin >> region_end;
    while ( ss >> clipLeft) 
      clipLeft_fn[clipLeft].push_back(region_begin+"_"+region_end);
  }
  fin.close();

  cout << "NB: private insertion is ignored !\n";
  ofstream fout(fn_output.c_str());
  fout << "clipLeft clipRight pn(count)\n";
  for (map< int, vector<string> >::iterator pi = clipLeft_fn.begin(); pi != clipLeft_fn.end(); pi++ ) {
    int clipLeft = pi->first;
    map <string, std::set<string> > pn_qName;
    for ( vector<string>::iterator fi = (pi->second).begin(); fi != (pi->second).end(); fi++ ) {
      fin.open( (path1 + *fi).c_str());
      assert(fin);
      string line, pn, tmp1, tmp2, sleft_right, qName;
      stringstream ss;
      int clipPos;
      while (getline(fin, line)) {
	ss.clear(); ss.str( line );
	ss >> pn >> tmp1 >> tmp2 >> sleft_right >> clipPos >> qName;
	int match_offset =  ( sleft_right[0] == 'L') ? CLIP_BP_LEFT : CLIP_BP_RIGHT ;
	if ( abs(clipPos - clipLeft ) <= match_offset )
	  pn_qName[pn].insert(qName);
      }
      fin.close();
    }
    if (pn_qName.size() <= 1) continue;
    fout << clipLeft << " -1" ;  // clipRight unknown for now
    for (map <string, std::set<string> >::iterator pi = pn_qName.begin(); pi != pn_qName.end(); pi++) 
      fout << " " << pi->first << "," << (pi->second).size() ;
    fout << endl;
  }
  fout.close();
}

void write_clipreads_fa(string pn, string chrn, vector< pair<int, int> > & insert_pos, BamFileHandler *bam_fh, string file_cons,  map < int, vector<ALUREAD_INFO> > & rid_alureads) {
  ofstream fout(file_cons.c_str());    
  seqan::BamAlignmentRecord record;
  string read_status;
  for (vector< pair<int, int> >::iterator pi = insert_pos.begin(); pi != insert_pos.end(); pi++) {
    int region_begin = (*pi).first - ALU_FLANK; // no need for larger flank region, if i only want to build consensus
    int region_end = (*pi).second + ALU_FLANK;	
    if (! bam_fh->jump_to_region(chrn, region_begin, region_end)) 
      continue;    
    while ( true ) {
      read_status = bam_fh->fetch_a_read(chrn, region_begin, region_end, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_insert_read(record)) continue;
      int thisEnd = record.beginPos + (int)getAlignmentLengthInRef(record);   
      if (record.rID != record.rNextId or abs(record.tLen) >= DISCORDANT_LEN) {
	if ( ( thisEnd < (*pi).first and !hasFlagRC(record) ) or 
	     ( record.beginPos > (*pi).second and hasFlagRC(record) ) ) {
	  ALUREAD_INFO one_aluRead = ALUREAD_INFO(record.qName, (*pi).first, record.pNext, pn, hasFlagRC(record)==hasFlagNextRC(record) );
	  rid_alureads[record.rNextId].push_back(one_aluRead);
	}
      }
      /// NB: one read can be both clipRead and alu mate 
      if ( ( has_soft_first(record, CLIP_BP) and abs( record.beginPos - (*pi).second) < 30 ) or    	    
	   ( has_soft_last(record, CLIP_BP) and abs( thisEnd - (*pi).first) < 30 ) )
	fout << (*pi).first << " " << record.seq  <<  " clipRead " << pn << endl;
    }
  }
  fout.close();
}

void add_alureads_fa(BamFileHandler *bam_fh, string file_cons,  map < int, vector<ALUREAD_INFO> > & rid_alureads) {
  fstream fout;
  fout.open(file_cons.c_str(), fstream::app|fstream::out);
  seqan::BamAlignmentRecord record;
  for ( map < int, vector<ALUREAD_INFO> >::iterator ri = rid_alureads.begin(); ri != rid_alureads.end(); ri++) {
    for ( vector<ALUREAD_INFO>::iterator ai = (ri->second).begin(); ai != (ri->second).end(); ai++ ) {
      if ( !bam_fh->jump_to_region( ri->first, (*ai).pos - 20, (*ai).pos + 20) ) continue;
      while (bam_fh->fetch_a_read(record) and record.beginPos <= (*ai).pos + 3) {
	if (record.qName == (*ai).qName) {
	  fout << (*ai).clipLeft << " " << record.seq << " aluRead " << (*ai).pn << " " << (*ai).sameRC <<  endl;
	  break;
	}
      }
      //cout << "reads done " << (*ai).qName << " " << (*ai).pos << endl;
    }  
    cout << ri->first << " " << (ri->second).size() << endl;
    (ri->second).clear();
  }
  fout.close();
}

bool read_Tseq(string fn, map <int, vector <string> > & pos_seqs) {
  ifstream fin( fn.c_str() );
  string line, seq, readtype; 
  int pos;
  stringstream ss;
  if (!fin) return false;
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> pos;
    pos_seqs[pos].push_back(line);
  }
  fin.close();
  return true;
}


int main( int argc, char* argv[] )
{
  string config_file = argv[1];
  string opt = argv[2];
  if (argc < 3) exit(1);
  boost::timer clocki;    
  clocki.restart();  

  vector<string> chrns;
  for (int i = 1; i < 23; i++) chrns.push_back("chr" + int_to_string(i));
  chrns.push_back("chrX");
  chrns.push_back("chrY");
  ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
  string path1 = cf_fh.get_conf( "file_alu_insert1") ;    
  string fout_path = cf_fh.get_conf( "file_clip_reads");
  check_folder_exists(fout_path);  
  for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci ++) {
    string tmp_path;
    tmp_path = fout_path + *ci + "_pos/";
    check_folder_exists( tmp_path);
    tmp_path = fout_path + *ci + "/";
    check_folder_exists( tmp_path );      
  }  
  std::set <string> pns_used;
  read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
  string fn_suffix = get_name_suffix(seqan::lexicalCast<float> (cf_fh.get_conf("freq_min")), seqan::lexicalCast<float> (cf_fh.get_conf("freq_max")));

  if (opt == "clipReads_by_pn") { // clip reads for each pn
    assert (argc == 4);
    string pn = get_pn(cf_fh.get_conf( "file_pn"), seqan::lexicalCast<int> (argv[3]));      
    if ( pns_used.find(pn) == pns_used.end() ){
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bam_input+".bai"); 
    string file_fa = cf_fh.get_conf("file_fa_prefix");    
#ifdef DEBUG_MODE   // only look at chr1 for now 
    chrns.clear();
    chrns.push_back("chr1");
#endif             
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {      
      string fin_pos = path1 + "insert_pos." + *ci + fn_suffix;
      string fout_reads = cf_fh.get_conf( "file_clip_reads") + *ci + "/" + pn + fn_suffix;
      FastaFileHandler * fasta_fh = new FastaFileHandler(file_fa + *ci + ".fa", *ci);
      clipreads_at_insertPos(pn, *ci, bam_fh, fasta_fh, fin_pos, fout_reads);
      delete fasta_fh;
    }
    delete bam_fh;
    cout << "output to " << cf_fh.get_conf( "file_clip_reads") << "chr*/" << pn << fn_suffix << endl;
    

  } else if ( opt == "clipReads_by_pos" ) { 
#ifdef DEBUG_MODE   // only look at chr1 for now 
    chrns.clear();
    chrns.push_back("chr1");
#endif
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
      combine_clipreads_by_pos(pns_used, cf_fh, fn_suffix, *ci);
      string path_clip_pos = cf_fh.get_conf( "file_clip_reads") + *ci + "_pos/";
      string file_insert_pos = path1 + "clip/" + *ci + fn_suffix;
      regions_pos_vote(path_clip_pos, file_insert_pos + ".tmp");
      exactpos_pns(file_insert_pos + ".tmp", path_clip_pos, file_insert_pos);
    }

    
  } else if ( opt == "cons_reads_pn" ) {   // write fasta reads for consensus, use clip reads first 
    string pn = get_pn(cf_fh.get_conf( "file_pn"), seqan::lexicalCast<int> (argv[3]));
    assert (argc == 4);

#ifdef DEBUG_MODE   // only look at chr1 for now 
    chrns.clear();
    chrns.push_back("chr1");
#endif

    if ( pns_used.find(pn) == pns_used.end() ){
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bam_input+".bai"); 
    string path_cons = cf_fh.get_conf("file_ins_cons");  // for consensus sequence
    check_folder_exists( path_cons);
    for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
      check_folder_exists( path_cons + *ci);
      string file_cons = path_cons + *ci + "/" + pn;
      vector< pair<int, int> > insert_pos;
      read_first2col( path1 + "clip/" + *ci + fn_suffix, insert_pos );      
      map < int, vector<ALUREAD_INFO>  > rid_alureads;
      write_clipreads_fa(pn, *ci, insert_pos, bam_fh, file_cons, rid_alureads);         
      add_alureads_fa(bam_fh, file_cons, rid_alureads); // get seq from bam file instead of genome.fa
      sort_file_by_col<int> (file_cons, 1, false);
    }   

  }  else if ( opt == "cons_reads_pos" ) { 
    
    string path_cons = cf_fh.get_conf("file_ins_cons");  
    string bin_path = cf_fh.get_conf("bin_path");   
    for (vector <string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++){
      string path2 = path_cons + *ci + "_pos/";
      check_folder_exists( path2 );
      map <int, vector<string> > pos_seqs;
      for (std::set <string>::iterator pi = pns_used.begin(); pi != pns_used.end(); pi ++ ) 
	read_Tseq(path_cons + *ci + "/" + *pi, pos_seqs);
      if (pos_seqs.empty()) continue;      

      string file_pos =  path_cons + *ci + ".pos";
      ofstream fout2(file_pos.c_str());
      cout << "cmds to run, written to " << file_pos << endl;
      for (map <int,  vector<string> >::iterator pi = pos_seqs.begin(); pi != pos_seqs.end() ; pi++ ) {
	string file_cons1 = path2 + int_to_string( pi->first);
	ofstream fout1(file_cons1.c_str());
	for (  vector< string > ::iterator si = (pi->second).begin(); si != (pi->second).end(); si++ )
	  fout1 << *si << endl;
	sort_file_by_col<string> (file_cons1, 4, false); // sort by pn
	fout1.close();		
	fout2 << pi->first << endl;
	/////fout2 << bin_path << "debug/alu_insdel "<< bin_path << "config.dk cons_reads_build " << file_cons1 << endl; 	
      }
      fout2.close();
    }

  } 

  cout << "done time spent " << clocki.elapsed()<< endl;
  return 0;
}
  
