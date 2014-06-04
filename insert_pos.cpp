#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"

void read_file_pn_used(string fn, set <string> & pns_used) {
  ifstream fin(fn.c_str()) ;
  assert(fin);
  string pn;
  while (fin >> pn) pns_used.insert(pn);
  fin.close();
}

bool align_to_ref(char left_right, int clipPos,  int align_len, seqan::BamAlignmentRecord &record, FastaFileHandler *fasta_fh, ofstream &fout, string  header) {
  // modified based on pn_with_insertion() from previous commits
  int score;
  seqan::CharString ref_fa;
  bool realign_ok = false;
  int read_len = length(record.seq);
  if (left_right == 'R') {
    fasta_fh->fetch_fasta_upper(clipPos, clipPos + align_len, ref_fa);    
    if (global_align_insert( hasFlagRC(record), infix(record.seq, read_len - align_len, read_len), ref_fa, score, ALIGN_END_CUT, 0.7))
      fout << header ;
    else {
      cout << "#R " << get_cigar(record) << endl;
      debug_print_read(record, cout);
      global_align_insert( hasFlagRC(record), infix(record.seq, read_len - align_len, read_len), ref_fa, score, ALIGN_END_CUT, 0.7, true);
    }
  }
  if (left_right == 'L') {
    fasta_fh->fetch_fasta_upper(clipPos - align_len, clipPos, ref_fa);
    if (global_align_insert( hasFlagRC(record), infix(record.seq, read_len - align_len, read_len), ref_fa, score, ALIGN_END_CUT, 0.7))
      fout << header ;
    else {
      cout << "#L " << get_cigar(record) << endl;
      debug_print_read(record, cout);
      global_align_insert( hasFlagRC(record), infix(record.seq, read_len - align_len, read_len), ref_fa, score, ALIGN_END_CUT, 0.7, true);
    }
  }
  
//	     << get_cigar(record) << " " << hasFlagRC(record) << " " << record.seq << endl;	    

  return realign_ok;
}

bool clipreads_at_insertPos(string pn, string chrn, BamFileHandler *bam_fh, FastaFileHandler *fasta_fh,  string fin_pos, string fout_reads) {
  ifstream fin( (fin_pos).c_str());
  if (!fin) return false;
  string line, numAluRead, _pn;
  getline(fin, line); // skip header
  ofstream fout(fout_reads.c_str());
  fout << "chrn region_begin region_end left_right adj_pos qName beginPos endPos cigar pNext seq_len seq\n" ;
  stringstream ss, ss_header;
  int region_begin, region_end, refBegin, refEnd;
  
  while (getline(fin, line)) {
    ss.clear(); ss.str( line );
    ss >> numAluRead >> region_begin >> region_end;
    bool pn_has_alumate = false;
    while ( ss >> _pn ) 
      if ( _pn == pn) { pn_has_alumate=true; break;}
    if (!pn_has_alumate) continue;
    
    refPos_by_region(region_begin, region_end, refBegin, refEnd);  
    if (! bam_fh->jump_to_region(chrn, refBegin, refEnd) ) continue;

    ss_header.clear(); ss_header.str("");
    ss_header << chrn << " " << region_begin << " " << region_end << " " ;    

    seqan::BamAlignmentRecord record;
    while ( true ) {
      string read_status = bam_fh->fetch_a_read(chrn, refBegin, refEnd, record);
      if (read_status == "stop" ) break;
      if (read_status == "skip" or !QC_delete_read(record)) continue;  // search for clip reads 
      int beginPos = record.beginPos;
      int endPos = beginPos + getAlignmentLengthInRef(record);
      int refBegin = region_begin - 200;
      int refEnd = region_end + 200;

      int clipPos, align_len;
      seqan::CharString ref_fa;
      list <char> cigar_opts;
      list <int> cigar_cnts;
      
      if (has_soft_last(record, CLIP_BP) and endPos >= region_begin - FLANK_REGION_LEN and endPos < region_end + FLANK_REGION_LEN) {
	clipPos = endPos;
	if (length(ref_fa) == 0) {
	  fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
	  parse_cigar(get_cigar(record), cigar_opts, cigar_cnts);
	}
	if (clipLeft_move_right( record.seq, ref_fa, cigar_cnts, refBegin, clipPos, align_len) and 
	    align_to_ref('L', clipPos, align_len, record, fasta_fh, fout, ss_header.str()) )
	  continue;
      }
      
      if ( has_soft_first(record, CLIP_BP) and beginPos >= region_begin - FLANK_REGION_LEN and beginPos < region_end + FLANK_REGION_LEN) {
	clipPos = beginPos;
	if (length(ref_fa) == 0) {
	  fasta_fh->fetch_fasta_upper(refBegin, refEnd, ref_fa);
	  parse_cigar(get_cigar(record), cigar_opts, cigar_cnts);
	}
	if (clipRight_move_left( record.seq, ref_fa, cigar_cnts, refBegin, clipPos, align_len) and 
	    align_to_ref('R', clipPos, align_len, record, fasta_fh, fout, ss_header.str()) )
	  continue;
      }
    }
  }
  fin.close();
  fout.close();
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
  
  set <string> pns_used;
  read_file_pn_used(cf_fh.get_conf( "file_pn_used"), pns_used); // some pn is ignored, due to too many reads
  
  if (opt == "clipReads_by_pn") { // clip reads for each pn
    assert (argc == 4);
    string pn = get_pn(cf_fh.get_conf( "file_pn"), seqan::lexicalCast<int> (argv[3]));      
    if ( pns_used.find(pn) == pns_used.end() ){
      cerr << pn << " is not used due to high(strange) coverage in potential regions\n";
      return 0;
    }
    cerr << "reading " << pn << endl;
    string bam_input = cf_fh.get_conf( "file_bam_prefix") + pn + ".bam";
    BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bam_input+".bai"); 
    string file_fa = cf_fh.get_conf("file_fa_prefix");    
#ifdef DEBUG_MODE   // only look at chr1 for now 
    chrns.clear();
    chrns.push_back("chr1");
#endif             
    string fn_suffix = get_name_suffix(seqan::lexicalCast<float> (cf_fh.get_conf("freq_min")), seqan::lexicalCast<float> (cf_fh.get_conf("freq_max")));
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {      
      string fin_pos = path1 + "insert_pos." + *ci + fn_suffix;
      string fout_reads = cf_fh.get_conf( "file_clip_reads") + *ci + "/" + pn + fn_suffix;
      FastaFileHandler * fasta_fh = new FastaFileHandler(file_fa + *ci + ".fa", *ci);
      clipreads_at_insertPos(pn, *ci, bam_fh, fasta_fh, fin_pos, fout_reads);
      delete fasta_fh;
    }
    delete bam_fh;
    
  } else if ( opt == "combine_clipReads" ) { // for each insertion region, combine reads from different pns 

    
  } else if ( opt == "debug" ) {  
    
  }
  return 0;
}
