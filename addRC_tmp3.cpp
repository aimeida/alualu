// add extra RC columns to *tmp3 files 
#define SEQAN_HAS_ZLIB 1
#include "utils.h"

bool find_RC_info(BamFileHandler *bam_fh, int rid, int pos, const string & qname, bool & isRC, bool find_by_pair) {
  seqan::BamAlignmentRecord record;
  if ( !bam_fh->jump_to_region( bam_fh->rID_chrn[rid], pos - 20, pos + 20) )
    return false;
  while (bam_fh->fetch_a_read(record) and record.beginPos <= pos + 3) {
    if (record.qName == qname) {
      isRC = find_by_pair ? hasFlagNextRC(record) : hasFlagRC(record);
      // debug_print_read(record, cout);
      return true;
    }
  }  
  return false;
}


int main( int argc, char* argv[] )
{
  //debug/addRC_tmp3 config.dk 1390-05 chr6 /home/qianyuxx/faststorage/AluDK/outputs/backup_insert_alu0/tmp3s/1390-05.tmp3.chr6  /home/qianyuxx/faststorage/AluDK/outputs/backup_insert_alu0/tmp3s/1390-05.tmp3.chr6.RC  
  boost::timer clocki;    
  clocki.restart();

  string config_file = argv[1];
  string pn = argv[2];
  string fn_input = argv[3];
  string fn_output = argv[4]; 
  assert( argc == 5 );

  ConfigFileHandler cf_fh = ConfigFileHandler(config_file);  
  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  chrns.push_back("chrX");
  chrns.push_back("chrY");
  string bam_input = cf_fh.get_conf("file_bam_prefix") + pn + ".bam";
  BamFileHandler *bam_fh = new BamFileHandler(chrns, bam_input, bam_input + ".bai");
					      
  ifstream fin(fn_input.c_str());
  assert(fin);
  ofstream fout(fn_output.c_str());
  string line, type_flag, qname;
  int this_chr, this_pos, bad_chr, bad_pos;
  stringstream ss;
  getline(fin, line);
  bool isRC;
  fout << line << " this_is_RC\n";
  while ( getline(fin, line) ) {
    ss.clear(); ss.str( line );
    ss >> type_flag >> qname >> this_chr >> this_pos >> bad_chr >> bad_pos;
    bool read_found = false;
    if ( find_RC_info(bam_fh, this_chr, this_pos, qname, isRC, false) ) {
      fout << line << " " << isRC << endl;
      read_found = 1;
    } else if ( find_RC_info(bam_fh, bad_chr, bad_pos, qname, isRC, true) ) {
      fout << line << " " << isRC << endl;
      read_found = 1;
    }
    if (!read_found) 
      cout << "read not found " << qname << endl;
  }
  fin.close();
  fout.close();
  delete bam_fh;   
  cout << "time used " << clocki.elapsed() << endl;
  return 0;					      
}
