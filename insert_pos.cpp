#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

typedef seqan::Dna5String TSeq;
#define NNN_LEN 70;

class SPLIT_INFO {
public:
  int split_to_refBegin; // dist to ref_begin
  bool alu_is_left;
  SPLIT_INFO( size_t p, bool ail) : split_to_refBegin(p), alu_is_left(ail) {}
};

void parse_fa_name(string const &fa_name, size_t &ref_begin, size_t &ref_end, size_t &alu_len){ // eg.nl_95410890_95411310_AluY_299
  stringstream ss(fa_name);
  string token;
  getline(ss, token, '_');
  getline(ss, token, '_');
  seqan::lexicalCast2(ref_begin, token);
  getline(ss, token, '_');
  seqan::lexicalCast2(ref_end, token);
  getline(ss, token, '_');
  getline(ss, token, '_');
  seqan::lexicalCast2(alu_len, token);
}

void get_filename_fastq(string fin_read_map, map <string, set< pair<int, int> > > & pos_of_clippedReads){
  int pre_pos = 0;
  string chrn, line;
  int region_begin, region_end, idx_pn_first;
  stringstream ss;
  ifstream fin(fin_read_map.c_str());
  while (getline(fin, line) ) {
    ss.clear(); ss.str(line);
    ss >> chrn >> region_begin >> region_end >> idx_pn_first;
    if ( (!pre_pos) or (region_begin != pre_pos ) ) 
      pos_of_clippedReads[chrn].insert( make_pair(region_begin, region_end) );
    pre_pos = region_begin;
  }
  fin.close();
}

bool findSplit(size_t & beginPos,size_t &jump_len, seqan::BamAlignmentRecord &record) {
  char operation;
  size_t count;
  if ( length(record.cigar) < 3) return false;
  beginPos = record.beginPos;
  jump_len = 0;
  for (size_t i = 0; i < length(record.cigar); i++) {
    operation = record.cigar[i].operation;
    count = record.cigar[i].count;
    if ( operation == 'M' or operation == 'D') //  ??? or operation == 'S') 
      beginPos += count;
    //fixme:: if ( operation == 'S') cerr << "soft ?? " << record.qName << endl;    
    if ( operation == 'N') {
      jump_len = count;
      break;
    }
  }
  return jump_len;
}

// move split pos further to ref_begin
void moveSplit_toRefBegin(size_t  & split_begin, TSeq const & ref, size_t nnn_begin, size_t jump_len){  
  /*
  cout << "moving ";
  for (size_t i = split_begin; i <= split_begin + 10; i++) cout << ref[i];
  cout << " " ;
  for (size_t i = split_begin + jump_len; i < split_begin + jump_len + 10; i++) cout << ref[i];
  cout << endl;
  */
  while (split_begin < nnn_begin and ref[split_begin] == ref[ split_begin + jump_len] ) split_begin++;
} 

void read_sam(string fin_sam, string fin_ref, vector< SPLIT_INFO *> & split_pos) {
  string line, qName, col2, fa_name, ref_fasta;
  string strand_lr, alu_type;
  stringstream ss;
  vector <string> fa_names;
  ifstream fin(fin_sam.c_str());
  assert(fin);
  while ( getline(fin, line) ) {
    ss.clear(); ss.str(line);
    ss >> qName >> col2 >> fa_name;
    fa_names.push_back(fa_name);
  }  
  fin.close();
  size_t split_begin, jump_len;  // relative to ref (hg18 + NNN + Alu)
  size_t ref_begin, ref_end, alu_len; // relative to hg18
  parse_fa_name(fa_name, ref_begin, ref_end, alu_len);
  size_t nnn_begin = 0, nnn_end = 0;
  size_t ref_len_aluN = alu_len + NNN_LEN;

  vector <string>::iterator fi = fa_names.begin();
  map<string, TSeq> name_to_ref;
  seqan::BamStream bamStream(fin_sam.c_str());  
  seqan::BamAlignmentRecord record;
  while (!atEnd(bamStream)){
    readRecord(record, bamStream);
    fa_name = *fi++;
    bool alu_is_left = (fa_name[1] == 'l');
    nnn_begin = alu_is_left ? alu_len : (ref_end - ref_begin);
    nnn_end = nnn_begin + NNN_LEN;
    if (!findSplit(split_begin, jump_len, record)) continue;  // split_end = split_begin + jump_len
    if ( split_begin > nnn_begin || split_begin + jump_len < nnn_end) continue;

    if ( name_to_ref.find(fa_name) == name_to_ref.end() ) 
      read_fasta_by_name(name_to_ref[fa_name], fin_ref, fa_name);
    moveSplit_toRefBegin(split_begin, name_to_ref[fa_name], nnn_begin, jump_len);    
    
    if ( alu_is_left )  split_pos.push_back( new SPLIT_INFO(split_begin + jump_len - ref_len_aluN, true) ) ;
    else  split_pos.push_back( new SPLIT_INFO(split_begin, false) );     
  }
  seqan::close(bamStream);
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  int opt;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];

  string path1 = read_config(config_file, "file_alu_mate1") ;    
  boost::timer clocki;    
  clocki.restart();

  map<int, string> ID_pn;
  get_pn(read_config(config_file, "file_pn"), ID_pn);
  map<int, seqan::CharString> rID_chrn;

  if (opt == 1) {    
    string path_sam = path1 + "split_mapping_splazer/";  
    string path_map = path1 + "split_mapping_clip/";
    string path_ref = path1 + "split_mapping/";

    map <string, set< pair<int, int> > > pos_of_clippedReads;    
    ifstream fin(read_config(config_file, "file_pnIdx_used").c_str());
    int idx_pn;
    while (fin >> idx_pn) {
      string pn = ID_pn[idx_pn];
      get_filename_fastq(path_map + pn + ".map", pos_of_clippedReads); // only to update pos_of_clippedReads
    }
    fin.close();    
    map<string, set< pair<int, int> > >::iterator pi;
    string file_sam_prefix, fin_sam, fin_ref, fout_fa, fout_vcf; 
    vector< SPLIT_INFO *> split_pos;
    vector< SPLIT_INFO *>::iterator si;
    for ( pi = pos_of_clippedReads.begin(); pi != pos_of_clippedReads.end(); pi++ ) {
      string chrn = pi->first;
      for (set< pair<int, int> >::iterator ri = (pi->second).begin(); ri != (pi->second).end(); ri++) {	
	file_sam_prefix = path_sam + chrn + "_" + int_to_string( (*ri).first ) + "_" + int_to_string ((*ri).second);
	fin_ref = path_ref + chrn + "_" + int_to_string( (*ri).first ) + "_" + int_to_string ((*ri).second) + ".fa";
	// clear previous counts
	for ( si = split_pos.begin(); si != split_pos.end(); si++) delete *si;
	split_pos.clear();

	fin_sam = file_sam_prefix + "_f.sam";
	if (is_nonempty_file(fin_sam) )
	  read_sam(fin_sam, fin_ref, split_pos);
	fin_sam = file_sam_prefix + "_r.sam";
	if (is_nonempty_file(fin_sam) )
	  read_sam(fin_sam, fin_ref, split_pos);
	if (split_pos.empty()) continue;

	if (split_pos.size() > 2) {
	  cout << "done with " << file_sam_prefix  << endl;  
	  for ( si = split_pos.begin(); si != split_pos.end(); si++ )
	    cout << (*si)->split_to_refBegin << " " << (*si)->alu_is_left << endl;	
	  exit(0);
	}
	// a split position that is supported by at least 60% of the split reads.	
	if (split_pos.size() > 25) 
	  exit(0);
      }
    }
  }
  
}
