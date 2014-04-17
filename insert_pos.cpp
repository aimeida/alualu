#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

bool findSplit(int & beginPos,int &jump_len, seqan::BamAlignmentRecord &record) {
  char operation;
  int count;
  if ( length(record.cigar) < 3) return false;
  beginPos = record.beginPos;
  jump_len = 0;
  for (size_t i = 0; i < length(record.cigar); i++) {
    operation = record.cigar[i].operation;
    count = record.cigar[i].count;
    if ( operation == 'M' or operation == 'D') //  ??? or operation == 'S') 
      beginPos += count;
    if ( operation == 'N') {
      jump_len = count;
      break;
    }
  }
  return jump_len;
}

// move split pos further to ref_begin
void moveSplit_toRefBegin(int  & split_begin, TSeq const & ref, int nnn_begin, int jump_len){  
  while (split_begin < nnn_begin and split_begin + jump_len < (int)length(ref) and ref[split_begin] == ref[ split_begin + jump_len] ) 
    split_begin++;
} 

// write all reads that is useful for split mapping or inferring insert sequence !
bool readbam_loci(string chrn, seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, size_t region_begin, int region_end, ofstream &frecord) {
  
  int ref_begin, ref_end; // range for broken reads
  refPos_by_region(region_begin, region_end, ref_begin, ref_end);  
  bool hasAlignments = false;  
  if (!jumpToRegion(inStream, hasAlignments, context, rID, ref_begin, ref_end, baiIndex)) return false;
  if (!hasAlignments) return false;
  seqan::BamAlignmentRecord record;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( record.rID != rID || record.beginPos >= ref_end) break;
    int read_end = record.beginPos + getAlignmentLengthInRef(record);
    if ( read_end <= ref_begin) continue; 
    if ( ! QC_insert_read(record) ) continue;
    if ( has_soft_last(record, CLIP_BP) or has_soft_first(record, CLIP_BP) ) 
      frecord << chrn << " " << region_begin << " " << region_end << " " << record.beginPos << " " << read_end 
	      << " " << get_cigar(record) << " " << hasFlagRC(record) << " " << record.seq << endl;
  } 
  return true;
}

void reads_insert_loci(int idx_pn, vector<string> & chrns, string bam_input, string bai_input, string fin_pos, string fout_reads_fa){
  seqan::BamIndex<seqan::Bai> baiIndex;
  assert (!read(baiIndex, bai_input.c_str()));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  seqan::Stream<seqan::Bgzf> inStream;
  assert(open(inStream, bam_input.c_str(), "r"));
  assert(!readRecord(header, context, inStream, seqan::Bam()));    
  map<int, seqan::CharString> rID_chrn;
  get_rID_chrn(bam_input, chrns, rID_chrn);

  ofstream frecord(fout_reads_fa.c_str());
  frecord << "chrn region_begin region_end beginPos endPos cigar hasFlagRC seq\n" ;
  int region_begin, region_end;
  string line;
  for (map<int, seqan::CharString>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) {
    string chrn = toCString(rc->second);
    ifstream fin( (fin_pos + chrn).c_str());
    if (!fin) continue;  // ignore random chr
    while (getline(fin, line)) {
      if ( coverage_idx_pass(line, region_begin, region_end, idx_pn) )
	readbam_loci(chrn, inStream, baiIndex, context, rc->first, region_begin, region_end, frecord);
    }
    fin.close();
  }
  seqan::close(inStream);     
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  int opt;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];

  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  string path1 = read_config(config_file, "file_alu_insert1") ;    

  boost::timer clocki;    
  clocki.restart();  
  map<int, string> ID_pn;
  get_pn(read_config(config_file, "file_pn"), ID_pn);  
  
  if (opt == 1) { 
    int idx_pn;
    seqan::lexicalCast2(idx_pn, argv[3]);
    string fin_pos = path1 + "insert_pos.";
    string fout_reads = read_config(config_file, "file_insert_reads");
    if ( access( fout_reads.c_str(), 0 ) != 0 ) system( ("mkdir " + fout_reads).c_str() );          
    string pn = get_pn(read_config(config_file, "file_pn"), idx_pn);
    string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
    string bai_input = bam_input + ".bai";  
    string fout_reads_fa = fout_reads + pn + ".clip" ; 
    // takes 2 min
    reads_insert_loci(idx_pn, chrns, bam_input, bai_input, fin_pos, fout_reads_fa);
    //writeRecord(fout, record.qName, record.seq, record.qual, seqan::Fastq());
    //writeRecord(fout, record.qName, record.seq, seqan::Fasta());

  }
  return 0;
}
