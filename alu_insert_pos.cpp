// util functions to estimate insert position
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#include "common.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

#define CLIP_BP 10
#define DISCORDANT_LEN 1200
#define WINLEN_ONE_SIDE 500
#define DEFAULT_READ_LEN 120
typedef seqan::Dna5String TSeq;
#define DIF_CHR 1
#define LONG_LEN 2
#define SAME_RC 3
#define SOFT_CLIP 4

// Concatenates a genome infix with a contig separated by 70Ns.
void get_alu_consensus( TSeq &alu_seq, string aluFile, seqan::CharString contig_name){
  // fasta index for alu_consensus
  seqan::FaiIndex fai_alu;
  if (read(fai_alu, seqan::toCString(aluFile))) {
    build(fai_alu, seqan::toCString(aluFile));
    seqan::CharString aluFaiFile = aluFile;
    aluFaiFile += ".fai";
    write(fai_alu, toCString(aluFaiFile));
  }
  unsigned aluIndex = 0;  
  assert (getIdByName(fai_alu, contig_name, aluIndex));
  readSequence(alu_seq, fai_alu, aluIndex);
  //fa_seq += alu_seq;
}

int is_discordant(seqan::BamAlignmentRecord &record){
  if ( record.rID != record.rNextId) return DIF_CHR;
  if ( abs(record.tLen) > 2 * DISCORDANT_LEN)  return LONG_LEN;
  if ( hasFlagNextRC(record) == hasFlagRC(record)) return SAME_RC;
  if ( has_soft_last(record, CLIP_BP) or has_soft_first(record, CLIP_BP)) return SOFT_CLIP;
  return 0;
}

bool check_region_fastq(string bam_input, string bai_input, map<int, seqan::CharString> &rID_chrn, seqan::Stream<seqan::Bgzf> &inStream, seqan::SequenceStream &fastq_f, seqan::SequenceStream &fastq_r, seqan::BamIndex<seqan::Bai> &baiIndex,TBamIOContext &context, int rID, int beginPos, int endPos, seqan::BamHeader &header){  
  bool hasAlignments = false;

  seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
  bamStreamOut.header = header;
  
  if (!jumpToRegion(inStream, hasAlignments, context, rID, beginPos, endPos, baiIndex)) return 0;
  if (!hasAlignments) return 0;
  seqan::BamAlignmentRecord record, record2;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( record.rID != rID || record.beginPos >= endPos) break;
    if ( record.beginPos < beginPos) continue;            
    if ( hasFlagQCNoPass(record) or hasFlagDuplicate(record) or hasFlagUnmapped(record) or hasFlagNextUnmapped(record) or (not hasFlagMultiple(record))) continue;
    int read_type = is_discordant(record);
    if ( read_type == DIF_CHR or read_type == LONG_LEN) {  // ignore read_type == SAME_RC for now	
      /*
      string chrn2 = seqan::toCString(rID_chrn[record.rNextId]);
      string qname2 = seqan::toCString(record.qName);
      find_read(bam_input, bai_input, chrn2, qname2, record.pNext, record2, 0);
      //// don't reverse !!! cut sequence here !
      if ( hasFlagRC(record2) ) {    // fastq always forward
	reverseComplement(record2.seq);
	reverse(record2.qual);
      } 
      if ( hasFlagRC(record) ) writeRecord(fastq_f, record2.qName, record2.seq, record2.qual);
      else writeRecord(fastq_r, record2.qName, record2.seq, record2.qual);
      */
      continue;
    }
    if (read_type == SOFT_CLIP )  {   // useful !! fixme: infix --> cut off end of read ????
      cerr << "write "<< record.qName << endl;
      writeRecord(bamStreamOut, record);
      if ( hasFlagRC(record) ) { 
	reverseComplement(record.seq); reverse(record.qual);
	writeRecord(fastq_r, record.qName, record.seq, record.qual);
      } else { 
	writeRecord(fastq_f, record.qName, record.seq, record.qual);
      }
    }
  }
  seqan::close(bamStreamOut);
  return 1;
}

void call_splazers(string &cmd, string &bin_splazers, string param, string fn1, string fn2, string fn3) {
  cmd = bin_splazers + " " + param + " " + fn1 + " " + fn2 + " -o " + fn3;
}

void write_ref(string seq_ori, string fn, string fa_name, TSeq &nnn, TSeq &fa_seq, TSeq alu_seq, ofstream &fout) {
  TSeq fa_seq2;
  if (seq_ori == "pr") { // Positive chain, alu Right 
    fa_seq2 = fa_seq;
    fa_seq2 += nnn;
    fa_seq2 += alu_seq;
  } else if (seq_ori == "pl") { 
    fa_seq2 = alu_seq;
    fa_seq2 += nnn;
    fa_seq2 += fa_seq;
  } else if (seq_ori == "nr") {
    seqan::reverseComplement(alu_seq);
    fa_seq2 = fa_seq;
    fa_seq2 += nnn;
    fa_seq2 += alu_seq;
  } else if (seq_ori == "nl") {
    seqan::reverseComplement(alu_seq); 
    fa_seq2 = alu_seq;
    fa_seq2 += nnn;
    fa_seq2 += fa_seq;
  }
  writeRecord(fout, seq_ori + "_" + fa_name, fa_seq2, seqan::Fasta()); 
}


int main( int argc, char* argv[] )
{

  if (argc < 2) exit(1);
  int opt;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/config.properties";
  string pn_file = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/PN_all";
  string path_output = "/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/20_chr0/insert_alu1/split_mapping/";

  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  chrns.push_back("chrX");
  chrns.push_back("chrY");

    string chrn = "chr1";
    int region_begin = 100766793;
    int region_end = 100766944;
    int ref_begin = region_begin - WINLEN_ONE_SIDE;
    int ref_end = region_end + 2 * DEFAULT_READ_LEN;
    int check_begin = region_begin - WINLEN_ONE_SIDE;
    int check_end = region_end + WINLEN_ONE_SIDE;


  if ( opt == 1 ) {  // ref.fa with alu consensus inserted
    int idx_pn;  // start from 0
    seqan::lexicalCast2(idx_pn, argv[2]);
    string pn = get_pn(pn_file, idx_pn);    
    string bam_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/bamFile/" + pn + ".bam";
    string bai_input = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/bamFile/" + pn + ".bam.bai"; 
    seqan::BamAlignmentRecord record;  
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);
    seqan::BamHeader header;
    seqan::Stream<seqan::Bgzf> inStream;
    assert (open(inStream, bam_input.c_str(), "r"));
    seqan::BamIndex<seqan::Bai> baiIndex;
    assert (!read(baiIndex, bai_input.c_str()));
    assert(!readRecord(header, context, inStream, seqan::Bam()) );
    int rID = 0;
    assert (getIdByName(nameStore, chrn, rID, nameStoreCache));
    map<int, seqan::CharString> rID_chrn;
    get_rID_chrn(bam_input, chrns, rID_chrn);    
    
    string file_prefix = path_output + chrn + "_" + int_to_string(ref_begin);        
    seqan::SequenceStream fastq_f(seqan::toCString(file_prefix + "_f.fastq"), seqan::SequenceStream::WRITE); // to match forward
    seqan::SequenceStream fastq_r(seqan::toCString(file_prefix + "_r.fastq"), seqan::SequenceStream::WRITE); 
    check_region_fastq(bam_input, bai_input, rID_chrn, inStream, fastq_f, fastq_r, baiIndex, context, rID, check_begin, check_end, header);
    seqan::close(fastq_f);
    seqan::close(fastq_r);
    cerr << "output to " << file_prefix << endl;

  } else if (opt == 2) {  // call splazers for split mapping 
    // step 1, write ref fa file 
    string file_fa_prefix = read_config(config_file, "file_fa_prefix");
    seqan::FaiIndex faiIndex;
    unsigned fa_idx = 0;
    assert (!read(faiIndex, (file_fa_prefix + chrn + ".fa").c_str()) );
    assert (getIdByName(faiIndex, chrn, fa_idx));
    TSeq nnn;
    resize(nnn, 70, 'N');  
    TSeq fa_seq = (TSeq)fasta_seq(faiIndex, fa_idx, ref_begin, ref_end, true);

    //string aluFile = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/alu_jon.fa";
    string aluFile = "/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/alu_cons.fa";
    TSeq alu_seq;
    seqan::CharString contig_name = "AluY";
    get_alu_consensus(alu_seq, aluFile, contig_name);    
    stringstream ss;
    ss << contig_name << "_" << ref_end << "_" << length(alu_seq);
    
    // step 2, mapping
    string file_prefix = path_output + "chr1_100766293";
    //string file_prefix = argv[2];
    //string param1 = "-id -i 80 -m 1 -tr 95"; // fixme: some read length is 120 bp 
    //string param1 = "-id -i 80 -m 1 -tr 85"; 
    string param1 = "-id -i 80 -m 1";
    string param2 = "-sm 20 -ep 3 -es 3 -minG 5 -of 4"; 

    //  -id, allow indels.
    //  -i, --percent-identity 
    //  -a, dump alignment, -f only forward match 
    //  -m 1, only output best hits 

    string bin_splazers = read_config(config_file , "bin_splazers"); 
    string file_ref = file_prefix+"_ref.fa";
    ofstream fout(file_ref.c_str());
    char *seq_oris[] = {"pr", "pl", "nr", "nl"};
    for (int i = 0; i < 4; i ++ ) 
      write_ref(seq_oris[i], file_ref, ss.str(), nnn, fa_seq, alu_seq, fout);
    fout.close();

    string cmd;
    // split mapping
    call_splazers(cmd, bin_splazers, param1 + " -f " + param2, file_ref, file_prefix + "_f.fastq", file_prefix + "_f.sam");
    system( cmd.c_str() );
    call_splazers(cmd, bin_splazers, param1 + " -r " + param2, file_ref, file_prefix + "_r.fastq", file_prefix + "_r.sam");
    //cout << cmd << endl;
    system( cmd.c_str() );
    
    /* // full mapping 
    call_splazers( bin_splazers, param1 + " -f -sm 0 -of 1", file_prefix+"_ref.fa", file_prefix + "_f.fastq", file_prefix + "_f.fa");
    call_splazers( bin_splazers, param1 + " -r -sm 0 -of 1", file_prefix+"_ref.fa", file_prefix + "_r.fastq", file_prefix + "_r.fa");    

    // debug mode, -a print out alignment
    call_splazers( bin_splazers, param1 + " -a -f -sm 0 -of 3", file_prefix+"_p.ref", file_prefix + "_f.fastq", file_prefix + "_pf.gff");
    */

  }
}
