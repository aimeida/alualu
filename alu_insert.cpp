#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "diststat.h"
#include "utils.h"
#include <sys/time.h>
#include <boost/timer.hpp>

typedef std::map<std::string, std::ofstream*> Map;
#define ALU_MINLEN 15
#define MAPPED_ALU_LEN 100
#define FLAG_L2000 -1

void match_alupos(){
  /*
  map<int, AluRefPos *> rID_alurefpos;
  map<int, AluRefPos *>::iterator ra;
  ////for (i = 0; i < length( header.sequenceInfos); i++) rID_chrx[i] = header.sequenceInfos[i].i1;  
  for ( vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++) {
    assert(getIdByName(nameStore, *ci, rID, nameStoreCache));
    rID_chrx[rID] = *ci;
    rID_alurefpos[rID] = new AluRefPos(file_alupos_prefix + *ci, true);    
  }
  for (ra = rID_alurefpos.begin(); ra != rID_alurefpos.end(); ra++) delete ra->second;
  rID_alurefpos.clear();  
  */
  /*
    if (record.rID < record.rNextId) { // only check it once
      if ( ( (ra = rID_alurefpos.find(record.rNextId)) != rID_alurefpos.end() ) and (ra->second)->insideAlu(record.pNext, record.pNext + MAPPED_ALU_LEN, ALU_MINLEN) ) {
	fout1 << "dif_chr " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rNextId << " " << record.pNext << "1\n";
      } else if ( ( (ra = rID_alurefpos.find(record.rID)) != rID_alurefpos.end() ) and (ra->second)->insideAlu(record.beginPos, record.beginPos + getAlignmentLengthInRef(record), ALU_MINLEN) ) {
	fout1 << "dif_chr " << record.qName << " " << record.rNextId << " " << record.pNext  << " " << record.rID << " " << record.beginPos << "1\n";
      } 
    */
}

bool empty_map(map<seqan::CharString, pair<int,int> > &map_chr, ofstream & fout,int rID, seqan::CharString &chrn){
  if (!map_chr.size()) return false;
  cerr << chrn << " size " << map_chr.size() << endl;
  map<seqan::CharString, pair<int,int> >::iterator mc;
  for (mc = map_chr.begin(); mc!= map_chr.end(); mc++)
    fout << "left " << mc->first << " " << rID << " " << (mc->second).first << " na -1 " << (mc->second).second << endl;
  map_chr.clear();
  return true;
}

//search alu_mate by flag, write coordiate file
bool alu_mate_flag( string & bam_input, map<int, seqan::CharString> &rID_chrx, ofstream &fout1){
  seqan::Stream<seqan::Bgzf> inStream;
  assert (open(inStream, bam_input.c_str(), "r"));
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;
  assert(!readRecord(header, context, inStream, seqan::Bam()) );

  int rID= -1, ia = 0;
  map<seqan::CharString, pair<int,int> > same_chr_left; // value: (pos, count)
  map<seqan::CharString, pair<int,int> >::iterator sc;

  //bool hasAlignments = false;
  //jumpToRegion(inStream, hasAlignments, context, 10, 100000, 100000, baiIndex); // debug chr19
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( hasFlagQCNoPass(record) or hasFlagDuplicate(record) or hasFlagUnmapped(record) or hasFlagNextUnmapped(record) or (not hasFlagMultiple(record))) continue;
    //if (ic++ > 6000) break;
    //// this FLAG is broken ==> if (hasFlagSecondary(record)) writeRecord(bamStreamOut, record);
    ia++;
    if (ia % 1000000 == 0) cerr << ia << " " << same_chr_left.size() << endl;
    if (record.rID != record.rNextId) { // NB: redundant
      fout1 << "dif_chr " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rNextId << " " << record.pNext << " " << numOfBestHits(record) <<  endl;      
      continue;
    } else {
      if ( rID_chrx.find(record.rID) == rID_chrx.end() ) continue;
      if ( record.rID != rID) {
	empty_map(same_chr_left, fout1, rID, rID_chrx[rID]);
	cerr << "done with " << rID_chrx[rID] << " " << rID << endl;
	rID = record.rID;	
      }       
      int hits = numOfBestHits(record);      
      if (record.beginPos < record.pNext) {	
	if ( hits > 1) { 
	  same_chr_left[record.qName] = make_pair(record.beginPos, hits);
	} else if ( (abs(record.tLen) > 2000) or (hasFlagNextRC(record) == hasFlagRC(record))) { 
	  same_chr_left[record.qName] = make_pair(record.beginPos, FLAG_L2000); // non-redundant
	}	
      } else if (record.beginPos > record.pNext) {
	sc = same_chr_left.find(record.qName);
	if (sc == same_chr_left.end()) {
	  if (hits > 1) fout1 << "mul_map " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << hits << endl;
	} else {
	  if ( (sc->second).second == FLAG_L2000) {
	    fout1 << "i2000_RC " << record.qName << " " << record.rID << " " << record.pNext << " " << record.rID << " " << record.beginPos << " " << hits << endl;
	  } else if (hits == 1) { // at least one pair is unique mapped
	    fout1 << "mul_map " << record.qName << " " << record.rID << " " << record.beginPos << " " << record.rID << " " << record.pNext << " " << (sc->second).second << endl;
	  }
	  same_chr_left.erase(sc);
	}	      
      }
    }
  }  
  seqan::close(inStream);     
  return 0;
}

int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);

  int opt, idx_pn;
  seqan::lexicalCast2(opt, argv[1]);
  string config_file = argv[2];
  seqan::lexicalCast2(idx_pn, argv[3]);

//  unsigned coverage_max;
//  seqan::lexicalCast2(coverage_max, (read_config(config_file, "coverage_max")));

  string pn = get_pn(read_config(config_file, "file_pn"), idx_pn);
  cerr << "reading pn: " << idx_pn << " " << pn << "..................\n";
  string bam_input = read_config(config_file, "file_bam_prefix") + pn + ".bam";
  string bai_input = bam_input + ".bai";  
  string file_alu_mate_prefix = read_config(config_file, "file_alu_mate_prefix");  

  string file_alupos_prefix = read_config(config_file, "file_alupos_prefix"); 
  vector<string> chrns;
  for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
  chrns.push_back("chrX");
  chrns.push_back("chrY");
  map<int, seqan::CharString> rID_chrx;
  get_rID_chrx(bam_input, chrns, rID_chrx);
  
  if (opt == 1) {
    ofstream fout1( (file_alu_mate_prefix + pn).c_str() );
    fout1 << "flag qname this_chr this_pos bad_chr bad_pos num_hits\n"; // bad_*: potential, need check
    alu_mate_flag(bam_input, rID_chrx, fout1);
    fout1.close();
  } else if (opt == 2) {
    string file_alu_mate_filter = read_config(config_file, "file_alu_mate_filter");  
    string f_input = file_alu_mate_prefix + pn + ".alumate_flag";
    for (map<int, seqan::CharString>::iterator rc = rID_chrx.begin(); rc != rID_chrx.end(); rc++) {
      
    }
    
    //
    //check_type(f_input, "dif_chr");
  }

  return 0;  
}
