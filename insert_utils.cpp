// utils function for alu_insert 
#include "insert_utils.h"

bool alu_mate_flag( BamFileHandler * bam_fh, MapFO &fileMap ){
  seqan::BamAlignmentRecord record;
  int rID_pre = -1, ia = 0;
  while ( bam_fh -> fetch_a_read(record) ) { 
    if ( !QC_insert_read(record) ) continue;
    if ( bam_fh->rID_chrn.find(record.rID) == bam_fh->rID_chrn.end() ) continue;
    if ( record.rID != rID_pre) {
      if (rID_pre > -1) {
	cerr << "done with " << rID_pre << endl; 
      }
      rID_pre = record.rID;	
    }       
    //if (ia++ % 1000000 == 0) cerr << ia << " " << endl;
    int len_seq = length(record.seq);
    if ( getAlignmentLengthInRef(record) < (size_t)(len_seq - 10)) continue; // want very good quality reads 
    if (record.rID < record.rNextId and bam_fh->rID_chrn.find(record.rNextId) != bam_fh->rID_chrn.end() ) { 
      *(fileMap[record.rNextId]) << "dif_chr " << record.qName << " " << record.rNextId << " " << record.pNext << " " 
				 << record.rID << " " << record.beginPos << " " << len_seq << endl;      
      *(fileMap[record.rID]) << "dif_chr " << record.qName << " " << record.rID << " " << record.beginPos << " " 
				 << record.rNextId << " " << record.pNext << " " << len_seq << endl;            
      continue;
    } 

    if (record.rID == record.rNextId and record.beginPos < record.pNext) { 
      if (abs(record.tLen) > DISCORDANT_LEN) {
	*(fileMap[record.rID]) << "ilong " << record.qName << " " << record.rID << " " << record.pNext << " " 
			       << record.rID << " " << record.beginPos << " " << len_seq << endl;      
	*(fileMap[record.rID]) << "ilong " << record.qName << " " << record.rID << " " << record.beginPos << " " 
			       << record.rID << " " << record.pNext << " " << len_seq << endl;      
      }    
    }
  }
  return 0;
}

// only support 4 types for now 
string parse_alu_type(string alu_name){
  assert ( !alu_name.empty() );
  if ( alu_name.substr(0,4) == "AluY") return "AluY";
  if ( alu_name.substr(0,4) == "AluS") return "AluSx";
  if ( alu_name.substr(0,4) == "AluJ") {
    if ( alu_name.substr(0,5) == "AluJo") return "AluJo";
    if ( alu_name.substr(0,5) == "AluJb") return "AluJb";
    return "AluJo";
  }
  return "AluY";
}

