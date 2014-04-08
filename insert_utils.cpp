// utils function for alu_insert 
#include "insert_utils.h"

bool readbam_clip(seqan::Stream<seqan::Bgzf> &inStream, seqan::BamIndex<seqan::Bai> &baiIndex, TBamIOContext &context, int rID, size_t ref_begin, int ref_end, seqan::BamStream &bamStreamOut, vector<seqan::CharString> & qnames, size_t offset){
  bool hasAlignments = false;
  if (!jumpToRegion(inStream, hasAlignments, context, rID, ref_begin, ref_end, baiIndex)) return false;
  if (!hasAlignments) return false;
  seqan::BamAlignmentRecord record;
  while (!atEnd(inStream)) {
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if ( record.rID != rID || record.beginPos >= ref_end) break;
    if ( record.beginPos + getAlignmentLengthInRef(record) <= ref_begin) continue; // has overlap with this region           
    if ( hasFlagQCNoPass(record) or hasFlagDuplicate(record) or hasFlagUnmapped(record) or hasFlagNextUnmapped(record) or (not hasFlagMultiple(record)) ) continue;
    if ( has_soft_last(record, offset) or has_soft_first(record, offset) ) {
      writeRecord(bamStreamOut, record);
      qnames.push_back(record.qName);
    }
  }
  return true;
}

