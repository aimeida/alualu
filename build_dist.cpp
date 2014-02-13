// build distribution file based on flow cells (reading group)
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"

template <class K, class V> 
void addKey(map <K,V> &m, K key, int cnt=1)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) == m.end()) {
    m[key] = cnt;
  } else (it->second) += cnt;
}

void write_counts(map <seqan::CharString, map<int, int> > &rg_lenCounts, string &rg_lenCounts_file){
  map <seqan::CharString, map<int, int> >::iterator sii;
  map<int, int>::iterator si;
  cerr << "output to: " << rg_lenCounts_file << endl;
  for (sii = rg_lenCounts.begin(); sii != rg_lenCounts.end() ; sii++) {
    ofstream fout((rg_lenCounts_file + toCString(sii->first)).c_str());
    for (si = (sii->second).begin(); si != (sii->second).end(); si++) 
      fout << si->first << " "  << si->second << endl;
    fout.close();
  }
}

void read_pn(string &rg_lenCounts_file, string &bamInput_file, string &chrx){
  typedef seqan::StringSet<seqan::CharString> TNameStore;
  typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
  typedef seqan::BamIOContext<TNameStore>     TBamIOContext;  
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);
  seqan::BamHeader header;
  seqan::BamAlignmentRecord record;  
  seqan::Stream<seqan::Bgzf> inStream;
  open(inStream, bamInput_file.c_str(), "r");
  readRecord(header, context, inStream, seqan::Bam());
  seqan::close(inStream);
  int rID = 0;
  if(!getIdByName(nameStore, chrx, rID, nameStoreCache)) {
    cerr << "ERROR: Reference sequence named "<< chrx << " not known.\n";
  }
  //cerr << "rID " << rID << endl;
  
  seqan::BamStream bamStreamIn(bamInput_file.c_str());
  map <seqan::CharString, map<int, int> > rg_lenCounts; // strata by reading group 
  size_t i = 0;
  unsigned idx_RG;
  size_t true_counts = 0, false_counts = 0;
  while (!atEnd(bamStreamIn)) {
    readRecord(record, bamStreamIn);
    if (record.rID != rID) {
      //      continue;
      if (!i) continue;
      else break;
    }
    i++;
    //////if (i++ > 50000 ) break;    
    if ((not hasFlagQCNoPass(record) ) and hasFlagAllProper(record) and (not hasFlagDuplicate(record)) and hasFlagMultiple(record) ) {
      if (hasFlagFirst(record)) continue; // look at only one end of the pair      
      //cerr << record.beginPos << " " << record.pNext << " " << length(record.seq) << " " <<  getAlignmentLengthInRef(record) << endl; 
      seqan::BamTagsDict tags(record.tags);
      if (!findTagKey(idx_RG, tags, "RG")) continue;
      seqan::CharString rg = getTagValue(tags, idx_RG);
      addKey(rg_lenCounts[rg], abs(record.tLen));    
    }    
  }  
  seqan::close(bamStreamIn);
  if (true_counts or false_counts) cerr << "true counts: " << true_counts << ", false counts: " << false_counts << endl;
  write_counts(rg_lenCounts, rg_lenCounts_file); 
}


int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  string pn_file = argv[1];
  string bam_in_path = argv[2];
  string out_path = argv[3];
  string chrx = argv[4]; 
  int idx_pn;
  seqan::lexicalCast2(idx_pn, argv[5]);
  
  ifstream fin( pn_file.c_str());
  string pn;
  int i = 0;
  while (fin >> pn) {
    if (i == idx_pn ) {
      cerr << "reading pn: " << pn << "..................\n";
      string rg_lenCounts_file = out_path + pn + ".count." + chrx + "."; /* + RG*/
      string bamInput_file = bam_in_path + pn + "/" + pn + ".bam";
      read_pn(rg_lenCounts_file, bamInput_file, chrx);
    }
    i++;
  }
  fin.close();  
}
