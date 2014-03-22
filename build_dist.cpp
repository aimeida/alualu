// build distribution file based on flow cells (reading group)
//#define USING_MAIN_1
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"

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

void read_pn_chr(string &rg_lenCounts_file, string &bamInput_file, string &chrn){
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
  if(!getIdByName(nameStore, chrn, rID, nameStoreCache)) 
    cerr << "ERROR: Reference sequence named "<< chrn << " not known.\n";
  seqan::BamStream bamStreamIn(bamInput_file.c_str());
  map <seqan::CharString, map<int, int> > rg_lenCounts; // strata by reading group 
  size_t i = 0;
  unsigned idx_RG;
  size_t true_counts = 0, false_counts = 0;
  while (!atEnd(bamStreamIn)) {
    readRecord(record, bamStreamIn);
    if (record.rID != rID) {
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

void read_pn(string &rg_lenCounts_file, string &bamInput_file){
  seqan::BamAlignmentRecord record;  
  seqan::BamStream bamStreamIn(bamInput_file.c_str());
  map <seqan::CharString, map<int, int> > rg_lenCounts; // strata by reading group 
  unsigned idx_RG;
  while (!atEnd(bamStreamIn)) {
    readRecord(record, bamStreamIn);
    //if (i++ > 50000 ) break;    
    if ((not hasFlagQCNoPass(record) ) and hasFlagAllProper(record) and (not hasFlagDuplicate(record)) and hasFlagMultiple(record) ) {
      if (hasFlagFirst(record)) continue; // look at only one end of the pair      
      seqan::BamTagsDict tags(record.tags);
      if (!findTagKey(idx_RG, tags, "RG")) continue;
      seqan::CharString rg = getTagValue(tags, idx_RG);
      addKey(rg_lenCounts[rg], abs(record.tLen));    
    }    
  }  
  seqan::close(bamStreamIn);
  write_counts(rg_lenCounts, rg_lenCounts_file); 
}


void write_pdf(string &f_count, string &f_prob, int len_min, int len_max, int bin_width){
  int insertlen, count, counts = 0;
  map <int, int> bin_counts;
  ifstream fin(f_count.c_str());
  assert(fin);
  while (fin >> insertlen >> count) {
    if (insertlen < len_min) continue;
    if (insertlen >= len_max) break;    
    int id_bin = (insertlen - len_min)/bin_width;
    addKey(bin_counts, id_bin, count);
    counts += count;
  }
  fin.close();

  int half_bin_width = bin_width / 2;
  ofstream fout(f_prob.c_str());
  for (map<int, int>::iterator bc = bin_counts.begin(); bc != bin_counts.end(); bc++)
    fout << len_min + (bc->first) * bin_width + half_bin_width << " " << (bc->second) / (float)counts << endl;
  fout.close();
  cerr << "written into " << f_prob << endl;
}


int main( int argc, char* argv[] )
{
  if (argc < 2) exit(1);
  int opt;
  seqan::lexicalCast2(opt, argv[1]);
  string pn_file = argv[2];  
  int idx_pn;  // start from 0, run diff pn parallel 
  seqan::lexicalCast2(idx_pn, argv[3]);  
  string in_path = argv[4];
  string out_path = argv[5];
  string pn = get_pn(pn_file, idx_pn);
  
  if (opt == 1) {
    string chrn = argv[6]; 
    string rg_lenCounts_file, bamInput_file;
    cerr << "reading pn: " << pn << "..................\n";
    bamInput_file = in_path + pn + "/" + pn + ".bam";
    if (chrn != "chr0") {
      rg_lenCounts_file = out_path + pn + ".count." + chrn + "."; /* + RG*/
      read_pn_chr(rg_lenCounts_file, bamInput_file, chrn);
    } else {
      rg_lenCounts_file = out_path + pn + ".count."; /* + RG*/
      read_pn(rg_lenCounts_file, bamInput_file);
    }
  }  else if (opt == 2) {
    int len_min = 100, len_max = 1000, bin_width = 5; // set this as default    
    if (argc > 4) { 
      seqan::lexicalCast2(len_min, argv[4]);
      seqan::lexicalCast2(len_max, argv[5]);
      seqan::lexicalCast2(bin_width, argv[6]);
    }
    string rg;
    string suffix = int_to_string(len_min) + "_" + int_to_string(len_max) + "_" + int_to_string(bin_width); 
    ifstream fin( (out_path + "RG." + pn).c_str());
    assert(fin);
    while (fin >> rg) {
      string f_count = in_path + pn + ".count." + rg;
      string f_prob = out_path + pn + ".count." + rg + "." + suffix;
      cerr << "Reading " << f_count << "\nWriting " << f_prob << endl;
      write_pdf(f_count, f_prob, len_min, len_max, bin_width);    
    }
    fin.close();
  }
  
  return 0;
}

