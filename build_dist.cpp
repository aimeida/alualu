// build distribution file based on flow cells (reading group)
//#define USING_MAIN_1
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "common.h"
#include "utils.h"

void write_counts(map <seqan::CharString, map<int, int> > &rg_lenCounts, string &rg_lenCnt_file){
  map <seqan::CharString, map<int, int> >::iterator sii;
  map<int, int>::iterator si;
  cerr << "output to: " << rg_lenCnt_file << endl;
  for (sii = rg_lenCounts.begin(); sii != rg_lenCounts.end() ; sii++) {
    ofstream fout((rg_lenCnt_file + toCString(sii->first)).c_str());
    for (si = (sii->second).begin(); si != (sii->second).end(); si++) 
      fout << si->first << " "  << si->second << endl;
    fout.close();
  }
}

void read_pn_chr(string &rg_lenCnt_file, string &bamInput_file, string &chrn){
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

  write_counts(rg_lenCounts, rg_lenCnt_file); 
}

void read_pn(map <seqan::CharString, map<int, int> > &rg_lenCnt, string &bamInput_file, int jump_first, int max_read_num){
  seqan::BamAlignmentRecord record;  
  seqan::BamStream bamStreamIn(bamInput_file.c_str());
  unsigned idx_RG;
  int i = 0;
  while (!atEnd(bamStreamIn)) {
    i++;

    readRecord(record, bamStreamIn);
    if ( i < jump_first) continue;
    if ( max_read_num and i > max_read_num ) break;    
    
    if ((not hasFlagQCNoPass(record) ) and hasFlagAllProper(record) and (not hasFlagDuplicate(record)) and hasFlagMultiple(record) ) {
      if (hasFlagFirst(record)) continue; // look at only one end of the pair      
      seqan::BamTagsDict tags(record.tags);
      if (!findTagKey(idx_RG, tags, "RG")) continue;
      seqan::CharString rg = getTagValue(tags, idx_RG);
      addKey(rg_lenCnt[rg], abs(record.tLen));    
    }    
  }  
  seqan::close(bamStreamIn);
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
  string config_file = argv[1];
  int idx_pn;  // start from 0, run diff pn parallel 
  seqan::lexicalCast2(idx_pn, argv[2]);  
  string chrn = argv[3]; 

  map<int, string> ID_pn;
  get_pn(read_config(config_file, "file_pn"), ID_pn);
  string pn = ID_pn[idx_pn];
  
  string path_count = read_config(config_file,"file_insertlen_count");
  string path_prob = read_config(config_file,"file_dist_prefix");
  string bamInput_file = read_config(config_file, "file_bam_prefix") + pn + ".bam";  
  string rg_lenCnt_file = path_count + pn + ".count."; /* + RG*/

  //cout << "counting: " << pn << " now\n";
  if (chrn != "chr0") {
    rg_lenCnt_file += chrn + "."; /* + RG*/
    read_pn_chr(rg_lenCnt_file, bamInput_file, chrn);
  } else {
    map <seqan::CharString, map<int, int> > rg_lenCnt; // strata by reading group 
    read_pn(rg_lenCnt, bamInput_file, 5e3, 5e7);  // use only 5% reads to estimate

    cout << rg_lenCnt_file << endl;
    write_counts(rg_lenCnt, rg_lenCnt_file); 
    // write down rg groups 
    ofstream fout( (path_prob + "RG." + pn).c_str()); 
    for (map <seqan::CharString, map<int, int> >::iterator ri = rg_lenCnt.begin(); ri != rg_lenCnt.end(); ri++ )
      fout << ri->first << endl; 
    fout.close(); 
  }
  
  cerr << "counting: " << pn << " done\n";
  stringstream ss;
  string pdf_param = read_config(config_file, "pdf_param");
  ss.str(pdf_param);
  ///int len_min = 100, len_max = 1000, bin_width = 5; 
  int len_min, len_max, bin_width; 
  string token;
  getline(ss, token, '_');
  seqan::lexicalCast2(len_min, token);
  getline(ss, token, '_');
  seqan::lexicalCast2(len_max, token);
  getline(ss, token, '_');
  seqan::lexicalCast2(bin_width, token);
  string rg;
  ifstream fin( (path_prob + "RG." + pn).c_str());
  assert(fin);
  while (fin >> rg) {
    string f_count = path_count + pn + ".count." + rg;
    string f_prob = path_prob + pn + ".count." + rg + "." + pdf_param;
    cerr << "Reading " << f_count << "\nWriting " << f_prob << endl;
    write_pdf(f_count, f_prob, len_min, len_max, bin_width);    
    }
  fin.close();
  
  return 0;
}

