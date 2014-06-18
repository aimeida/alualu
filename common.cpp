#include "common.h"

ConfigFileHandler::ConfigFileHandler(string config_file) {
  ifstream fin(config_file.c_str());
  if (!fin) 
    try {
      throw 1;    
    } catch(int e) {
      cerr << "#### ERROR #### file: "<< config_file << " not exists!!!" << endl;
    }   
  stringstream ss;
  string line, key, value;
  while (getline(fin, line)) {
    if (line[0] == '#')  continue;
    ss.clear(); ss.str( line );  
    ss >> key >> value;
    configs[key] = value;
  }
  fin.close();
}

string ConfigFileHandler::get_conf(string key) {
  map <string, string>::iterator ci;
  if ( (ci = configs.find(key)) == configs.end() ) {
    cerr << "#### ERROR #### key: " << key <<" doesn't exist\n";
    exit(1);
  }
  return ci->second;
}

string int_to_string(int i) {
  stringstream ss;
  ss << i;
  return ss.str();
}

string get_pn(string pn_file, int idx_pn){
  ifstream fin( pn_file.c_str());
  string pn;
  bool pn_exist = false;
  int i = 0;
  while (fin >> pn) {
    if (i == idx_pn ) {
      fin.close();  
      pn_exist = true;
      break;
    }
    i++;
  }
  assert(pn_exist);
  return pn;
}

void get_pn(string pn_file, map<int, string> &ID_pn){
  ID_pn.clear();
  ifstream fin( pn_file.c_str());
  assert(fin);
  string pn;
  int i = 0;
  while (fin >> pn)  ID_pn[i++] = pn;
  fin.close();
}

void read_file_pn_used(string fn, std::set <string> & pns_used) {
  ifstream fin(fn.c_str()) ;
  assert(fin);
  string pn;
  while (fin >> pn) pns_used.insert(pn);
  fin.close();
}

void read_file_pn_used(string fn, vector <string> & pns_used) {
  ifstream fin(fn.c_str()) ;
  assert(fin);
  string pn;
  while (fin >> pn) pns_used.push_back(pn);
  fin.close();
}

int is_nonempty_file(string fn){
  FILE * pFile = fopen(fn.c_str(), "r");
  if (pFile == NULL) return -1; // not exist
  else {
    fseek (pFile, 0, SEEK_END);   // non-portable
    long size = ftell(pFile);
    fclose (pFile);
    return size;  // empty OR >0
  }
}


void nonempty_files_ls(string path1, vector <string> fns){
  DIR *d;
  struct dirent *dir;
  d = opendir(path1.c_str());
  fns.clear();
  if (d) {
    while ( (dir = readdir(d))!=NULL )  {
      string fn = dir->d_name;
      if (is_nonempty_file(fn)) fns.push_back(fn);
    }
    closedir(d);
  }
}

