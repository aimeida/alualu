// for header files and debug print 
#ifndef COMMON_H
#define COMMON_H

#include <map>
#include <set>
#include <queue>
#include <vector>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <math.h> 
#include <stdio.h>
#include <assert.h>
using namespace std;

typedef map<string, vector <int> >  ReadsPosStore;

string int_to_string(int i);
void print_vec(vector<int> &m);
string get_pn(string pn_file, int idx_pn);
void get_pn(string pn_file, map<int, string> &ID_pn);
int is_nonempty_file(string fn);

template <class V>
void clear_vec_of_ptr( vector < V *> & m)
{
  typename vector < V *>::iterator it;
  for (it = m.begin(); it != m.end(); it++) delete *it;
  m.clear();
};


template <class K, class V> 
void addKey(map <K,V> &m, K key, int cnt=1)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) == m.end()) {
    m[key] = cnt;
  } else (it->second) += cnt;
};

template <class K, class V>
void print_map(map <K,V> &m, size_t n_max = 0) {
  cout << "size of map " << m.size() << endl;
  typename map <K,V>::iterator it;
  if (n_max == 0) n_max = m.size();
  size_t ni = 0;
  for (it = m.begin(); it != m.end(),ni < n_max; it++, ni++) cerr << it->first << ":" << it->second << " ";
  cerr << endl;  
};

template <class K>
void print_vec(vector <K> &m) 
{
  typename vector <K>::iterator it;
  for (it = m.begin(); it != m.end(); it++) cerr << *it << " ";
  cerr << endl;
};


#endif /*COMMON_H*/
