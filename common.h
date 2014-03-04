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
#include <assert.h>
using namespace std;

typedef map<string, vector <int> >  ReadsPosStore;

string int_to_string(int i);
void print_vec(vector<int> &m);

template <class K, class V> 
void addKey(map <K,V> &m, K key, int cnt=1)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) == m.end()) {
    m[key] = cnt;
  } else (it->second) += cnt;
};

template <class K, class V>
void print_map(map <K,V> &m) 
{
  typename map <K,V>::iterator it;
  for (it = m.begin(); it != m.end(); it++) cerr << it->first << " " << it->second << endl;  
};

template <class K>
void print_vec(vector <K> &m) 
{
  typename vector <K>::iterator it;
  for (it = m.begin(); it != m.end(); it++) cerr << *it << " ";
  cerr << endl;
};


#endif /*COMMON_H*/
