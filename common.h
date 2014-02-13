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
#include <math.h> 
using namespace std;

typedef map<string, vector <int> >  ReadsPosStore;
//typedef map<string, int>  ReadsPosStore;

string int_to_string(int i);
void print_vec(vector<int> &m);

#endif /*COMMON_H*/
