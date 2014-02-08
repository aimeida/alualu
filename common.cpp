#include "common.h"

string int_to_string(int i) {
  stringstream ss;
  ss << i;
  return ss.str();
}
