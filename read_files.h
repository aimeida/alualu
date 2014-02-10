// template used in header filer, to avoid linking error
//template <typename Type> Type read_config(string config_file, string key);

string read_config(string config_file, string key);

class AluRefPos
{
  int flanking;
  queue<int> beginP, endP;
 public:
  AluRefPos(string file_alupos, int i);
  int updatePos(int &beginPos, int &endPos);
  ~AluRefPos(void);
};

void printtest();
