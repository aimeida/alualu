#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include "delete_utils.h"
#include "gtest/gtest.h"

int GetGp(string & config_file, string  line, int errCode = -1) {
  ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
  string file_pn_used = cf_fh.get_conf( "file_pn_used");
  std::set <string> pns_used;
  read_file_pn_used( file_pn_used, pns_used);
  string pn = *(pns_used.begin());  // use the first pn
  string pdf_param = cf_fh.get_conf("pdf_param"); // 100_1000_5   
  string file_dist_prefix = cf_fh.get_conf("file_dist_prefix");
  map <int, EmpiricalPdf *> pdf_rg;    
  read_pdf_pn(file_dist_prefix, pn, pdf_param, pdf_rg);
  float *gp = new float[3];
  parseline_ins(line, cout, pdf_rg, 3.5, 300, errCode, true, gp); 
  if ( gp[0] >= max(gp[1], gp[2]) ) return 2;
  else if ( gp[1] >= max(gp[0], gp[2]) ) return 1;
  return 0;
}

TEST(GenoCallTest, Input1) {
  string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
  EXPECT_EQ(1, GetGp(config_file, "chr1 188438210 188438213,188438207 3 50 0", 1));
}

// TEST(GenoCallTest, Input2) {
//   string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
//   EXPECT_EQ(2, GetGp(config_file, "chr1 188438210 188438213,188438207 10 2 0", 3));
// }
// 
// TEST(GenoCallTest, Input22) {
//   string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
//   EXPECT_EQ(1, GetGp(config_file, "chr1 188438210 188438213,188438207 10 2 0") );
// }

// TEST(GenoCallTest, Input3) {
//   string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
//   EXPECT_EQ(2, GetGp(config_file, "chr21 22739230 22739230,0 7 0 3 1:286 0:505 7:196"));
// }
// 
// TEST(GenoCallTest, Input4) {
//   string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
//   EXPECT_EQ(1, GetGp(config_file, "chr21 9721759 0,9721759 1 12 7 0:626 0:364 0:509 0:324 0:553 0:335 0:315"));
// }
// 
// TEST(GenoCallTest, Input5) {
//   string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
//   EXPECT_EQ(2, GetGp(config_file, "chr21 31706099 31706106,31706092 3 0 1 2:443"));
// }
// 
// TEST(GenoCallTest, Input6) {
//   string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
//   EXPECT_EQ(2, GetGp(config_file, "chr21 37926423 37926430,37926416 10 1 1 2:526"));
// }

