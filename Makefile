#-*-makefile-*-
#OD=debug
#make debug/alu_delete
OD=opt3
include Makefile.vars
ARCH	= $(shell uname -i)
BOOST   = /nfs_mount/bioinfo/apps-$(ARCH)/boost/1.48.0
B_INC	= -I$(BOOST)/include
B_LIB	= -L$(BOOST)/lib -Wl,-rpath,$(BOOST)/lib

HDF5=/nfs_mount/bioinfo/apps-x86_64/hdf5-dg/BUILD_3_51
HDF5_include=$(HDF5)/include
HDF5_lib=$(HDF5)/lib
HDF5_link = -Wl,-rpath,$(HDF5_lib)

INCLUDE = -I/home/qianyuxx/local/include/ -I/katla/groups/statistics/bjarnih/code/seqan-trunk/core/include -I/katla/groups/statistics/bjarnih/code/seqan-trunk/extras/include -I/home/stat/include -I/home/stat/lib64/R/include -Il/katla/groups/statistics/bjarnih/code/stat/libraries   -I/home/stat/include/cppunit -I$(HDF5_include) $(B_INC)
WARN= -W -Wall

CPPFLAGS = $(WARN) $(FLAG) $(INCLUDE) 
CC = g++ -fno-merge-constants -fopenmp

#LIB = -lz -lbz2 $(B_LIB)
LIB = -lz $(B_LIB)

ALU_DELETE_FILES_CPP = utils.cpp alu_delete.cpp common.cpp delete_utils.cpp
ALU_DELETE_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_DELETE_FILES_CPP))

ALU_INSERT_FILES_CPP = alu_insert.cpp common.cpp utils.cpp insert_utils.cpp delete_utils.cpp
ALU_INSERT_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_INSERT_FILES_CPP))

BUILD_DIST_FILES_CPP = build_dist.cpp common.cpp utils.cpp
BUILD_DIST_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(BUILD_DIST_FILES_CPP))

FILES_CPP = *cpp *h 
FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(FILES_CPP))

TEST_FILES_CPP = alu_delete.cpp common.h 
TEST_FILES_O = $(patsubst %.cc,%.o,$(TEST_FILES_CPP))

EXE ?= alu_delete build_dist

$(OD) :
	mkdir $(OD)

$(OD)/%.o : %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

$(OD)/alu_delete: $(OD) $(ALU_DELETE_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_DELETE_FILES_O)

$(OD)/alu_insert: $(OD) $(ALU_INSERT_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_INSERT_FILES_O)

$(OD)/build_dist: $(OD) $(BUILD_DIST_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(BUILD_DIST_FILES_O)

test.dep:
	$(CC) -c -MM $(INCLUDE) $(TEST_FILES_CPP) > test.dep	

depend:
	make test.dep 

clean:
	/bin/rm -rf $(OD)/* $(addprefix $(OD)/,$(EXE))

##fixme: add all.dep