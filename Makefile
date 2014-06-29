#-*-makefile-*-
##make OD=debug debug/alu_delete
##make OD=opt3 opt3/alu_delete
include Makefile.vars

INCLUDE = -I/home/qianyuxx/local/include/
WARN= -W -Wall

CPPFLAGS = $(WARN) $(FLAG) $(INCLUDE) 
CC = g++ -fno-merge-constants -fopenmp

LIB = -lz -lbz2 $(B_LIB)

ALU_DELETE_FILES_CPP = utils.cpp alu_delete.cpp common.cpp delete_utils.cpp
ALU_DELETE_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_DELETE_FILES_CPP))

ALU_INSERT_FILES_CPP = alu_insert.cpp common.cpp utils.cpp insert_utils.cpp delete_utils.cpp
ALU_INSERT_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_INSERT_FILES_CPP))

ALU_INSERT2_FILES_CPP = alu_insert2.cpp common.cpp utils.cpp insert_utils.cpp
ALU_INSERT2_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_INSERT2_FILES_CPP))

ALU_HG18_FILES_CPP = utils.cpp alu_hg18.cpp common.cpp 
ALU_HG18_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_HG18_FILES_CPP))

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

$(OD)/alu_insert2: $(OD) $(ALU_INSERT2_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_INSERT2_FILES_O)

$(OD)/alu_hg18: $(OD) $(ALU_HG18_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_HG18_FILES_O)

$(OD)/build_dist: $(OD) $(BUILD_DIST_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(BUILD_DIST_FILES_O)

test.dep:
	$(CC) -c -MM $(INCLUDE) $(TEST_FILES_CPP) > test.dep	

depend:
	make test.dep 

clean:
	/bin/rm -rf $(OD)/* $(addprefix $(OD)/,$(EXE))

##fixme: add all.dep