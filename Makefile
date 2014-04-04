#-*-makefile-*-
OD=debug
## make debug/trybam
include Makefile.vars

INCLUDE = -I../libraries $(B_INC) 
WARN= -W -Wall

CPPFLAGS = $(WARN) $(FLAG) $(INCLUDE) 
CC = g++ -fno-merge-constants -fopenmp

LIB = -lz -lbz2 $(B_LIB)

BAM2FASTQ_FILES_CPP = bam2fastq.cpp
BAM2FASTQ_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(BAM2FASTQ_FILES_CPP))

ALU_DELETE_FILES_CPP = utils.cpp alu_delete.cpp common.cpp diststat.cpp
ALU_DELETE_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_DELETE_FILES_CPP))

ALU_INSERT_FILES_CPP = utils.cpp alu_insert.cpp common.cpp 
ALU_INSERT_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_INSERT_FILES_CPP))

INSERT_POS_FILES_CPP = utils.cpp insert_pos.cpp common.cpp 
INSERT_POS_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(INSERT_POS_FILES_CPP))

ALU_HG18_FILES_CPP = utils.cpp alu_hg18.cpp common.cpp 
ALU_HG18_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_HG18_FILES_CPP))

ALU_NOW_FILES_CPP = utils.cpp alu_now.cpp common.cpp diststat.cpp 
ALU_NOW_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_NOW_FILES_CPP)) 

UTILS_DEBUG_FILES_CPP = utils_debug.cpp common.cpp utils.cpp
UTILS_DEBUG_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(UTILS_DEBUG_FILES_CPP))

BUILD_DIST_FILES_CPP = build_dist.cpp common.cpp 
BUILD_DIST_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(BUILD_DIST_FILES_CPP))

FILES_CPP = *cpp *h 
FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(FILES_CPP))

TEST_FILES_CPP = alu_delete.cpp common.h diststat.h
TEST_FILES_O = $(patsubst %.cc,%.o,$(TEST_FILES_CPP))

EXE ?= alu_delete build_dist

$(OD) :
	mkdir $(OD)

$(OD)/%.o : %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

$(OD)/bam2fastq: $(OD) $(BAM2FASTQ_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(BAM2FASTQ_FILES_O)

$(OD)/alu_delete: $(OD) $(ALU_DELETE_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_DELETE_FILES_O)

$(OD)/alu_insert: $(OD) $(ALU_INSERT_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_INSERT_FILES_O)

$(OD)/insert_pos: $(OD) $(INSERT_POS_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(INSERT_POS_FILES_O)

$(OD)/alu_hg18: $(OD) $(ALU_HG18_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_HG18_FILES_O)

$(OD)/alu_now: $(OD) $(ALU_NOW_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_NOW_FILES_O)

$(OD)/build_dist: $(OD) $(BUILD_DIST_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(BUILD_DIST_FILES_O)

$(OD)/utils_debug: $(OD) $(UTILS_DEBUG_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(UTILS_DEBUG_FILES_O)


test.dep:
	$(CC) -c -MM $(INCLUDE) $(TEST_FILES_CPP) > test.dep	

depend:
	make test.dep 

clean:
	/bin/rm -rf $(OD)/*.o $(addprefix $(OD)/,$(EXE))

##fixme: add all.dep