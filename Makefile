#-*-makefile-*-
OD=debug
## make debug/trybam
include Makefile.vars

INCLUDE = -I../libraries $(B_INC) 
WARN= -W -Wall

CPPFLAGS = $(WARN) $(FLAG) $(INCLUDE) 
CC = g++ -fno-merge-constants -fopenmp

LIB = -lz -lbz2 $(B_LIB)

TRYBAM_FILES_CPP = trybam.cpp
TRYBAM_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(TRYBAM_FILES_CPP))

ALU_DELETE_FILES_CPP = utils.cpp alu_delete.cpp common.cpp diststat.cpp
ALU_DELETE_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_DELETE_FILES_CPP))

ALU_MULTI_FILES_CPP = alu_multi.cpp common.cpp utils.cpp
ALU_MULTI_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_MULTI_FILES_CPP)) 

ALU_NOW_FILES_CPP = alu_now.cpp common.cpp utils.cpp
ALU_NOW_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_NOW_FILES_CPP)) 

BUILD_DIST_FILES_CPP = build_dist.cpp common.cpp 
BUILD_DIST_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(BUILD_DIST_FILES_CPP))

FILES_CPP = *cpp *h 
FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(FILES_CPP))

TEST_FILES_CPP = alu_delete.cpp common.h diststat.h
TEST_FILES_O = $(patsubst %.cc,%.o,$(TEST_FILES_CPP))

EXE ?= alu_delete trybam build_dist

$(OD) :
	mkdir $(OD)

$(OD)/%.o : %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

$(OD)/alu_delete: $(OD) $(ALU_DELETE_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_DELETE_FILES_O)

$(OD)/alu_multi: $(OD) $(ALU_MULTI_FILES_O) 	
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_MULTI_FILES_O)	

$(OD)/alu_now: $(OD) $(ALU_NOW_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_NOW_FILES_O)

$(OD)/trybam: $(OD) $(TRYBAM_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(TRYBAM_FILES_O)

$(OD)/build_dist: $(OD) $(BUILD_DIST_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(BUILD_DIST_FILES_O)

test.dep:
	$(CC) -c -MM $(INCLUDE) $(TEST_FILES_CPP) > test.dep	

depend:
	make test.dep 

clean:
	/bin/rm -rf $(OD)/*.o $(addprefix $(OD)/,$(EXE))

##fixme: add all.dep