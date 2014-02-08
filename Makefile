#-*-makefile-*-
OD=debug
## make debug/trybam
include Makefile.vars

INCLUDE = -I../libraries
WARN= -W -Wall

CPPFLAGS = $(WARN) $(FLAG) $(INCLUDE) 
CC = g++ -fno-merge-constants -fopenmp

LIB = -lz -lbz2

TRYBAM_FILES_CPP = trybam.cpp
TRYBAM_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(TRYBAM_FILES_CPP))

ALU_DELETE_CPP = read_files.cpp alu_delete.cpp common.cpp
ALU_DELETE_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_DELETE_CPP))

FILES_CPP = *cpp *h 
FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(FILES_CPP))

TEST_FILES_CPP = alu_delete.cpp common.h diststat.h
TEST_FILES_O = $(patsubst %.cc,%.o,$(TEST_FILES_CPP))

EXE ?= alu_delete trybam 

$(OD) :
	mkdir $(OD)

$(OD)/%.o : %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

$(OD)/alu_delete: $(OD) $(ALU_DELETE_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_DELETE_O)

$(OD)/trybam: $(OD) $(TRYBAM_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(TRYBAM_FILES_O)


test.dep:
	$(CC) -c -MM $(INCLUDE) $(TEST_FILES_CPP) > test.dep	

depend:
	make test.dep 

clean:
	/bin/rm -rf $(OD)/*.o $(addprefix $(OD)/,$(EXE))

##fixme: add all.dep