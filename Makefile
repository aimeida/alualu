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

ALU_DELETE_CPP = alu_delete.cpp
ALU_DELETE_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_DELETE_CPP))

FILES_CPP = *cpp *h 
FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(FILES_CPP))

TEST_FILES_CPP = alu_delete.cpp common.h configfile.h
TEST_FILES_O = $(patsubst %.cc,%.o,$(TEST_FILES_CPP))

EXE ?= alu_delete trybam 

$(OD) :
	mkdir $(OD)

$(OD)/%.o : %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

$(OD)/alu_delete: $(OD) $(ALU_DELETE_O) all.dep
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_DELETE_O)

$(OD)/trybam: $(OD) $(TRYBAM_FILES_O) all.dep
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(TRYBAM_FILES_O)

all.dep:
	$(CC) -c -MM $(INCLUDE) $(FILES_CPP) | sed 's/^\(.*\.o:\)/$$(OD)\/\1/g;' > all.dep

-include all.dep

test.dep:
	$(CC) -c -MM $(INCLUDE) $(TEST_FILES_CPP) > test.dep	

depend:
	rm all.dep
	make all.dep
	make test.dep 

clean:
	/bin/rm -rf $(OD)/*.o $(addprefix $(OD)/,$(EXE))

