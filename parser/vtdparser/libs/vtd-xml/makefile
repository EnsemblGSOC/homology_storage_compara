CC=g++

DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_EXE = ./bin

INC=-I./include -I./include/XPathParser

SRC += $(wildcard $(addsuffix /*.cpp, $(DIR_SRC)))
SRC += $(wildcard $(addsuffix /*.cc, $(DIR_SRC)))

OBJ := $(subst $(DIR_SRC), $(DIR_OBJ), $(patsubst %.cpp, %.o, $(SRC)))
OBJ := $(patsubst %.cc, %.o, $(OBJ))

VTD_PARSER_OBJ = vtdparser.o

EXE = benchmark_vtdxml
 
CFLAGS = -O3 -Wall -Winline -fomit-frame-pointer -fforce-addr -fexpensive-optimizations -fmerge-all-constants -fno-branch-count-reg -funroll-loops -falign-functions -falign-loops -falign-jumps -freorder-blocks -fprefetch-loop-arrays --param inline-unit-growth=250 --param max-inline-recursive-depth=2 --param large-function-growth=650 -fpermissive

CFLAGS2 = -c -g
LDFLAGS = -O3 -fomit-frame-pointer -fforce-addr -fexpensive-optimizations -lm 
LDFLAGS2  = -ggdb

all: $(DIR_OBJ)/$(VTD_PARSER_OBJ)

$(DIR_OBJ)/%.o: $(DIR_SRC)/%.cpp
	${CC} $(CFLAGS2) $(INC) -c $< -o $@

$(DIR_OBJ)/%.o: $(DIR_SRC)/%.cc
	${CC} $(CFLAGS2) $(INC) -c $< -o $@
	
#$(DIR_EXE)/$(EXE): $(DIR_OBJ)/$(VTD_PARSER_OBJ)
#	$(CC) $(CFLAGS) $(DIR_OBJ)/$(VTD_PARSER_OBJ) -o $@

$(DIR_OBJ)/$(VTD_PARSER_OBJ): $(OBJ)
	ld -r -o $@ $(OBJ)

clean:
	-rm $(DIR_OBJ)/*.o