CC = g++

DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_EXE = ./bin

INC = -I./include -I./libs/vtd-xml/include -I./libs/vtd-xml/include/XPathParser

SRC = $(DIR_SRC)/genetree.cpp $(DIR_SRC)/genetree_node.cpp $(DIR_SRC)/indexed_genetree_node.cpp $(DIR_SRC)/genetree_index.cpp
OBJ := $(subst $(DIR_SRC), $(DIR_OBJ), $(patsubst %.cpp, %.o, $(SRC)))

VTD_PARSER_OBJ = ./libs/vtd-xml/obj/vtdparser.o

CFLAGS = -g -Wall -std=c++17 -Wno-unused-private-field
CFLAGS_DBG = -Wall -g -std=c++17 -Wno-unused-private-field

PYTHON_PATH = $(shell which python | rev | cut -d'/' -f2- | rev)/..
PYTHON_VERSION = $(shell python -c "import sys; print(str(sys.version_info[0]) + '.' + str(sys.version_info[1]))")
PYTHON_INC = -I$(PYTHON_PATH)/lib/python$(PYTHON_VERSION)/site-packages/numpy/core/include -I$(PYTHON_PATH)/include/python$(PYTHON_VERSION) -L $(PYTHON_PATH)/lib -rpath $(PYTHON_PATH)/lib -lpython$(PYTHON_VERSION)

MPL_INC = -I./include/matplotlibcpp.h

all: gtutils

gtutils: $(DIR_EXE)/gtutils

benchmark: $(DIR_EXE)/benchmark

# benchmark: OBJ $(DIR_EXE)/gtutils.cpp
# 	$(CC) $(CFLAGS) $(OBJ) -o $(DIR_EXE)/$@

$(DIR_EXE)/gtutils: $(OBJ) $(DIR_SRC)/gtutils.cpp $(VTD_PARSER_OBJ)
	$(CC) $(CFLAGS) $(INC) $^ -o $@

$(DIR_EXE)/benchmark: $(OBJ) $(DIR_SRC)/benchmark.cpp $(VTD_PARSER_OBJ)
	$(CC) $(CFLAGS) $(INC) $(MPL_INC) $(PYTHON_INC) $^ -o $@

$(DIR_OBJ)/%.o: $(DIR_SRC)/%.cpp
	${CC} $(CFLAGS) $(INC) -c $< -o $@

clean:
	-rm $(DIR_OBJ)/*.o