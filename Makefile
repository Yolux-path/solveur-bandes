# general definitions
CXX=g++
CXXFLAGS=-Wall -Wextra -std=c++20 -pedantic -O2
LDFLAGS=$(CXXFLAGS)
LDLIBS=

SRC       = mat_vec_io.cpp mat_vec_norm.cpp \
			linalg_band.cpp table.cpp
OBJ       = $(SRC:.cpp=.o)
HPP       = $(SRC:.cpp=.h) mat_vec.h

SRC_EXE   = mini_prog.cpp
OBJ_EXE   = $(SRC_EXE:.cpp=.o)
EXE       = $(SRC_EXE:.cpp=)

all: $(EXE)

# Make tp1
mini_prog: mini_prog.o $(OBJ)
	$(CXX) -o $@ $(LDFLAGS) $^ 
# Make projet

# Dépendances (manuelles)
mini_prog.o test_linalg.o $(OBJ): $(HPP)

test: test_linalg
	./test_linalg

test_linalg: test_linalg.o $(OBJ)
	$(CXX) -o $@ $(LDFLAGS) $^ 


# General rule to compile cpp files -- do not change anything below
%.o: %.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $< 

clean:
	rm -rf $(OBJ) $(OBJ_EXE) test_linalg.o test_linalg
mrproper: clean
	rm -rf $(EXE)
