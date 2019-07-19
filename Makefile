EXEC=UnitedAtom
SRC=src

CXX=g++
CXXFLAGS=-std=c++17 -O3 -lboost_system -lboost_filesystem -fopenmp -lpthread -march=native -mtune=native -Wall -Wextra -Wdisabled-optimization

HEADERS := $(shell find $(SRC) -name "*.h")


all : $(EXEC) pdb2cg1

$(EXEC) : $(SRC)/main.cpp $(HEADERS)
	$(CXX) -o $(EXEC) $(SRC)/main.cpp $(CXXFLAGS)


pdb2cg1 : $(SRC)/pdb2cg1.cpp
	$(CXX) -o pdb2cg1 $(SRC)/pdb2cg1.cpp $(CXXFLAGS)


clean :
	rm -f $(EXEC) core
