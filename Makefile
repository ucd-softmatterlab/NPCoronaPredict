EXEC=UnitedAtom
SRC=src

CXX=g++
CXXFLAGS=-std=c++17 -O3 -lboost_system -lboost_filesystem -fopenmp -lpthread -ftree-vectorizer-verbose=2 -Wall -Wextra -Wdisabled-optimization

HEADERS := $(shell find $(SRC) -name "*.h")

$(EXEC) : $(SRC)/main.cpp $(HEADERS)
	$(CXX) -o $(EXEC) $(SRC)/main.cpp $(CXXFLAGS)

clean :
	rm -f $(EXEC) core
