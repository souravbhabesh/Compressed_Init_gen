SRC = $(wildcard *.cpp)
OBJS = $(SRC:.cpp=.o)
DEPS = lattice_variables.h
CXX = g++
DEBUG = -g
CXXFLAGS = -Wall -c $(DEBUG) -std=c++11 -mcmodel=medium
LFLAGS = $(DEBUG) -O2 -Wall

$COMPRESS : $(OBJS)
	$(CXX) -o COMPRESS $(OBJS) $(LFLAGS)

lattice.o : $(DEPS) lattice.h  lattice.cpp 
	$(CXX) $(CXXFLAGS) lattice.cpp

streamCompaction.o : $(DEPS) streamCompaction.h  streamCompaction.cpp
	$(CXX) $(CXXFLAGS) streamCompaction.cpp

main.o : $(DEPS) lattice.h streamCompaction.h main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

clean:
	\rm *.o *~ COMPRESS

