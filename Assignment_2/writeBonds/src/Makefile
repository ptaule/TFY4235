CXX = g++
CXXFLAGS = -Wall -O3 --std=c++11
LDLIBS = -lblas -llapack -larmadillo

# Reminder:
# $@ - target
# $^ - all dependencies
# $< - first dependency

all: bonds.prog

run: all
	time ./bonds.prog

bonds.prog: main.o writeBonds.o
	$(CXX) -o $@ $^ $(LDLIBS)

main.o: main.cpp writeBonds.hpp
	$(CXX) -c $(CXXFLAGS) $<

writeBonds.o: writeBonds.cpp writeBonds.hpp
	$(CXX) -c $(CXXFLAGS) $<

clean:
	rm -f *.o
	rm -f *.prog
