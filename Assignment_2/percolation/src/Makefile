CXXFLAGS = -Wall -Wextra -Wpedantic -O3 -std=c++11 -flto -pthread
LDLIBS = -lblas -llapack -larmadillo

# Reminder:
# $@ - target
# $^ - all dependencies
# $< - first dependency

all: percolation.prog

percolation.prog: main.o percolation.o utilities.o
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LDLIBS)

main.o: percolation.hpp utilities.hpp
utilities.o: utilities.hpp percolation.hpp
percolation.o: percolation.hpp

run: all
	time ./percolation.prog

clean:
	rm -f *.o
	rm -f *.prog
