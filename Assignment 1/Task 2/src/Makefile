CXX = g++
FLAGS = -Wall -Wextra -Wpedantic -O3 --std=c++11
LIBS = -lblas -llapack -larmadillo
DEBUG_FLAGS = -ggdb

# Reminder:
# $@ - target
# $^ - all dependencies
# $< - first dependency

plot: run
	gnuplot diffusion.gpi
	# epstopdf ../res/diffusion.eps

run: all
	time ./diffusion.prog > diffusion.log

all: diffusion.prog

diffusion.prog: main.o diffusion.o
	${CXX} $^ -o $@ ${LIBS}

main.o: main.cpp diffusion.hpp
	${CXX} $< -c ${FLAGS}

diffusion.o: diffusion.cpp diffusion.hpp
	${CXX} $< -c ${FLAGS}

debug: main.cpp diffusion.cpp
	${CXX} $^ -o debug.prog ${FLAGS} ${LIBS} ${DEBUG_FLAGS}

clean:
	rm -f *.o
