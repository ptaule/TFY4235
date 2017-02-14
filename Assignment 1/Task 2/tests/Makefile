CXX = g++
FLAGS = -Wall -O3 --std=c++11
LIBS = -lblas -llapack -larmadillo
DEBUG_FLAGS = -ggdb

# Reminder:
# $@ - target
# $^ - all dependencies
# $< - first dependency

all: tests.prog

run: tests.prog
	time ./tests.prog

tests.prog: main.cpp tests.o ../src/diffusion.o
	${CXX} $^ -o $@ ${LIBS}

tests.o: tests.cpp ../src/diffusion.hpp
	${CXX} $< -c ${FLAGS}

../src/diffusion.o: ../src/diffusion.cpp
	cd ../src/ && ${CXX} $^ -c ${FLAGS} ${LIBS}

clean:
	rm -f *.o
