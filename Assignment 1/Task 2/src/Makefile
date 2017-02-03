CXX = g++
FLAGS = -Wall -O2 --std=c++11
LIBS = -lblas -llapack -larmadillo
DEBUG_FLAGS = -ggdb

all: simulator.prog

run: all
	./simulator.prog

simulator.prog: main.o diffusion.o
	${CXX} $^ -o $@

main.o: main.cpp
	${CXX} $^ -c ${FLAGS} ${LIBS}

diffusion.o: diffusion.cpp
	${CXX} $^ -c ${FLAGS} ${LIBS}

debug:
	${CXX} $^ -o $@ ${FLAGS} ${LIBS} ${DEBUG_FLAGS}

clean:
	rm -f *.o
