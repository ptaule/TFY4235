CXX = g++
FLAGS = -Wall -O2 --std=c++11
LIBS = -lblas -llapack -larmadillo
DEBUG_FLAGS = -ggdb

all: randomWalk

randomWalk: randomWalk.cpp
	${CXX} $^ -o $@ ${FLAGS} ${LIBS}

debug:
	${CXX} $^ -o $@ ${FLAGS} ${LIBS} ${DEBUG_FLAGS}

