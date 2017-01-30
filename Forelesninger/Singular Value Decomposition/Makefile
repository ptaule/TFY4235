CXX = g++
FLAGS = -Wall -O2 --std=c++11
LIBS = -lblas -llapack -larmadillo
DEBUG_FLAGS = -ggdb

all: svd

svd: svd.cpp
	${CXX} $^ -o $@ ${FLAGS} ${LIBS}

debug: svd.cpp
	${CXX} $^ -o $@ ${FLAGS} ${LIBS} ${DEBUG_FLAGS}

