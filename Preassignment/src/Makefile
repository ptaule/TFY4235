CC = g++
FLAGS = -Wall -O2 --std=c++11
LIBS = -lblas -llapack -larmadillo
DEBUG_FLAGS = -ggdb

randomWalk: randomWalk.cpp
	${CC} $^ -o $@ ${FLAGS} ${LIBS}

# build:
# 	${CC} ${FLAGS} ${SOURCE}

# debug:
# 	${CC} ${FLAGS} ${DEBUG_FLAGS} ${SOURCE}

