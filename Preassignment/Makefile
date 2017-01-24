SOURCE = code.cpp
CC = g++
FLAGS = -Wall -O2 --std=c++11
LIBS = -lblas -llapack -larmadillo
DEBUG_FLAGS = -ggdb

.PHONY: code
code:
	${CC} ${SOURCE} -o $@ ${FLAGS} ${LIBS}
	./$@

# build:
# 	${CC} ${FLAGS} ${SOURCE}

# debug:
# 	${CC} ${FLAGS} ${DEBUG_FLAGS} ${SOURCE}

