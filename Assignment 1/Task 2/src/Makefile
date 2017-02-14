CXX = g++
FLAGS = -Wall -Wextra -Wpedantic -O3 --std=c++11
LIBS = -lblas -llapack -larmadillo
DEBUG_FLAGS = -ggdb

# Reminder:
# $@ - target
# $^ - all dependencies
# $< - first dependency

plot: run
	cat diffusion.log gnuPlotSettings.txt > diffusion.gpi
	gnuplot5 diffusion.gpi

run: all
	time ./diffusion.prog > diffusion.log

all: diffusion.prog

diffusion.prog: main.o diffusion.o calculator.o
	${CXX} $^ -o $@ ${LIBS}

main.o: main.cpp diffusion.hpp calculator.hpp
	${CXX} $< -c ${FLAGS}

diffusion.o: diffusion.cpp diffusion.hpp
	${CXX} $< -c ${FLAGS}

calculator.o: calculator.cpp calculator.hpp diffusion.hpp
	${CXX} $< -c ${FLAGS}

debug: main.cpp diffusion.cpp
	${CXX} $^ -o debug.prog ${FLAGS} ${LIBS} ${DEBUG_FLAGS}

plot_only:
	gnuplot5 diffusion.gpi

clean:
	rm -f *.o
