#!/usr/bin/env bash

cd src/

make randomWalk
./randomWalk

gnuplot randomWalk.gpi

cd ../res

epstopdf plot.eps
latexmk randomWalk.tex
