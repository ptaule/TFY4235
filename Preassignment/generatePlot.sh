#!/usr/bin/env bash

cd src/

make randomWalk
printf "Running program..."
./randomWalk
printf "done."

gnuplot randomWalk.gpi

cd ../res

epstopdf plot.eps
latexmk randomWalk.tex
