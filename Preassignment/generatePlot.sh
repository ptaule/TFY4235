#!/usr/bin/env bash

cd src/

make randomWalk

gnuplot randomWalk.gpi

cd ../res

epstopdf plot.eps
latexmk randomWalk.tex
