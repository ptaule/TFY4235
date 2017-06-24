/*
   percolation.hpp

   Created by Petter Taule on 28.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#ifndef percolation_hpp
#define percolation_hpp

#include <armadillo>
#include <string>
#include <vector>
#include <random>

using realNum = double;
using long_realNum = long double;

// Number of nodes and bonds, set in readToMatrix-function
extern int NODES;
extern int BONDS;

const int REALIZATIONS = 1000; // Number of system realizations
const int P_VALS       = 10000; // Resolution of probability

// Numerical threshold
const realNum epsilon = 1e-6;

// Number of threads
const int numThreads = 3;

void percolation(arma::Mat<int>& bonds, std::mt19937_64& generator, const std::string& lattice);
void convoluteResults(int start, int stop);

#endif /* ifndef percolation_hpp */
