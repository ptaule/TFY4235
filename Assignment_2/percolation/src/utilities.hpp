/*
   utilities.hpp

   Created by Petter Taule on 27.03.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/
#ifndef utilities_hpp
#define utilities_hpp

#include <random>
#include <armadillo>
#include <string>

#include "percolation.hpp"

void readToMatrix(arma::Mat<int>& bonds, const std::string& filename);

void shuffleRows(arma::Mat<int>& bonds, std::mt19937_64& generator);
int findRoot(const arma::Col<int>& sites, int j);

/* realNum logBinomial(int n, int k, const arma::Col<realNum>& log_numbers); */
void logBinomialCoefficients(int M, arma::Col<realNum>& log_B);
/* int binomial(int n, int k); */
/* void writeLogBinCoefficients(int M); */
/* void readLogBinCoefficients(arma::Col<realNum>& log_B, const std::string& filename); */

void findLargestCluster(const arma::Col<int>& sites, const int largestClusterIdx);

void writeMeasurements(arma::Col<realNum>& p_inf, arma::Col<realNum>& s, arma::Col<realNum>& chi, const std::string& basename);
void writeLattice(const std::vector<bool>& partOfCluster, const std::string& basename);
#endif /* ifndef utilities_hpp */
