/*
   diffusion.hpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#ifndef diffusion_hpp
#define diffusion_hpp

#include<string>
#include<armadillo>

using realNum = double;

const int TOTTIME     = 10;   // Total time of simulation
const int N           = 10000; // Number of timesteps
const int LENGTH      = 10;   // Size of system (1 dimensjon)
const int INTERVALS   = 50;  // Partial intervals
const realNum D_0     = 1;     // Diffusion constant
const realNum U_0     = 1;     // Initial density

const realNum D_plus  = 0.5;
const realNum D_minus = 1.5;

realNum getALPHA();
void thomasAlgorithm(const arma::SpMat<realNum>& A, const arma::Col<realNum>& B, arma::Col<realNum>& X, const unsigned int size);



namespace dirichlet {
    void eulerExplicit(arma::Col<realNum>& U, const int size, const int timeSteps);
    void eulerImplicit(arma::Col<realNum>& U, const int size, const int timeSteps);
    // Constant diffusion parameter
    void crankNicolson(arma::Col<realNum>& U, const int size, const int timeSteps);
    // Varying diffusion parameter
    void crankNicolson(arma::Col<realNum>& U, const arma::Col<realNum>& D, const int size, const int timeSteps);
    void writeResults(arma::Col<realNum>& U, const std::string& basename, const int number, const double time);
}



namespace neumann {
    void eulerExplicit(arma::Col<realNum>& U, int size, const int timeSteps);
    void eulerImplicit(arma::Col<realNum>& U, int size, const int timeSteps);
    // Constant diffusion parameter
    void crankNicolson(arma::Col<realNum>& U, int size, const int timeSteps);
    // Varying diffusion parameter
    void crankNicolson(arma::Col<realNum>& U, const arma::Col<realNum>& D, const int size, const int timeSteps);
    void writeResults(arma::Col<realNum>& U, const std::string& basename, const int number, const double time);
}

#endif /* ifndef diffusion_hpp */
