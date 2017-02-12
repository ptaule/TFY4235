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

const int TOTTIME   = 10;   // Total time of simulation
const int N         = 10000; // Number of timesteps
const int LENGTH    = 10;   // Size of system (1 dimensjon)
const int INTERVALS = 1000;  // Partial intervals
const realNum D     = 1;     // Diffusion constant
const realNum U_0   = 1;     // Initial density

realNum getALPHA();
void thomasAlgorithm(const arma::Mat<realNum>& A, const arma::Col<realNum>& B, arma::Col<realNum>& X, const unsigned int size);



namespace dirichlet {
    void eulerExplicit(arma::Col<realNum>& U, const int size, const int timeSteps);
    void eulerImplicit(arma::Col<realNum>& U, const int size, const int timeSteps);
    void crankNicolson(arma::Col<realNum>& U, const int size, const int timeSteps);
    void writeResults(arma::Col<realNum>& U, const std::string& basename, const int number, const int timeStep);
}



// Note that these functions change the size of U;
// element 0 and N+1 is added to create appropriate
// boundary conditions
namespace neumann {
    void eulerExplicit(arma::Col<realNum>& U, int size, const int timeSteps);
    void eulerImplicit(arma::Col<realNum>& U, int size, const int timeSteps);
    void crankNicolson(arma::Col<realNum>& U, int size, const int timeSteps);
    void writeResults(arma::Col<realNum>& U, const std::string& basename, const int number, const int timeStep);
}

#endif /* ifndef diffusion_hpp */
