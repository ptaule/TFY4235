/*
   calculator.hpp

   Created by Petter Taule on 14.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#ifndef calculator_hpp
#define calculator_hpp

#include<armadillo>

void solveReflectiveBCs(arma::Col<realNum>& U, const unsigned int iterations, const int time);
void solveAbsorbingBCs(arma::Col<realNum>& U, const unsigned int iterations, const int time);
void writeResults(arma::Col<realNum>& U, const std::string& basename);

#endif /* ifndef calculator_hpp */
