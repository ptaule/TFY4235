/*
   calculator.cpp

   Created by Petter Taule on 14.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<string>
#include<cmath>
#include"diffusion.hpp"
#include"calculator.hpp"

using namespace arma;



void solveReflectiveBCs(arma::Col<realNum>& U, const unsigned int iterations, const int time) {
    U = zeros<Col<realNum>>(U.size());
    const unsigned int size = U.size();
    Col<realNum> X(size);
    for (unsigned int i = 0; i < size; ++i) {
        X(i) = static_cast<realNum>(i)/size;
    }

    U += U_0 * std::sqrt(1.0/LENGTH);

    for (unsigned int i = 0; i < iterations; i += 2) {
        realNum iterFactor = std::pow(-1,i) * U_0 * std::sqrt(1.0/(2*LENGTH)) * std::exp(-std::pow(i*M_PI/LENGTH,2)*D_0*time);
        for (unsigned int j = 0; j < size; ++j) {
            U(j) += iterFactor * std::cos(i*M_PI*X(j));
        }
    }
}



void solveAbsorbingBCs(arma::Col<realNum>& U, const unsigned int iterations, const int time) {
    U = zeros<Col<realNum>>(U.size());
    const unsigned int size = U.size();
    Col<realNum> X(size);
    for (unsigned int i = 0; i < size; ++i) {
        X(i) = static_cast<realNum>(i)/size;
    }

    U += U_0 * std::sqrt(1.0/LENGTH);

    for (unsigned int i = 0; i < iterations; i += 2) {
        realNum iterFactor = std::pow(-1,i+1) * U_0 * std::sqrt(1.0/(2*LENGTH)) * std::exp(-std::pow(i*M_PI/LENGTH,2)*D_0*time);
        for (unsigned int j = 0; j < size; ++j) {
            U(j) += iterFactor * std::sin(i*M_PI*X(j));
        }
    }
}

void writeResults(Col<realNum>& U, const std::string& basename) {
    FILE* results;
    std::string filename = "./data/" + basename + ".txt";
    results = fopen(filename.c_str(),"w");

    const unsigned int size = U.size();
    realNum x = 0;
    for (unsigned int i = 0; i < size; ++i) {
        x = static_cast<realNum>(i)/size;
        fprintf(results,"%f\t%f\n",x,U(i));
    }
    fclose(results);
}
