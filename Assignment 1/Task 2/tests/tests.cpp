/*
   tests.cpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<iostream>
#include"../src/diffusion.hpp"

using namespace arma;

const realNum epsilon = 0.001;

void cflConstant() {
    const realNum alpha = getALPHA();
    if (alpha > 0.5) {
        std::cout << "Warning! Alpha = " << alpha << " > 0.5." << std::endl;
    }
    else {
        std::cout << "Alpha = " << alpha << std::endl;
    }

}


void sumConcentration() {
    Col<realNum> U = zeros<Col<realNum>>(INTERVALS);
    U(INTERVALS/2) = 1;
    eulerExplicit(U, INTERVALS, N);
    std::cout << "Summen av U er " << sum(U) << std::endl;
    std::cout << "U_0 = " << U_0 << std::endl;
}


int main(int argc, char *argv[]) {
    cflConstant();
    sumConcentration();
    return 0;
}
