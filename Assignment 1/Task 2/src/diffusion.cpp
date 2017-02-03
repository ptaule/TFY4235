/*
   diffusion.cpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<cmath>
#include<armadillo>
#include"diffusion.hpp"

using namespace arma;

const realNum ALPHA = (D*static_cast<realNum>(TOTTIME)/N)/std::pow(static_cast<realNum>(LENGTH)/INTERVALS,2);

void euler_explicit(const int size, const int timeSteps) {
    SpMat<int>     A = speye<SpMat<int>>(size,size);
    SpMat<realNum> B = zeros<SpMat<realNum>>(size,size);
    Mat<realNum>   U = zeros<Mat<realNum>>(size,size);

    // Initializing B matrix
    B.diag() += 1 - 2*ALPHA;
    B.diag(1) += ALPHA;
    B.diag(-1) += ALPHA;

    // Initializing U matrix
    U((size/2),0) = 1;

    for (int i = 1; i < timeSteps; ++i) {
        U.col(i) = B*(U.col(i-1));
    }

    U.print();

}


realNum getALPHA() {
    return ALPHA;
}
