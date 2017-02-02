/*
   diffusion.cpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<cmath>
#include"diffusion.hpp"

using namespace arma;

const realNum ALPHA = (D*static_cast<realNum>(TOTTIME)/N)/std::pow(static_cast<realNum>(LENGTH)/INTERVALS,2);

void euler_explicit(const int size) {
    SpMat<int>     A = speye<SpMat<int>>(size,size);
    SpMat<realNum> B = zeros<SpMat<realNum>>(size,size);
    Col<realNum> U = zeros<Col<realNum>>(size);

    B.diag() += 1 - 2*ALPHA;
    B.diag(1) += ALPHA;
    B.diag(-1) += ALPHA;

}

realNum getALPHA() {
    return ALPHA;
}
