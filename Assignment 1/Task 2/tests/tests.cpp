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


void sumTest_neumannExplicit() {
    Col<realNum> U = zeros<Col<realNum>>(INTERVALS+2);
    U(((INTERVALS+2)/2)) = U_0;

    neumann::eulerExplicit(U, INTERVALS, N);
    realNum diff = std::abs(sum(U) - U_0);
    if (diff > epsilon) {
        std::cout << "Eksplisitt metode: Avvik i masse: " << diff
            << " er større enn epsilon: " << epsilon << std::endl;
    }
}

void sumTest_neumannImplicit() {
    Col<realNum> U = zeros<Col<realNum>>(INTERVALS);
    U((INTERVALS/2)) = U_0;

    neumann::eulerImplicit(U, INTERVALS, N);
    realNum diff = std::abs(sum(U) - U_0);
    if (diff > epsilon) {
        std::cout << "Implisitt metode:Avvik i masse: " << diff
            << " er større enn epsilon: " << epsilon << std::endl;
    }
}

void sumTest_neumannCrankNicolson() {
    Col<realNum> U = zeros<Col<realNum>>(INTERVALS);
    U((INTERVALS/2)) = U_0;

    neumann::crankNicolson(U, INTERVALS, N);
    realNum diff = std::abs(sum(U) - U_0);
    if (diff > epsilon) {
        std::cout << "crankNicolson metode: Avvik i masse: " << diff
            << " er større enn epsilon: " << epsilon << std::endl;
    }
}

void thomasAlgorithm_identityTest() {
    const unsigned int size = 10;
    SpMat<realNum> I = speye<SpMat<realNum>>(size,size);
    Col<realNum> b = ones<Mat<realNum>>(size);
    Col<realNum> X = zeros<Mat<realNum>>(size);

    thomasAlgorithm(I,b,X,size);

    if (!approx_equal(X,(I*b),"absdiff", epsilon)) {
        std::cout << "identitytest failed. " << std::endl;
    }
}

void thomasAlgorithm_test1() {
    const unsigned int size = 10;
    SpMat<realNum> A = zeros<SpMat<realNum>>(size,size);
    Col<realNum> b = ones<Mat<realNum>>(size);
    Col<realNum> X = zeros<Mat<realNum>>(size);

    A.diag() += 3;
    A.diag(-1) += 1;
    A.diag(1) += 1;

    thomasAlgorithm(A,b,X,size);
}

void thomasAlgorithm_test2() {
    const unsigned int size = 10;
    SpMat<realNum> A = zeros<SpMat<realNum>>(size,size);
    Col<realNum> b = ones<Mat<realNum>>(size);
    Col<realNum> X = zeros<Mat<realNum>>(size);

    A.diag() += 3;
    A.diag(-1) += -1;
    A.diag(1) += 2;

    thomasAlgorithm(A,b,X,size);

}
