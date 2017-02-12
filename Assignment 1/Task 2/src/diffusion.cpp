/*
   diffusion.cpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<stdio.h>
#include<cmath>
#include"diffusion.hpp"

using namespace arma;

const realNum ALPHA = (D*static_cast<realNum>(TOTTIME)/N)/std::pow(static_cast<realNum>(LENGTH)/INTERVALS,2);

realNum getALPHA() {
    return ALPHA;
}

void thomasAlgorithm(const Mat<realNum>& A, const Col<realNum>& B, Col<realNum>& X, const unsigned int size) {
    const int n = size - 1;

    if (B.size() != size) {
        printf("Ulik størrelse!");
        return;
    }

    Col<realNum> C = zeros<Mat<realNum>>(size);
    Col<realNum> D = zeros<Mat<realNum>>(size);

    C(0) = A(0,1)/A(0,0);
    D(0) = B(0)/A(0,0);

    // Calculating intermediate coefficients
    for (int i = 1; i < n; ++i) {
        C(i) = A(i,i+1) / (A(i,i) - A(i,i-1)*C(i-1));
        D(i) = (B(i) - A(i,i-1)*D(i-1)) / (A(i,i) - A(i,i-1)*C(i-1));
    }
    D(n) = (B(n) - A(n,n-1)*D(n-1)) / (A(n,n) - A(n,n-1)*C(n-1));

    // Back substitution
    X(n) = D(n);
    for (int i = n - 1; i >= 0; --i) {
        X(i) = D(i) - C(i)*X(i+1);
    }
}

void dirichlet::eulerExplicit(Col<realNum>& U, const int size, const int timeSteps) {
    SpMat<int>     A = speye<SpMat<int>>(size,size);
    SpMat<realNum> B = zeros<SpMat<realNum>>(size,size);

    B.diag() += 1 - 2*ALPHA;
    B.diag(1) += ALPHA;
    B.diag(-1) += ALPHA;

    int j = 1;
    for (int i = 1; i < timeSteps; ++i) {
        U = B*U;
        if (i % (timeSteps/4) == 0) {
            dirichlet::writeResults(U, "dirichlet_eulerExplicit", j, i);
            ++j;
        }
    }
    dirichlet::writeResults(U, "dirichlet_eulerExplicit", j, timeSteps);
}


void dirichlet::eulerImplicit(Col<realNum>& U, const int size, const int timeSteps) {
    Mat<realNum> A = zeros<Mat<realNum>>(size,size);
    Mat<int>     B = eye<Mat<int>>(size,size);

    A.diag() += 1 + 2*ALPHA;
    A.diag(1) += -ALPHA;
    A.diag(-1) += -ALPHA;

    int j = 1;
    for (int i = 1; i < timeSteps; ++i) {
        thomasAlgorithm(A,(B*U),U,size);
        if (i % (timeSteps/4) == 0) {
            dirichlet::writeResults(U, "dirichlet_eulerImplicit", j, i);
            ++j;
        }
    }
    dirichlet::writeResults(U, "dirichlet_eulerImplicit", j, timeSteps);
}

void dirichlet::crankNicolson(Col<realNum>& U, const int size, const int timeSteps) {
    Mat<realNum> A = zeros<Mat<realNum>>(size,size);
    Mat<realNum> B = zeros<Mat<realNum>>(size,size);

    A.diag() += 1 + ALPHA;
    A.diag(1) += -ALPHA/2;
    A.diag(-1) += -ALPHA/2;

    B.diag() += 1 - ALPHA;
    B.diag(1) += ALPHA/2;
    B.diag(-1) += ALPHA/2;

    int j = 1;
    for (int i = 1; i < timeSteps; ++i) {
        thomasAlgorithm(A,(B*U),U,size);
        if (i % (timeSteps/4) == 0) {
            dirichlet::writeResults(U, "dirichlet_crankNicolson", j, i);
            ++j;
        }
    }
    dirichlet::writeResults(U, "dirichlet_crankNicolson", j, timeSteps);
}

void dirichlet::writeResults(Col<realNum>& U, const std::string& basename, const int number, const int timeStep) {
    FILE* results;
    std::string filename = basename + "_" + std::to_string(number) + ".txt";
    results = fopen(filename.c_str(),"w");

    realNum deltaX = static_cast<realNum>(LENGTH)/INTERVALS;
    realNum x = -deltaX;

    for (unsigned int i = 0; i < U.size(); ++i) {
        x += deltaX;
        fprintf(results,"%f\t%f\n",x,U(i));
    }
    printf("Wrote %s at timestep %i\n",filename.c_str(),timeStep);
    fclose(results);
}

void neumann::eulerExplicit(Col<realNum>& U, int size, const int timeSteps) {
    size += 2;
    U.set_size(size);
    SpMat<int>     A = speye<SpMat<int>>(size,size);
    SpMat<realNum> B = zeros<SpMat<realNum>>(size,size);

    B.diag() += 1 - 2*ALPHA;
    B.diag(1) += ALPHA;
    B.diag(-1) += ALPHA;
    B(0,1) = 2*ALPHA;
    B(size-1,size-2) = 2*ALPHA;

    int j = 1;
    for (int i = 1; i < timeSteps; ++i) {
        U = B*U;
        if (i % (timeSteps/4) == 0) {
            neumann::writeResults(U, "neumann_eulerExplicit", j, i);
            ++j;
        }
    }
    neumann::writeResults(U, "neumann_eulerExplicit", j, timeSteps);
}

void neumann::eulerImplicit(Col<realNum>& U, int size, const int timeSteps) {
    size += 2;
    U.set_size(size);
    Mat<realNum> A = zeros<Mat<realNum>>(size,size);
    Mat<int>     B = eye<Mat<int>>(size,size);

    A.diag() += 1 + 2*ALPHA;
    A.diag(1) += -ALPHA;
    A.diag(-1) += -ALPHA;
    A(0,1) = -2*ALPHA;
    A(size-1,size-2) = -2*ALPHA;

    int j = 1;
    for (int i = 1; i < timeSteps; ++i) {
        thomasAlgorithm(A,(B*U),U,size);
        if (i % (timeSteps/4) == 0) {
            neumann::writeResults(U, "neumann_eulerImplicit", j, i);
            ++j;
        }
    }
    neumann::writeResults(U, "neumann_eulerImplicit", j, timeSteps);
}

void neumann::crankNicolson(Col<realNum>& U, int size, const int timeSteps) {
    size += 2;
    U.set_size(size);
    Mat<realNum> A = zeros<Mat<realNum>>(size,size);
    Mat<realNum> B = zeros<Mat<realNum>>(size,size);

    A.diag() += 1 + ALPHA;
    A.diag(1) += -ALPHA/2;
    A.diag(-1) += -ALPHA/2;
    A(0,1) = -ALPHA;
    A(size-1,size-2) = -ALPHA;

    B.diag() += 1 - ALPHA;
    B.diag(1) += ALPHA/2;
    B.diag(-1) += ALPHA/2;
    B(0,1) = ALPHA;
    B(size-1,size-2) = ALPHA;

    int j = 1;
    for (int i = 1; i < timeSteps; ++i) {
        thomasAlgorithm(A,(B*U),U,size);
        if (i % (timeSteps/4) == 0) {
            neumann::writeResults(U, "neumann_crankNicolson", j, i);
            ++j;
        }
    }
    neumann::writeResults(U, "neumann_crankNicolson", j, timeSteps);
}

void neumann::writeResults(Col<realNum>& U, const std::string& basename, const int number, const int timeStep) {
    FILE* results;
    std::string filename = basename + "_" + std::to_string(number) + ".txt";
    results = fopen(filename.c_str(),"w");

    realNum deltaX = static_cast<realNum>(LENGTH)/INTERVALS;
    realNum x = 0;
    for (unsigned int i = 1; i < U.size() - 1; ++i) {
        x += deltaX;
        fprintf(results,"%f\t%f\n",x,U(i));
    }
    printf("Wrote %s at timestep %i\n",filename.c_str(),timeStep);
    fclose(results);
}
