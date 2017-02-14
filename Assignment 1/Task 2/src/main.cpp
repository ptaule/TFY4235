/*
   main.cpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<iostream>
#include"diffusion.hpp"
#include"calculator.hpp"

using namespace arma;

int main() {
    printf("TOTTIME   = %i\n", TOTTIME   );
    printf("N         = %i\n", N         );
    printf("LENGTH    = %i\n", LENGTH    );
    printf("INTERVALS = %i\n", INTERVALS );
    printf("D_0       = %f\n", D_0       );
    printf("U_0       = %f\n", U_0       );
    printf("ALPHA     = %f\n", getALPHA());

    Col<realNum> U = zeros<Col<realNum>>(INTERVALS);
    U(INTERVALS/2) = U_0;

    Col<realNum> D = zeros<Col<realNum>>(INTERVALS);
    D.fill(D_plus);

    for (int i = 0; i < INTERVALS/2; ++i) {
        D(i) = D_minus;
    }


    /* dirichlet::eulerExplicit(U, INTERVALS, N); */
    /* dirichlet::eulerImplicit(U, INTERVALS, N); */
    /* dirichlet::crankNicolson(U, INTERVALS, N); */
    dirichlet::crankNicolson(U, D, INTERVALS, N);

    /* neumann::eulerExplicit(U, INTERVALS, N); */
    /* neumann::eulerImplicit(U, INTERVALS, N); */
    /* neumann::crankNicolson(U, INTERVALS, N); */

    return 0;
}
