/*
   main.cpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<iostream>
#include"diffusion.hpp"

using namespace arma;

int main() {
    printf("TOTTIME   = %i\n", TOTTIME   );
    printf("N         = %i\n", N         );
    printf("LENGTH    = %i\n", LENGTH    );
    printf("INTERVALS = %i\n", INTERVALS );
    printf("D         = %f\n", D         );
    printf("U_0       = %f\n", U_0       );
    printf("ALPHA     = %f\n", getALPHA());

    Col<realNum> U = zeros<Col<realNum>>(INTERVALS);
    U((INTERVALS/2)) = U_0;

    /* neumann::eulerExplicit(U, INTERVALS, N); */
    /* neumann::eulerImplicit(U, INTERVALS, N); */
    neumann::crankNicolson(U, INTERVALS, N);

    return 0;
}


