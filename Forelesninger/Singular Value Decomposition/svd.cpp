/*
   svd.cpp

   Created by Petter Taule on 30.01.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/
#include<iostream>
#include<armadillo>

using namespace arma;

int main(int argc, char *argv[]) {

    mat X(2,3);
    X(0,0) = 3;
    X(1,0) = -1;
    X(0,1) = 1;
    X(1,1) = 3;
    X(0,2) = 1;
    X(1,2) = 1;

    mat U;
    mat V;
    vec s;

    svd(U,s,V,X);

    U.print();
    V.print();
    s.print();

    return 0;
}
