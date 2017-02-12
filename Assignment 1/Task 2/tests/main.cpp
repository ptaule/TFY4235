/*
   main.cpp

   Created by Petter Taule on 12.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<iostream>
#include"tests.hpp"


int main() {
    /* sumTest_neumannExplicit(); */
    /* sumTest_neumannImplicit(); */
    sumTest_neumannCrankNicolsohn();
    thomasAlgorithm_identityTest();
    thomasAlgorithm_test1();
    thomasAlgorithm_test2();
    return 0;
}
