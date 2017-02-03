/*
   main.cpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<iostream>
#include"diffusion.hpp"

int main(int argc, char *argv[]) {
    std::cout << getALPHA();
    euler_explicit(INTERVALS, N);

    return 0;
}
