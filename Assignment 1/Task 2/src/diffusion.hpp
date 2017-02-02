/*
   diffusion.hpp

   Created by Petter Taule on 02.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#ifndef diffusion_hpp
#define diffusion_hpp

#include<iostream>
#include<armadillo>

using realNum = double;

const int TOTTIME      = 100; // Total time of simulation
const int N            = 10;  // Number of timesteps

const int LENGTH       = 1;   // Size of system (1 dimensjon)
const int INTERVALS    = 10;  // Partial intervals

const realNum D        = 1;


void euler_explicit(const int size);
realNum getALPHA();

#endif /* ifndef diffusion_hpp */
