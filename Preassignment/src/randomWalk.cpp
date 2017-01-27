/*
   randomWalk.cpp

   Created by Petter Taule on 24.01.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<iostream>
#include<random>
#include<armadillo>
#include<ctime>
#include<thread>

using namespace arma;

const int N = 10000;
const int threads = 4;

Col<int> distribution = zeros<Col<int>>(N);

int main(int argc, char *argv[]) {
    clock_t time;
    time = clock();

    std::random_device seeder;
    std::default_random_engine generator(seeder());
    std::uniform_int_distribution<int> uniform(0,1);

    for (int i = 0; i < N/4; ++i) {
        int sum = 0;
        int steps = 0;

        do {
            int random = uniform(generator);
            random == 1 ? ++sum : --sum;
            steps++;
        }
        while (sum != 0);

        if (steps < N) {
            distribution(steps)++;
        }
    }

    if(!distribution.save("distribution.txt", raw_ascii)) {
        std::cout << "Error: kunne ikke lagre filen" << std::endl;
    }

    time = clock() - time;
    std::cout << "Programmet brukte " << static_cast<double>(time)/CLOCKS_PER_SEC << " sekunder" << std::endl;
    return 0;
}
