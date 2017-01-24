/*
   code.cpp

   Created by Petter Taule on 24.01.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include<iostream>
#include<random>
#include<armadillo>

using namespace arma;

const int N = 1000;

Col<int> distribution = zeros<Col<int>>(N);

int main(int argc, char *argv[]) {
    std::random_device seeder;
    std::default_random_engine generator(seeder());
    std::uniform_int_distribution<int> uniform(0,1);

    for (int i = 0; i < N; ++i) {
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

    distribution.save("distribution.txt", raw_ascii);
    return 0;
}
