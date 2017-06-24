/*
   main.cpp

   Created by Petter Taule on 28.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include <armadillo>
#include <cstdlib>
#include <istream>
#include <string>

#include "utilities.hpp"
#include "percolation.hpp"

int main() {
    std::random_device seeder;
    std::mt19937_64 generator(seeder());

    arma::Mat<int> bonds;
    int L = 0;
    std::cout << "Write length size (L):" << std::endl;
    std::cin >> L;
    const std::string bondsFile = "../../writeBonds/src/data/Triangular_"
        + std::to_string(L) + "x" + std::to_string(L) + ".txt";
    readToMatrix(bonds, bondsFile);

    printf("NODES = %i\n", NODES);
    printf("BONDS = %i\n", BONDS);
    /* printf("PMARK = \"%02i\"\n", (int)(P*10)); */
    percolation(bonds, generator,"triangular/triangular_");
}
