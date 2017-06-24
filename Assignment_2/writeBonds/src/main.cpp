/*
   main.cpp

   Created by Petter Taule on 28.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include <string>
#include "writeBonds.hpp"


int main() {
    const int L = 1000;
    const std::string L_string = std::to_string(L);

    /* writeSquareBonds(L,"Square_" + L_string + "x" + L_string); */
    writeTriangularBonds(L,"Triangular_" + L_string + "x" + L_string);
    /* writeHonneyBonds(L,"Honney_" + L_string + "x" + L_string); */
    return 0;
}
