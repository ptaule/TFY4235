/*
   writeBonds.hpp

   Created by Petter Taule on 28.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#ifndef writeBonds_hpp
#define writeBonds_hpp

#include <string>

void writeSquareBonds(const int L, const std::string& basename);
void writeTriangularBonds(const int L, const std::string& basename);
void writeHonneyBonds(const int L, const std::string& basename);

#endif /* ifndef writeBonds_hpp */
