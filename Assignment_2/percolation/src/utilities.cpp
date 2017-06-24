/*
   utilities.cpp

   Created by Petter Taule on 27.03.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <cmath>
#include <armadillo>
#include "utilities.hpp"
#include "percolation.hpp"

using namespace arma;

// Change to c-read i.e. fscan etc??
void readToMatrix(Mat<int>& bonds, const std::string& filename) {
    std::ifstream input(filename);
    if (input.fail()) {
        printf("Could not open: %s\n", filename.c_str());
        return;
    }
    std::string line;

    // First line in file is number of nodes
    // The number is separated by an = sign
    getline(input,line);
    std::size_t pos = line.find_last_of("=");
    NODES = stoi(line.substr(pos+1));
    // Second line in file is number of bonds
    getline(input,line); pos = line.find_last_of("=");
    BONDS = stoi(line.substr(pos+1));

    bonds.set_size(BONDS,2);
    for (int i = 0; i < BONDS; ++i) {
        getline(input,line);
        std::stringstream ss(line);
        for (int j = 0; j < 2; ++j) {
            ss >> bonds(i,j);
        }
    }
    input.close();
}


void shuffleRows(Mat<int>& bonds, std::mt19937_64& generator) {
    const int rows = bonds.n_rows;

    for (int i = 0; i < rows-1; ++i) {
        std::uniform_int_distribution<int> distribution(i+1,rows - 1);
        int swpIndex = distribution(generator);
        bonds.swap_rows(i,swpIndex);
    }
}


int findRoot(const Col<int>& sites, int j) {
    if (sites(j) < 0) return j;
    return findRoot(sites,sites(j));
}


void logBinomialCoefficients(int M, Col<realNum>& log_B) {
    log_B = zeros<Col<realNum>>(M+1);

    for (int i = 1; i <= M; ++i) {
        log_B(i) = std::log(static_cast<realNum>(M-i+1)/i) + log_B(i-1);
    }
}


void writeMeasurements(Col<realNum>& p_inf, Col<realNum>& s, Col<realNum>& chi, const std::string& basename) {
    const std::string filename = "./data/" + basename + ".txt";
    FILE* f = fopen(filename.c_str(),"w");
    printf("MEASUREMENTSFILE = \"%s\"\n",filename.c_str());

    fprintf(f,"# p\t\t\tp_inf(p)\ts(p)\t\tchi(p)\n");
    for (int i = 0; i < P_VALS; ++i) {
        const realNum p = i * (1.0/P_VALS);
        fprintf(f,"%f\t%.10e\t%.10e\t%.10e\n", p, p_inf(i), s(i), chi(i));
    }
}


void findLargestCluster(const Col<int>& sites, const int largestClusterIdx) {
    std::vector<bool> partOfCluster(NODES);
    for (int i = 0; i < NODES; ++i) {
        if (findRoot(sites,i) == largestClusterIdx) partOfCluster[i] = true;
    }
    // writeLattice(partOfCluster,"partOfCluster_" + std::to_string(P));
}


void writeLattice(const std::vector<bool>& partOfCluster, const std::string& basename) {
    FILE* f;
    const std::string filename = "./data/" + basename + ".txt";
    f = fopen(filename.c_str(),"w");
    printf("LATTICEFILE = \"%s\"\n",filename.c_str());

    // TODO: any better way to do this?
    const int L = std::sqrt(NODES);
    printf("L = %i\n", L);
    int x = 0;
    int y = 0;
    for (unsigned int i = 0; i < partOfCluster.size(); ++i) {
        if (partOfCluster[i]) {
            x = i % L;
            y = i / L;
            fprintf(f, "%i\t%i\n", x,y);
        }
    }
    fclose(f);
}

/* // The Col<realNum>& argument is a vector which stores */
/* // logarithms of natural numbers up to n, in order to */
/* // reduce computation time */
/* realNum logBinomial(int n, int k, const Col<realNum>& log_numbers) { */
/*     realNum result = 0; */
/*     for (int i = 2; i <= k; ++i) { */
/*         result -= log_numbers(i); */
/*     } */
/*     for (int i = n-k+1; i <= n; ++i) { */
/*         result += log_numbers(i); */
/*     } */
/*     return result; */
/* } */

/* int binomial(int n, int k) { */
/*     assert(n >= k); */
/*     if (k == 0) { */
/*         return 1; */
/*     } */
/*     else if (k == 1 || k == n) { */
/*         return n; */
/*     } */
/*     return binomial(n, k - 1) * (n-k+1) / k; */
/* } */

/* void writeLogBinCoefficients(int M) { */
/*     const std::string filename = "./constants/logBin_" + std::to_string(M) + ".txt"; */
/*     FILE* f = fopen(filename.c_str(),"w"); */

/*     Col<realNum> log_numbers = zeros<Col<realNum>>(M+1); */

/*     // The zero index value is of no interest */
/*     for (int n = 1; n <= M; ++n) { */
/*         log_numbers(n) = std::log(n); */
/*     } */
/*     for (int n = 0; n < M; ++n) { */
/*         fprintf(f, "%.10e\n", logBinomial(M,n,log_numbers)); */
/*     } */
/* } */

/* void readLogBinCoefficients(Col<realNum>& log_B, const std::string& filename) { */
/*     std::ifstream input(filename); */
/*     if (input.fail()) { */
/*         printf("Could not open: %s\n", filename.c_str()); */
/*         return; */
/*     } */
/*     for (int i = 0; i < BONDS; ++i) { */
/*         input >> log_B(i); */
/*     } */
/* } */
