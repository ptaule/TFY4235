/*
   percolation.cpp

   Created by Petter Taule on 28.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <thread>
#include <armadillo>
#include "percolation.hpp"
#include "utilities.hpp"

using namespace arma;

// Initializing extern variables
int NODES = 0;
int BONDS = 0;

// Defining global vectors:

// Column vectors with indecies n=1..M
Col<realNum> p_inf_n;
Col<realNum> p_inf2_n;
Col<realNum> s_n;
Col<realNum> log_B;

// Column vectors with indecies i=1..P_VALS
Col<realNum> p_inf_i;
Col<realNum> p_inf2_i;
Col<realNum> s_i;
Col<realNum> chi;


void percolation(Mat<int>& bonds, std::mt19937_64& generator, const std::string& lattice) {

    // Initializing column vectors
    p_inf_n  = zeros<Col<realNum>>(BONDS);
    p_inf2_n = zeros<Col<realNum>>(BONDS);
    s_n      = zeros<Col<realNum>>(BONDS);

    p_inf_i  = zeros<Col<realNum>>(P_VALS);
    p_inf2_i = zeros<Col<realNum>>(P_VALS);
    s_i      = zeros<Col<realNum>>(P_VALS);
    chi      = zeros<Col<realNum>>(P_VALS);

    Col<int> sites(NODES);

    // Launch thread calculating logarithm of binomial
    // coefficients
    std::thread t1(logBinomialCoefficients, BONDS, std::ref(log_B));

    printf("Start looping over system realizations\n");
    // Looping over system realizations
    for (int k = 0; k < REALIZATIONS; ++k) {
        shuffleRows(bonds, generator);

        // Initialize sites
        sites.fill(-1);

        long int average_s = NODES;
        // int largestClusterIdx = 0;
        int largestClusterSize = 1;

        for (int n = 0; n < BONDS; ++n) {
            const int root1 = findRoot(sites,bonds(n,0));
            const int root2 = findRoot(sites,bonds(n,1));

            if (root1 != root2){
                // Subtract squared cluster sizes from average_s
                average_s -= (std::pow(sites(root1),2) + std::pow(sites(root2),2));

                // If cluster1 is larger than cluster 2 (sizes are negative)
                if (sites(root1) < sites(root2)) {
                    sites(root1) += sites(root2);
                    sites(root2) = root1;
                    // Add new merged cluster to average_s
                    average_s += std::pow(sites(root1),2);
                    if (-sites(root1) > largestClusterSize) {
                        largestClusterSize = -sites(root1);
                        // largestClusterIdx = root1;
                    }
                }
                else {
                    sites(root2) += sites(root1);
                    sites(root1) = root2;
                    average_s += std::pow(sites(root2),2);
                    if (-sites(root2) > largestClusterSize) {
                        largestClusterSize = -sites(root2);
                        // largestClusterIdx = root2;
                    }
                }
            }

            const realNum p_inf0 = static_cast<realNum>(largestClusterSize)/NODES;
            // Avoiding diverging s_n(n) when p_inf_n -> N
            if (p_inf0 < 1 - epsilon) {
                s_n(n) += static_cast<realNum>(average_s -
                        std::pow(largestClusterSize,2))/(NODES - largestClusterSize);
            }
           p_inf_n(n)  += p_inf0;
           p_inf2_n(n) += p_inf0*p_inf0;
        }
    }

    // Averaging quantities
    p_inf_n /= REALIZATIONS;
    p_inf2_n /= REALIZATIONS;
    s_n /= REALIZATIONS;

    printf("Done looping over system realizations\n");

    // Joining thread t1
    t1.join();

    printf("Start looping over p-values\n");
    std::thread threads[numThreads];
    const int interval = P_VALS/numThreads;

    for (int i = 0; i < numThreads; ++i) {
        threads[i] = std::thread(convoluteResults, i*interval, (i+1)*interval);
    }

    // Joining threads
    for (auto& it: threads) {
        it.join();
    }
    printf("Done looping over p-values\n");

    // Using abs to avoid negatie numbers due to "floating fluctuations"
    // when the expression inside the sqrt() is very small
    for (int i = 0; i < P_VALS; ++i) {
        chi(i) = NODES * std::sqrt(std::abs(p_inf2_i(i) - p_inf_i(i)*p_inf_i(i)));
    }

    // Write results
    printf("Writing results.\n");
    writeMeasurements(p_inf_i, s_i, chi, lattice + std::to_string(NODES));
}


void convoluteResults(int start, int stop) {
    // Avoiding p=0:
    if (start == 0) {
        start++;
    }

    for (int i = start; i < stop; ++i) {
        const realNum p = i * (1.0/P_VALS);
        const realNum log_p = std::log(p);
        const realNum log_p_comp = std::log(1-p);
        for (int n = 1; n <= BONDS; ++n) {
            const realNum log_of_factor = log_B(n) + n*log_p + (BONDS - n)*log_p_comp;
            // Using n-1 since the other vectors are zero indexed
            p_inf_i(i)  += std::exp(log_of_factor) * p_inf_n(n-1);
            p_inf2_i(i) += std::exp(log_of_factor) * p_inf2_n(n-1);
            s_i(i)      += std::exp(log_of_factor) * s_n(n-1);
        }
    }
}
