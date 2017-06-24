/*
   writeBonds.cpp

   Created by Petter Taule on 27.02.2017
   Copyright (c) 2017 Petter Taule. All rights reserved.
*/

#include <string>
#include <cstdio>
#include <cassert>

void writeSquareBonds(const int L, const std::string& basename) {
    const int N = L*L;
    FILE* f;
    std::string filename = "./data/" + basename + ".txt";
    f = fopen(filename.c_str(),"w");
    printf("Writing to %s\n", filename.c_str());
    fprintf(f, "# Nodes=%i\n# Bonds=%i\n",N,2*N);

    for (int i = 0; i < N-L; ++i) {
        // Check if last column, if true connect to first column
        if ((i+1)% L == 0) {
            fprintf(f,"%i\t%i\n", i,(i - L + 1));
            fprintf(f,"%i\t%i\n", i,i+L);
        }
        else {
            fprintf(f,"%i\t%i\n", i,i+1);
            fprintf(f,"%i\t%i\n", i,i+L);
        }
    }
    // Last row, connect to first row
    for (int i = N-L; i < N - 1; ++i) {
        fprintf(f,"%i\t%i\n", i, i+1);
        fprintf(f,"%i\t%i\n", i, (i-N+L));
    }
    // Connect last element
    fprintf(f,"%i\t%i\n", N-1, (N-L));
    fprintf(f,"%i\t%i\n", N-1, (L-1));
    // Close file
    fclose(f);
}

void writeTriangularBonds(const int L, const std::string& basename) {
    const int N = L*L;
    FILE* f;
    std::string filename = "./data/" + basename + ".txt";
    f = fopen(filename.c_str(),"w");
    printf("Writing to %s\n", filename.c_str());
    fprintf(f, "# Nodes=%i\n# Bonds=%i\n",N,3*N);

    for (int i = 0; i < N-L; ++i) {
        // Check if last column, if true connect to first column
        if ((i+1)% L == 0) {
            fprintf(f,"%i\t%i\n", i,(i - L + 1));
            fprintf(f,"%i\t%i\n", i,(i+L));
            fprintf(f,"%i\t%i\n", i,(i+2));
        }
        else {
            fprintf(f,"%i\t%i\n", i,i+1);
            fprintf(f,"%i\t%i\n", i,i+L);
            fprintf(f,"%i\t%i\n", i,(i+L+1));
        }
    }
    // Last row, connect to first row
    for (int i = N-L; i < N - 1; ++i) {
        fprintf(f,"%i\t%i\n", i, i+1);
        fprintf(f,"%i\t%i\n", i, (i-N+L));
        fprintf(f,"%i\t%i\n", i, (i-N+L+1));
    }
    // Connect last element
    fprintf(f,"%i\t%i\n", N-1, N-L);
    fprintf(f,"%i\t%i\n", N-1, (L-1));
    fprintf(f,"%i\t%i\n", N-1, (0));
    fclose(f);
}

void writeHonneyBonds(const int L, const std::string& basename) {
    // TODO: possible with an odd value of L?
    assert(L % 2 == 0);
    const int N = L*L;
    FILE* f;
    std::string filename = "./data/" + basename + ".txt";
    f = fopen(filename.c_str(),"w");
    printf("Writing to %s\n", filename.c_str());
    fprintf(f, "# Nodes=%i\n# Bonds=%i\n",N,3*N/2);

    int n;

    /* Inner nodes */
    {
        // Every other row have two connections
        for (int i = 0; i < L - 1; i += 2) {
            for (int j = 0; j < L - 1; ++j) {
                n = i*L + j;
                fprintf(f,"%i\t%i\n", n,n+L);
                fprintf(f,"%i\t%i\n", n,(n+L+1));
            }
        }
        // The nodes left have one connection
        for (int i = 1; i < L - 1; i += 2) {
            for (int j = 0; j < L - 1; ++j) {
                n = i*L + j;
                fprintf(f,"%i\t%i\n", n,n+L);
            }
        }
    }

    /* Rightmost nodes */
    {
        // Nodes with two connections
        for (int i = 0; i < L -1; i+=2) {
            n = i*L + L - 1;
            fprintf(f,"%i\t%i\n", n,n+L);
            fprintf(f,"%i\t%i\n", n,(n+1));
        }
        // Nodes with one connection
        for (int i = 1; i < L -1; i+=2) {
            n = i*L + L - 1;
            fprintf(f,"%i\t%i\n", n,n+L);
        }

    }
    /* Buttom nodes */
    {
        // One connection
        for (int i = 0; i < L -1; ++i) {
            n = L*(L-1) + i;
            fprintf(f,"%i\t%i\n", n,(n - N +L));
        }
        // Right-buttom corner
        fprintf(f,"%i\t%i\n", N - 1, L - 1);
    }

    fclose(f);
}
