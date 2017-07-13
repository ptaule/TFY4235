#
#   eigenvalues.jl
#
#   Created by Petter Taule on 29.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("spectral.jl")

using PyPlot

const N      = 100  # Number of steps
const h      = 1/N # Step size
const n_max  = 100   # Number of eigenvalues/vectors computed

const exact = pi^2 * sort(vec([i^2+j^2 for i=1:10, j=1:10]))

const n = linspace(1,100,100)
const sizes = [10,30,50,70,90]

for i in sizes
    (eigvals, eigvecs) = @time helmholtz2D(i,n_max)
    semilogy(n,eigvals,".",label=L"$N=$"*"$i", ms=0.8)
end


semilogy(n,exact, "-", color="black", linewidth=0.6, label="Exact")
xlabel(L"$n$")
ylabel(L"$\lambda$")
legend()
savefig("test.svg")
