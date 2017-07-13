#
#   eigenvectors.jl
#
#   Created by Petter Taule on 01.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("spectral.jl")

using PyPlot

const N      = 70  # Number of steps (outer interval)
const h      = 1/N # Step size
const n_max  = 10  # Number of eigenvalues/vectors computed

const X      = linspace(0,1,N)
const Y      = linspace(0,1,N)

const n = 4

# Analytic solution

u_ana = - analyticHelmholtz(:+,2,2,N)

# Numeric solution
(_, eigvecs) = @time helmholtz2D(N-2,n_max)

u_num = reindexArrayTo2D(eigvecs[:,n])
u_num = normalizeAndAddBoundary!(u_num)

println("error = $(integrateDiff2D(u_num,u_ana))")

contourf(X,Y,u_num, cmap=ColorMap("winter"))
xlabel(L"x")
ylabel(L"y")
savefig("num_contour_$n.svg")
clf()
surf(X,Y,u_num, cmap=ColorMap("winter"))
xlabel(L"x")
ylabel(L"y")
savefig("num_3D_$n.svg")
clf()

contourf(X,Y,u_ana, cmap=ColorMap("winter"))
xlabel(L"x")
ylabel(L"y")
savefig("ana_contour_$n.svg")
clf()
surf(X,Y,u_ana, cmap=ColorMap("winter"))
xlabel(L"x")
ylabel(L"y")
savefig("ana_3D_$n.svg")
