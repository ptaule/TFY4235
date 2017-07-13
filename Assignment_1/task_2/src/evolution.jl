#
#   evolution.jl
#
#   Created by Petter Taule on 06.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("spectral.jl")
using PyPlot

const N        = 70  # Number of steps (outer interval)
const h        = 1/N # Step size
const n_max    = 30  # Number of eigenvalues/vectors computed

const time     = 0.15
const equation = wave # diffusion, schrodinger or wave

const X        = linspace(0,1,N)
const Y        = linspace(0,1,N)

const sigma    = 0.001
const x_0      = 0.5
const y_0      = 0.5

# Initial condition
f(x,y) = exp(-((x-x_0)^2 + (y-y_0)^2)/sigma)
U_in = zeros(N,N)

for j=1:N, i=1:N
    U_in[i,j] = f(X[i],Y[j])
end

# Numeric solution
(eigvals, eigvecs_) = @time helmholtz2D(N-2,n_max)

eigvecs = zeros(N,N,n_max)
for n=1:n_max
    Z = reindexArrayTo2D(eigvecs_[:,n])
    eigvecs[:,:,n] = normalizeAndAddBoundary!(Z)
end

# Evolve U

U = real(evolveSolution(U_in, eigvecs, eigvals, time, eq=equation))

# Plot initial field
# contourf(X,Y,U_in, cmap=ColorMap("winter"))
# xlabel(L"x")
# ylabel(L"y")
# savefig("init_contour.svg")
# clf()
# surf(X,Y,U_in, cmap=ColorMap("winter"))
# xlabel(L"x")
# ylabel(L"y")
# savefig("init_3D.svg")
# clf()

# Plot final field
contourf(X,Y,U, cmap=ColorMap("winter"))
xlabel(L"x")
ylabel(L"y")
savefig("final_contour_time$time"*"_$equation.svg")
clf()
surf(X,Y,U, cmap=ColorMap("winter"))
xlabel(L"x")
ylabel(L"y")
savefig("final_3D_time$time"*"_$equation.svg")
