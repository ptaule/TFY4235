#
#   reflecting.jl
#
#   Created by Petter Taule on 06.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("diffusion.jl")

using PyPlot

const N        = 50
const h        = 1/N
const δt       = 1e-3
const D        = 1
const u_0      = 1
const cfl      = D * δt/h^2

const x0_index = Int(round(N/2))

const X        = linspace(0,1,N)
const x_0      = X[x0_index]

# CFL number
println("CFL = $(cfl)")

# Initial concentration
U = zeros(N)
U[x0_index] = u_0
U[x0_index-1] = u_0
U[x0_index+1] = u_0

# Times
time = 0.01
const timeDiff = 0.03
const numSteps = Int(round(timeDiff/δt))
const adjFactor = 3

# First plot
U     = Neumann.crankNicolson(U,Int(round(time/δt)),h,δt,D)
U_ana = Neumann.analyticSol(collect(X),time,D,3*u_0*h,x_0,50)

plot(X, U     , "-", label=L"t="*"$(round(time , 2))")
plot(X, U_ana , ".", color="Black", ms=0.8)

# Plots
for i=1:3
    time += timeDiff
    U     = Neumann.crankNicolson(U,numSteps,h,δt,D)
    U_ana = Neumann.analyticSol(collect(X),time,D,3*u_0*h,x_0,50)

    plot(X, U     , "-" , label=L"t="*"$(round(time,2))")
    plot(X, U_ana , "." , color="Black", ms=0.8)
end

# Save plot
const basename = "../res/reflecting/" * "crank_Delta" * "_CFL$(cfl)"

legend()
xlabel(L"x")
savefig(basename*".svg")

# Write info
f = open(basename*".log", "w")
println(f,"N   = $N   " );
println(f,"δt  = $δt  " );
println(f,"D   = $D   " );
println(f,"u_0 = $u_0 " );
println(f,"cfl = $cfl " );
println(f,"deltafunction from three middle points")
close(f)
