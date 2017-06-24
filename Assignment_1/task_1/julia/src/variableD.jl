#
#   variableD.jl
#
#   Created by Petter Taule on 06.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("diffusion.jl")

using PyPlot

const N        = 100
const h        = 1/N
const δt       = 1e-5
const u_0      = 1

const x0_index = Int(round(N/2))

const X        = linspace(0,1,N)
const x_0      = X[x0_index]

# Diffusion constant

const D_plus   = 1.5
const D_minus  = 0.5

const D = [fill(D_minus,x0_index-5);
0.6;
0.7;
0.8;
0.9;
1.0;
1.1;
1.2;
1.3;
1.4;
fill(D_plus,x0_index-4)]

# Initial concentration
U = zeros(N)
U[x0_index] = u_0

# Times
time = 0.0005
const timeDiff = 0.001
const numSteps = Int(round(timeDiff/δt))
println("numSteps  = $numSteps ");

# First plot
U[2:end-1] = @time VariableD.crankNicolson(U[2:end-1],Int(round(time/δt)),h,δt,D)
U_ana      = VariableD.analyticStepD(collect(X),time,D_minus,D_plus,u_0*h/2,x_0)

plot(X, U     , "-", label=L"t="*"$(round(time , 2))")
plot(X, U_ana , ".", color="Black", ms=0.8)

# Plots
for i=1:3
    time += timeDiff
    U[2:end-1] = @time VariableD.crankNicolson(U[2:end-1],numSteps,h,δt,D)
    U_ana      = VariableD.analyticStepD(collect(X),time,D_minus,D_plus,u_0*h/2,x_0)

    plot(X, U     , "-" , label=L"t="*"$(round(time,3))")
    plot(X, U_ana , "." , color="Black", ms=0.8)
end

#Save plot
const basename = "../res/variableD/" * "halfD"

legend()
xlabel(L"x")
savefig(basename*".svg")

# Write info
f = open(basename*".log", "w")
println(f,"N   = $N   " );
println(f,"δt  = $δt  " );
println(f,"u_0 = $u_0 " );
println(f,"D_plus    = $D_plus   ");
println(f,"D_minus   = $D_minus  ");
close(f)
