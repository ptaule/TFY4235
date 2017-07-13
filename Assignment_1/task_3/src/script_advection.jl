#
#   script_advection.jl
#
#   Created by Petter Taule on 30.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("advection.jl")

using PyPlot

u_0(x) = max(1-16*(x-0.5).^2,0) # Initial function
u(x,t) = u_0(x - t)      # Exact solution

U = u_0(x)

const plotDiff = 0.1
const timeSteps = Int(round(plotDiff/Δt))

const basename = "../res/downwind_CFL<1"

# Plot some values
plot(x,U,label="t = 0")
plot(x,u(x,0),color="black",linestyle="dashed",linewidth=0.7)
for i = 1:3
    t = i*plotDiff
    U = downwind!(U,timeSteps)
    plot(x,U,label="t = $(round(t,2))")
    plot(x,u(x,t),color="black",linestyle="dashed",linewidth=0.7)
end

xlabel(L"$x$")
ylabel(L"$u$")
legend()
savefig(basename * ".svg")

# Write log

f = open(basename * ".info", "w")
println(f, "c  = $c")
println(f, "Δx = $Δx")
println(f, "Δt = $Δt")
println(f, "α  = $α ")
close(f)
