#
#   script_hopf.jl
#
#   Created by Petter Taule on 06.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("hopf_burgers.jl")

using PyPlot

const N = 100
const h = 1/N
const δt = 0.001


# u_0(x) = 1/2 + max(1/2 - 25*(x-1/3)^2,0)
u_0(x) = 1/2 * max(1 - 9*(x-1/2)^2,0)

U = zeros(N)
X = linspace(0,1,N)

for j=1:N
    U[j] = u_0(X[j])
end

const plotDiff = 0.1
const timeSteps = Int(round(plotDiff/δt))

plot(X,U,label=L"$t=0$")
plot(X,U,".", color="Black",ms=0.8)

for i=1:3
    t = i*plotDiff

    U = lax_Wendroff_hopf(U,δt,h,timeSteps)
    U_ana = @time iterative_hopf(collect(X),t,u_0)
    plot(X,U,label="t=$(round(t,2))")
    plot(X,U_ana,".", color="Black",ms=0.8)
end
for i=4:5
    t = i*plotDiff

    U = lax_Wendroff_hopf(U,δt,h,timeSteps)
    plot(X,U,label="t=$(round(t,2))")
end

legend()
xlabel(L"$x$")
savefig("test.svg")
clf()
