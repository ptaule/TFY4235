#
#   stepPsi.jl
#
#   Created by Petter Taule on 27.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("wavefunctions.jl")

import PyPlot
plt = PyPlot


# ------------------------------------------------ #
# Calculating first eigenvector

(lambda, eigenvecs) = @time boxPotential(N, n_max)

# Normalizing eigenvectors and eigenvalues

eigenvecs = normalizeEigVecs(eigenvecs)
lambda *= -1/h^2

# ------------------------------------------------ #

Psi = 1/sqrt(2) * (eigenvecs[:,1] + eigenvecs[:,2])
# Psi = zeros(N)
# Psi[Int(round(N/2))] = 1


time = 0
const timeDiff = 300
const numSteps = Int(round(timeDiff/Δt))
println("numSteps  = $numSteps ")

plt.clf()
plot(x,abs2(Psi),label=label=L"t="*"$(round(time,2))")

for i = 1:5
    time += timeDiff
    Psi = @time crankNicolson(Psi, numSteps)
    plot(x,abs2(Psi),label=label=L"t="*"$(round(time,2))")
end

plt.xlabel(L"$x$")
plt.legend()
plt.savefig("crank.svg")
