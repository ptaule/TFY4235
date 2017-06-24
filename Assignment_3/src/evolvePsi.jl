#
#   juliascript.jl
#
#   Created by Petter Taule on 24.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("wavefunctions.jl")

using PyPlot

# ------------------------------------------------ #
# Calculating eigenvectors and eigenvalues

(lambda, eigenvecs) = @time numericalHelmholtz1D(N, n_max)


# Want to rescale and reverse the
# order of the eigenvalues
lambda *= -1/h^2

eigenvecs = normalizeEigVecs(eigenvecs)

# x = linspace(0,1,N)
# f(n) = A * sin(n*pi*x)

# eigenvecError(eigenvecs,f)

# ------------------------------------------------ #

Psi_0 = 1/sqrt(2) * (eigenvecs[:,1] + eigenvecs[:,2])

t1 = pi/(lambda[2] - lambda[1])

T = linspace(0,t1,5)

for t in T
    Psi = evolveWavefunction(Psi_0, eigenvecs, lambda, t)
    plot(x,abs2(Psi),label=string(round(t,0)))
end

xlabel(L"$x$")
legend()
savefig("Psi.svg")
