#
#   findRoots.jl
#
#   Created by Petter Taule on 27.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

using PyPlot, Roots

const nu_0 = 1000


k(λ) = sqrt(λ)
κ(λ) = sqrt(nu_0 - λ)

f(λ) = exp(κ(λ)/3) .* ( κ(λ) .* sin(k(λ)/3) + k(λ) .* cos(k(λ)/3) ).^2 - exp(-κ(λ)/3) .* ( κ(λ) .* sin(k(λ)/3) - k(λ) .* cos(k(λ)/3) ).^2

root = @time fzeros(f,0,nu_0)
println("λ_0 = $root")
println("f(λ_0) = $(f(root))")

λ = linspace(0,nu_0,1000)

plot(λ, f(λ),label=L"f(\lambda)")
scatter(root,zeros(length(root)))
xlabel(L"$\lambda$")
legend()
savefig("test.svg")

