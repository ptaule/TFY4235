#
#   trunccone_script.jl
#
#   Created by Petter Taule on 09.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("constants.jl")
include("Rayleigh.jl")

using Rayleigh, PyPlot

# Incident wavevector, components in real space
const kx = 2*pi/λ * sin(θ) * cos(ϕ)
const ky = 2*pi/λ * sin(θ) * sin(ϕ)


# Diffraction efficiencies as function of angle and m

function plotDiffractionEfficiencies(basename::String)
    const angles = linspace(-pi/2 + 0.01, pi/2 - 0.01,N_n)
    const G = [Int(i) for i = -H : H]
    Eff = zeros(N_n,2H+1)
    U   = zeros(N_n)

    for i = 1:N_n
        kx = 2*pi/λ * sin(angles[i]) * cos(ϕ)
        ky = 2*pi/λ * sin(angles[i]) * sin(ϕ)

        (_,eff) = solveRayleigh(kx,ky,sp=TRUNCCOS,bc=NEUMANN)
        Eff[i,:] = eff[:,H+1]
        # println("U       = $(sum(eff))")
        U[i] = sum([x for x in eff if !isnan(x)])
        println("U[$i] = $(U[i])")
    end

    writedlm(basename*"_e.txt",Eff)
    writedlm(basename*"_U.txt",U)

    surf(G,2*angles/pi,Eff)
    ax = gca()
    xlabel(L"m")
    xlim(-10,10)
    ylabel(L"\theta_0/(\pi/2)")
    zlabel(L"e_m(\theta)")
    yticks([-1,-0.5,0,0.5,1])
    savefig(basename * "_e.svg")
    clf()

    plot(angles/pi,U)
    xlabel(L"\theta/\pi")
    ylabel(L"\mathcal{U}")
    savefig(basename * "_U.svg")
end

# # -------------------------------------------------------"

# Plotting r and e

# (r,eff) = solveRayleigh(kx,ky,sp=TRUNCCOS)

# println("sum = $(sum([x for x in eff if !isnan(x)]))");

# X = linspace(-H,H,2*H+1)

# plot_wireframe(X,X,eff,cmap=ColorMap("winter"))
# xlabel(L"x")
# ylabel(L"y")
# savefig("eff.svg")
# clf()

# surf(X,X,abs(r),cmap=ColorMap("winter"))
# xlabel(L"x")
# ylabel(L"y")
# savefig("r.svg")

plotDiffractionEfficiencies(basename)
