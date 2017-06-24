#
#   wavefunctions.jl
#
#   Created by Petter Taule on 24.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("utilities.jl")

using PyPlot, LaTeXStrings

# Constants and initial conditions

const N     = 100            # Number of discretization steps
const h     = 1/N            # Step length
const Δt    = 1e-2           # Time step
const x     = linspace(0,1,N) # x vector
const n_max = 20              # Number of eigenvalues/eigenvectors calculated
const A     = sqrt(2)         # Amplitude
const NU_0  = Int32(1e3)      # Barrier size
const L_A   = 1//3            # Left barrier boundary
const L_B   = 2//3            # Right barrier boundary


println("CLF = $(Δt/h^2)")


function numericalHelmholtz1D(size, numValues)

    M = spzeros(Int8,size,size)

    M += spdiagm(-2*ones(Int8,size))
    M += spdiagm((ones(Int8,size-1),ones(Int8,size-1)), (-1,1))

    return eigs(M, nev=numValues, which=:SM)
end


function boxPotential(size, numValues)

    M = spzeros(size,size)

    nu = zeros(size)
    interval = Int(ceil(L_A * size)) : Int(ceil(L_B * size))
    nu[interval] += NU_0*h^2

    M += spdiagm(-2*ones(size) - nu)
    M += spdiagm((ones(Int8,size-1),ones(Int8,size-1)), (-1,1))

    return eigs(M, nev=numValues, which=:SM)
end


function forwardEuler(Psi_0::Vector, timeSteps)
    const S = size(Psi_0,1)
    const α = Δt/h^2

    Psi = Psi_0

    d = fill(1 - 2*im*α, S)
    interval = Int(ceil(L_A * S)) : Int(ceil(L_B * S))
    d[interval] -= Δt * im * NU_0

    d1 = fill(im*α,S-1)

    M = SymTridiagonal(d,d1)

    for i = 1:timeSteps
        Psi = M*Psi
    end
    return Psi
end


function crankNicolson(Psi_0::Vector, timeSteps)
    const S = size(Psi_0,1)
    const α = Δt/h^2
    Psi = Psi_0

    d_M = fill(1 + im * α, S)
    d_B = fill(1 - im * α, S)

    interval = Int(ceil(L_A * N)) : Int(ceil(L_B * N))
    d_M[interval] += im * Δt * NU_0 /2
    d_B[interval] -= im * Δt * NU_0 /2

    d1 = fill(im/2 * α, S-1)

    M = Tridiagonal(-d1,d_M,-d1)
    B = Tridiagonal(d1,d_B,d1)

    A = inv(M) * B

    for i = 1:timeSteps
        Psi = A * Psi
    end
    return Psi
end


function evolveWavefunction{T<:Number, S<:Real}(Psi_0::Array{T,1}, eigenvecs::Array{S,2}, lambda::Array{S,1}, time)
    n = size(eigenvecs,2)

    Psi = zeros(N)

    for i=1:n
        alpha = innerProduct(Psi_0,eigenvecs[:,i])
        Psi += alpha * exp( - lambda[i] * time * im) * eigenvecs[:,i]
    end
    return Psi
end


function normalizeEigVecs(vecs::Array{RealNum,2})
    for i = 1:size(vecs,2)
        vecs[:,i] /= sqrt(trapz(vecs[:,i].^2))
    end
    return vecs
end


function plotEigVecs(vecs::Matrix, numPlots, basename::String)
    filename = "../res/eigvecs_" * basename * "_" * string(N) * ".svg"

    for n = 1:numPlots
        plot(x,vecs[:,n],label=(L"n = "*string(n)))
    end

    xlabel(L"$x$")
    legend()
    savefig(filename)
    clf()
end


function plotEigVals(vals::Vector, basename::String)
    filename = "../res/eigvals_" * basename * "_" * string(N) * ".svg"

    n = linspace(1,n_max,n_max)
    plot(n,vals, "b.", label="Numerical eigenvalues")
    xlabel(L"n")
    ylabel(L"E_n/E_0")
    legend()
    savefig(filename)
    clf()
end


function eigenvecError(vecs::Array{RealNum,2},f::Function)
    for n=1:n_max
        sigma = trapz(abs(f(n) - vecs[:,n]))
        @printf "Integrating error for vector %i: %f\n" n sigma
    end
end


function readAndPlotEigvals(filename::String, f::Function)
    vals = readdlm(filename)
    # Enforcing 1 dimension on eigenvalues
    vals = vals[:,1]

    n = linspace(1,length(vals),length(vals))

    plt = plot(n, vals, linewidth=0.5, "b", label="Numerical values")
    plot(n, f(n), linewidth=0.5, "g", label="Analytic solution")
    xlabel=(L"$n$")
    ylabel=(L"$E_n/E_0$")
    legend()
    savefig("../res/eigvals_box_"*string(N)*".svg")
    clf()
end
