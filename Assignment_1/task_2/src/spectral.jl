#
#   spectral.jl
#
#   Created by Petter Taule on 29.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("utilities.jl")

@enum EQUATION diffusion=1 schrodinger=2 wave=3

function helmholtz2D(S::Integer, numValues::Integer)
    # S = size
    # S2 = size^2
    const S2 = S^2
    M  = zeros(Int8     , S2    , S2)
    d1 = zeros(Int8     , S2-1)
    d3 =  fill(Int8(-1) , S2-S)

    for i = 1:S2 - 1
        # Short-circuit syntax
        (i % S != 0) && (d1[i] = -1)
    end

    M += diagm(fill(Int8(4),S2)) # Main diagonal
    M += diagm(d1,1)             # First upper subdiagonal
    M += diagm(d3,S)             # S'th  upper subdiagonal

    # Symmetric eigensolvers are faster
    M = Symmetric(M)

    #Finding eigenvalues and eigenvectors (using LAPACK behind the scenes)
    F = eigfact(M, 1:numValues)
    # Renormalizing eigenvalues
    eigenvals = F[:values] * S2

    return (eigenvals, F[:vectors])
end


function reindexArrayTo2D(arr::Vector)
    S = Int(sqrt(size(arr,1)))
    arr_2 = zeros(S,S)

    # Julia 2D arrays are stored by column
    for j=1:S, i=1:S
            arr_2[i,j] = arr[i + S*(j-1)]
    end
    return arr_2
end


function normalizeAndAddBoundary!(u::Matrix)
    u /= sqrt(innerProduct2D(u))

    # Adding zeros at the boundary
    const size_x = size(u,1)
    u = hcat(u,zeros(size_x))
    u = hcat(zeros(size_x), u)
    const size_y = size(u,2)
    u = vcat(u,zeros(1,size_y))
    u = vcat(zeros(1,size_y),u)
end


function analyticHelmholtz(
    op::Symbol,
    k_x::Integer,
    k_y::Integer,
    N_x::Integer,
    N_y::Integer = N_x;
    L_x::Real = 1,
    L_y::Real = 1
    )

    # Throw exception if k_x == k_y and op=minus
    ((op == :-) && (k_x == k_y)) && throw(ArgumentError("Invalid argument combination."))

    f(x,y,k_x,k_y) = 2* sin(k_x * pi * x / L_x) * sin(k_y * pi * y / L_y)
    const X      = linspace(0,L_x,N_x)
    const Y      = linspace(0,L_y,N_y)
    const factor = (k_x == k_y) ? (1/2) : (1/sqrt(2))
    result       = zeros(N_x,N_y)

    for j=1:N_y, i=1:N_x
        result[i,j] = eval(op)(f(X[i],Y[j],k_x,k_y), f(X[i],Y[j],k_y,k_x))
    end
    result *= factor
    return result
end


function evolveSolution{T<:Real}(U_in::Matrix{T}, eigenmodes::Array{T,3}, eigenvalues::Vector, time::Real; eq::EQUATION = diffusion)
    # Check that sizes are equal
    (size(U_in) == size(eigenmodes,1,2)) || throw(DimensionMismatch("Size mismatch between U_in and eigenmodes"))

    const n = size(eigenvalues,1)
    U       = zeros(size(U_in))

    for i=1:n
        const alpha = innerProduct2D(U_in,eigenmodes[:,:,i])
        exp_factor = 0
        if (eq == diffusion)
            exp_factor = exp(- eigenvalues[i] * time)
        elseif (eq == schrodinger)
            exp_factor = exp( -im * eigenvalues[i] * time)
        else # eq == wave
            exp_factor = exp( -im * sqrt(eigenvalues[i]) * time)
        end
        U += alpha * exp_factor * eigenmodes[:,:,i]
    end

    return U
end


function integrateDiff2D(A::Matrix, B::Matrix)
    return trapz2D(abs(A-B))
end
