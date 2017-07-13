#
#   advection.jl
#
#   Created by Petter Taule on 30.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

const N = 50
const Δx = 1/N
const Δt = 0.01
const c = 1
const α = c * Δt / (2 * Δx)

const x = linspace(Δx,1,N)

function upwind!{T<:Real}(U::Array{T,1}, timeSteps::Integer)
    const S = size(U,1)
    println("S          = $S");
    println("CFL-number = $(2*α)")

    M = zeros(S,S)

    d  = fill( 1-2α ,S   )
    d1 = fill( 2α   ,S-1 )

    M += diagm(d)
    M += diagm(d1,-1)
    M[1,S] = 2α

    for i=1:timeSteps
        U = M * U
    end
    return U

    # Matrix:
    # | 1-2α  0    ...       2α |
    # |  2α  1-2α  ...       0  |
    # |  0    2α   ...       0  |
    # |  .     .    .           |
    # |  0     .           1-2α |
end


function downwind!{T<:Real}(U::Array{T,1}, timeSteps::Integer)
    const S = size(U,1)
    M = zeros(S,S)
    println("S          = $S \n");
    println("CFL-number = $(2*α)")

    d  = fill( 1+2α ,S   )
    d1 = fill( -2α  ,S-1 )

    M += diagm(d)
    M += diagm(d1,-1)
    M[S,1] = 2α

    for i=1:timeSteps
        U = M * U
    end
    return U

    # Matrix:
    # | 1+2α -2α    ...       0  |
    # |  0  1+2α    ...       0  |
    # |  0    0     ...       0  |
    # |  .     .     .           |
    # |  2α    .            1+2α |
end



function explicitCentered!{T<:Real}(U::Array{T,1}, timeSteps::Integer)
    const S = size(U,1)
    M = zeros(S,S)
    println("S  = $S \n");

    d  = ones(S)
    d1 = fill(α,S-1)

    M += diagm(d)
    M += diagm(-d1,1)
    M += diagm(d1,-1)
    M[1,S] = α
    M[S,1] = -α

    for i=1:timeSteps
        U = M * U
    end
    return U

    # Matrix:
    # |  1  -α     ...      α |
    # |  α   1     ...      0 |
    # |  0   α     ...      0 |
    # |  .   .     .          |
    # |  -α  .              1 |
end



function implicitCentered!{T<:Real}(U::Array{T,1}, timeSteps::Integer)
    const S = size(U,1)
    M = zeros(S,S)
    println("S  = $S \n");

    d  = ones(S)
    d1 = fill(α,S-1)

    M += diagm(d)
    M += diagm(-d1,1)
    M += diagm(d1,-1)
    M[1,S] = -α
    M[S,1] = α

    M_inv = inv(M)

    for i=1:timeSteps
        U = M_inv * U
    end
    return U

    # Matrix:
    # |  1   α     ...     -α |
    # | -α   1     ...      0 |
    # |  0  -α     ...      0 |
    # |  .   .     .          |
    # |  α   .              1 |
end



function lax_Friedrich!{T<:Real}(U::Array{T,1}, timeSteps::Integer)
    const S = size(U,1)
    M = zeros(S,S)
    println("S  = $S \n");

    d  = zeros(S)
    d1u = fill(1/2 - α,S-1)
    d1l = fill(1/2 + α,S-1)

    M += diagm(d)
    M += diagm(d1u,1)
    M += diagm(d1l,-1)
    M[1,S] = 1/2 + α
    M[S,1] = 1/2 - α

    for i=1:timeSteps
        U = M * U
    end
    return U

    # Matrix:
    # |    0    0.5-α     ...      α |
    # |  0.5+α   0        ...      0 |
    # |    0    0.5+α     ...      0 |
    # |    .     .        .        0 |
    # |  0.5-α   .                 0 |
end



function lax_Wendroff!{T<:Real}(U::Array{T,1}, timeSteps::Integer)
    const S = size(U,1)
    const β = c * α / Δx
    M = zeros(S,S)
    println("S  = $S \n");

    d   = fill( 1 - 2*β , S   )
    d1u = fill(-α + β   , S-1 )
    d1l = fill( α + β   , S-1 )

    M += diagm(d)
    M += diagm(d1u,1)
    M += diagm(d1l,-1)
    M[1,S] =  α + β
    M[S,1] = -α + β

    for i=1:timeSteps
        U = M * U
    end
    return U

    # Matrix:
    # |  1-2β  -α+β        ...    α+β |
    # |  α+β   1-2β        ...      0 |
    # |    0    α+β        ...      0 |
    # |    .     .         .        0 |
    # |  -α+β    .               1-2β |
end

