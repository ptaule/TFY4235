#
#   diffusion.jl
#
#   Created by Petter Taule on 06.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

module Dirichlet

function eulerExplicit(
    U         ::Vector,
    timeSteps ::Integer,
    δx        ::Real,
    δt        ::Real,
    D         ::Real
    )

    const S = size(U,1)
    const α = D * δt / (δx)^2

    M = SymTridiagonal(fill(1-2α,S), fill(α,S-1))

    for i=1:timeSteps
        U = M*U
    end
    return U
end


function eulerImplicit(
    U         ::Vector,
    timeSteps ::Integer,
    δx        ::Real,
    δt        ::Real,
    D         ::Real
    )

    const S = size(U,1)
    const α = D * δt / (δx)^2

    M = SymTridiagonal(fill(1+2α,S), fill(-α,S-1))
    M_inv = inv(M)

    for i=1:timeSteps
        U = M_inv*U
    end
    return U
end


function crankNicolson(
    U         ::Vector,
    timeSteps ::Integer,
    δx        ::Real,
    δt        ::Real,
    D         ::Real
    )

    const S = size(U,1)
    const α = D * δt / (δx)^2

    A = SymTridiagonal(fill(1+α,S), fill(-α/2,S-1))
    B = SymTridiagonal(fill(1-α,S), fill( α/2,S-1))

    M = inv(A) * B

    for i=1:timeSteps
        U = M*U
    end
    return U
end


function analyticSol(
    X          ::Vector,
    t          ::Real,
    D          ::Real,
    u_0        ::Real,
    x_0        ::Real,
    iterations ::Integer,
    L          ::Real = 1
    )

    const S = size(X,1)
    U = zeros(S)

    v(x,n) = sqrt(2/L) * sin(n*pi*x/L)

    for n=1:iterations
        U += exp( - (n*pi/L)^2 * D * t) * v(X,n) * v(x_0,n)
    end
    U *= u_0
    return U
end
end


module Neumann

function eulerExplicit(
    U         ::Vector,
    timeSteps ::Integer,
    δx        ::Real,
    δt        ::Real,
    D         ::Real
    )

    const S = size(U,1)
    const α = D * δt / (δx)^2

    M = Tridiagonal(fill(α,S-1), fill(1-2α,S), fill(α,S-1))
    M[1,2]   += α
    M[S,S-1] += α

    for i=1:timeSteps
        U = M*U
    end
    return U
end


function eulerImplicit(
    U         ::Vector,
    timeSteps ::Integer,
    δx        ::Real,
    δt        ::Real,
    D         ::Real
    )

    const S = size(U,1)
    const α = D * δt / (δx)^2

    M = Tridiagonal(fill(-α,S-1), fill(1+2α,S), fill(-α,S-1))
    M[1,2]   -= α
    M[S,S-1] -= α

    M_inv = inv(M)

    for i=1:timeSteps
        U = M_inv*U
    end
    return U
end


function crankNicolson(
    U         ::Vector,
    timeSteps ::Integer,
    δx        ::Real,
    δt        ::Real,
    D         ::Real
    )

    const S = size(U,1)
    const α = D * δt / (δx)^2

    A = Tridiagonal(fill(-α/2,S-1), fill(1+α,S), fill(-α/2,S-1))
    B = Tridiagonal(fill( α/2,S-1), fill(1-α,S), fill( α/2,S-1))
    A[1,2]   -= α/2
    A[S,S-1] -= α/2
    B[1,2]   += α/2
    B[S,S-1] += α/2

    M = inv(A) * B

    for i=1:timeSteps
        U = M*U
    end
    return U
end


function analyticSol(
    X          ::Vector,
    t          ::Real,
    D          ::Real,
    u_0        ::Real,
    x_0        ::Real,
    iterations ::Integer,
    L          ::Real = 1
    )

    const S = size(X,1)
    U = zeros(S)

    v(x,n) = (n == 0) ? (1/sqrt(L)) : (sqrt(2/L) * cos(n*pi*x/L))

    for n=0:iterations
        U += exp( - (n*pi/L)^2 * D * t) * v(X,n) * v(x_0,n)
    end
    U *= u_0
    return U
end
end

module VariableD

function crankNicolson(
    U         ::Vector,
    timeSteps ::Integer,
    δx        ::Real,
    δt        ::Real,
    D         ::Vector
    )

    # D should have boundary points, and hence be
    # 2 elements larger
    (size(U,1) == (size(D,1) - 2)) || throw(DimensionMismatch())

    const S = size(U,1)
    const α = δt / (4*(δx)^2)
    println("α  = $α ");

    dl = zeros(S-1) # Lower diagonal
    dA = zeros(S)   # Main diagonal for A
    dB = zeros(S)   # Main diagonal for B
    du = zeros(S-1) # Upper diagonal

    for i=1:S-1
        # Shifting D by one, because of boundary point
        dl[i] = α * (D[i]   + D[i+1])
        du[i] = α * (D[i+1] + D[i+2])
    end
    for i=1:S
        dA[i] = 1 + α*(D[i] + 2 * D[i+1] + D[i+2])
        dB[i] = 1 - α*(D[i] + 2 * D[i+1] + D[i+2])
    end

    A = Tridiagonal(-dl,dA,-du)
    B = Tridiagonal(dl,dB,du)

    M = inv(A) * B

    for i=1:timeSteps
        U = M*U
    end
    return U
end


function analyticStepD(
    X          ::Vector,
    t          ::Real,
    D_minus    ::Real,
    D_plus     ::Real,
    u_0        ::Real,
    x_0        ::Real,
    )

    const S = size(X,1)
    U = zeros(S)

    A_plus = 2 * (1 + sqrt(D_minus/D_plus))^(-1)
    A_minus = A_plus * sqrt(D_minus/D_plus)

    for j=1:S
        if (X[j] >= x_0)
            U[j] += ( A_plus / sqrt(4*pi*D_plus*t)
                    * exp( - (X[j] - x_0)^2 / (4*D_plus*t))
                    )
        else
            U[j] += ( A_minus / sqrt(4*pi*D_minus*t)
                    * exp( - (X[j] - x_0)^2 / (4*D_minus*t))
                    )
        end
    end
    U *= u_0
    return U
end

end
