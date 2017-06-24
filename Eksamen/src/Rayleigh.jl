#
#   Rayleigh.jl
#
#   Created by Petter Taule on 08.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("Integrator.jl")


module Rayleigh

# Only solveRayleigh()-function is exported to external scope
export solveRayleigh

using Integrator

include("constants.jl")

# Enums for Surface profiles and boundary conditions
@enum SurfProfile DOUBLESIN=1 TRUNCCONE=2 TRUNCCOS=3
@enum BC          DIRICHLET=1 NEUMANN=2

# For loops exporting enum values from module
for s in instances(SurfProfile)
    @eval export $(Symbol(s))
end

for s in instances(BC)
    @eval export $(Symbol(s))
end

# Using variable precision with VarFloat
typealias VarFloat Float64


# Alpha function
# k2 and lambda2 corresponds to the squared quantities
α(k2) = ((k2 <= 4*pi^2/λ^2) ? (sqrt(4*pi^2/λ^2 - k2)) : (im*sqrt(k2 - 4*pi^2/λ)))


# M function
M(p1,p2,q1,q2) = (4*pi^2/λ^2 - p1*q1 - p2*q2) * 1/(α(sqrt(q1^2 + q2^2)))
# N function
N(p1,p2,q1,q2) =  - M(p1,p2,q1,q2)


# The matricies uses a 1D indexing to specify points in 2D
# reciprocal space, hence we use two functions to convert
# between indicies
function index2Dto1D(i::Integer, j::Integer, S::Integer)
    return (i + S*(j-1))
end

function index1Dto2D(k::Integer,S::Integer)
    return ((k-1) % S + 1 , div(k-1,S) + 1)
end

function reindexArray(v::Vector)
    const S = Int(sqrt(size(v,1)))
    v2 = zeros(Number,S,S)

    # Julia 2D arrays are stored by column
    for j=1:S, i=1:S
        v2[i,j] = v[index2Dto1D(i,j,S)]
    end
    return v2
end


function initializeLeftMatrix(
    kx :: Real,                    # Incoming (lateral) wavevector,
                                   # first component i real space
    ky :: Real;                    # .. second component in real space
    sp :: SurfProfile = DOUBLESIN, # Surface profile
    bc :: BC          = DIRICHLET  # Boundary condition
    )

    # Size of matrix
    const S = 2H+1

    const β = 2*pi / a
    # Using VarFloat which can be varied
    A = zeros(Complex{VarFloat},S^2,S^2)

    G_p_h1 = [Int(i) for i = -H : H] # G prime, first coordinate
    G_p_h2 = [Int(i) for i = -H : H] # G prime, second coordinate
    G_h1   = [Int(i) for i = -H : H] # G, first coordinate
    G_h2   = [Int(i) for i = -H : H] # G, second coordinate

    for j=1:S^2
        # G prime indecies are determined
        (primeidx_1, primeidx_2) = index1Dto2D(j,S)

        for i=1:S^2
            # G indicies are determined
            (idx_1, idx_2) = index1Dto2D(i,S)

            # K prime in real space
            K_p1 = kx + β * G_p_h1[primeidx_1]
            K_p2 = ky + β * G_p_h2[primeidx_2]
            # K prime squared
            K_p = K_p1^2 + K_p2^2

            h1 = G_h1[idx_1] - G_p_h1[primeidx_1]
            h2 = G_h2[idx_2] - G_p_h2[primeidx_2]

            G = β * sqrt(h1^2 + h2^2)

            if (sp == DOUBLESIN)
                A[i,j] = integrateDoubleSinusoidal(h1,h2,-α(K_p))
            elseif (sp == TRUNCCONE)
                A[i,j] = integrateTruncatedCone(G,-α(K_p))
            else
                A[i,j] = integrateTruncatedCosine(G,-α(K_p))
            end

            if (bc == NEUMANN)
                # K in real space
                K1 = kx + β * G_h1[idx_1]
                K2 = ky + β * G_h2[idx_2]
                A[i,j] *= M(K1,K2,K_p1,K_p2)
            end
        end
    end
    return A
end


function initializeRightVector(
    kx :: Real,                    # Incoming (lateral) wavevector,
                                   # first component i real space
    ky :: Real;                    # .. second component in real space
    sp :: SurfProfile = DOUBLESIN, # Surface profile
    bc :: BC          = DIRICHLET  # Boundary condition
    )

    # Vector size
    const S = 2H +1
    # Using VarFloat which can be varied
    b = zeros(Complex{VarFloat},S^2)

    const β = 2*pi / a
    # k squared in real space
    k2 = kx^2 + ky^2
    gamma = α(k2)

    G_h1 = [Int(i) for i = -H : H] # G, first coordinate
    G_h2 = [Int(i) for i = -H : H] # G, second coordinate

    for j=1:S^2
        (idx_1, idx_2) = index1Dto2D(j,S)

        h1 = G_h1[idx_1]
        h2 = G_h2[idx_2]

        G = β * sqrt(h1^2 + h2^2)

        if (sp == DOUBLESIN)
            b[j] = - integrateDoubleSinusoidal(h1,h2,gamma)
        elseif (sp == TRUNCCONE)
            b[j] = - integrateTruncatedCone(G,gamma)
        else
            b[j] = - integrateTruncatedCosine(G,gamma)
        end
        if (bc == NEUMANN)
            # K in real space
            K1 = kx + β * G_h1[idx_1]
            K2 = ky + β * G_h2[idx_2]
            b[j] *= N(K1,K2,kx,ky)
        end
    end
    return b
end


function diffractionEfficiency(
    kx :: Real,    # Incoming (lateral) wavevector, first component i real space
    ky :: Real,    # .. second component in real space
    r  :: Matrix   # Solution of Rayleigh equation (2D)
    )

    const S = size(r,1)
    # Check that sizes match
    ((2H+1) == S) || throw(DimensionMismatch("Size of r does not correspond to H"))
    eff = zeros(VarFloat,S,S)

    G_h1 = [Int(i) for i = -H : H] # G, first coordinate
    G_h2 = [Int(i) for i = -H : H] # G, second coordinate

    const β  = 2*pi / a
    const k2 = kx^2 + ky^2
    const c  = α(k2)

    # Keeping track of how many points are included
    # in the sum
    numEva = 0
    numPro = 0

    # Julia matricies are stored by column
    for j=1:S, i=1:S
        const K2 = (kx + β * G_h1[i])^2 + (ky + β * G_h2[j])^2

        if (K2 > (4*pi^2/λ^2))
            eff[i,j] = NaN
            numEva += 1
        else
            eff[i,j] = α(K2)/c * abs2(r[i,j])
            numPro += 1
        end
    end
    println("numEva  = $numEva ");
    println("numPro  = $numPro ");
    return eff
end


function solveRayleigh(
    kx :: Real,                    # Incoming (lateral) wavevector,
                                   # first component i real space
    ky :: Real;                    # .. second component in real space
    sp :: SurfProfile = DOUBLESIN, # Surface profile
    bc :: BC          = DIRICHLET  # Boundary condition
    )

    # Initialize matrix and vector
    println("Initializing right vector")
    @time b = initializeRightVector(kx,ky,sp=sp,bc=bc)
    println("Initializing left matrix")
    @time A = initializeLeftMatrix(kx,ky,sp=sp,bc=bc)
    println("Solving system")
    # Using LAPACK gesv function for solving linear system
    # The b vector is replaced with the solution
    @time Base.LinAlg.LAPACK.gesv!(A,b)

    r = reindexArray(b)

    eff = diffractionEfficiency(kx,ky,r)
    return (r,eff)
end

end
