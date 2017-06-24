#
#   Integrator.jl
#
#   Created by Petter Taule on 08.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

module Integrator

export integrateDoubleSinusoidal,
       integrateTruncatedCone,
       integrateTruncatedCosine

include("constants.jl")

function integrateDoubleSinusoidal(
    h1    :: Integer, # First component in reciprocal space
    h2    :: Integer, # Second component in reciprocal space
    gamma :: Number,  # Number can be both real and complex
    )

    # besselj is Bessel function of first kind (part of Julia base)
    x = gamma*ζ_0/2
    # besselj throws error for besselj(0,0*im), hence using real part of 0
    if (x == 0)
        x = real(x)
    end
    return ((-1.0*im)^(h1) * besselj(h1,x)
           *(-1.0*im)^(h2) * besselj(h2,x)
           )
end


function integrateTruncatedCone(
    G          :: Real,       # Absolute value of G
    gamma      :: Number,
    iterations :: Integer = 20 # Number of iterations
    )

    # Using zero function to ensure that julia infers
    # correct type
    result = zero(Float64)

    # Add sum/integral part first
    for n = 1:iterations
        f(x) = (ρ_b - (ρ_b - ρ_t)*x) * besselj0( G * (ρ_b - (ρ_b - ρ_t)*x)) * x^n
        # Using quadgk from julia base library
        integral = quadgk(f,0,1)
        # Print if error is large
        result += (-im*gamma*ζ_0)^n/(factorial(n)) * integral[1]
    end
    result *= 2*pi* (ρ_b - ρ_t)/a^2

    # Then, add two first terms
    # Check if G is zero (less than float64 epsilon)
    if (abs(G) < eps())
        result +=    1 + pi *(ρ_t/a)^2 * ( exp(-im*gamma*ζ_0) - 1)
    else
        result += (2*pi *(ρ_t/a)^2 * ( exp(-im*gamma*ζ_0) - 1)
                   * besselj1(G*ρ_t) / (G*ρ_t)
                  )
    end
    return result
end


function integrateTruncatedCosine(
    G          :: Real,    # Absolute value of G
    gamma      :: Number,
    iterations :: Integer = 20  # Number of iterations
    )

    # Truncated cosine function
    S(x) = ζ_0 * cos(pi*x/(2*ρ_0))

    result = zero(Float64)

    for n = 1:iterations
        f(x) = x * besselj0(G*x) * (S(x))^n

        # Using quadgk from julia base library
        integral = quadgk(f,0,ρ_0)
        result += (-im*gamma)^n/factorial(n) * integral[1]
    end
    result *= 2*pi/a^2

    # Add one if G < eps()
    (abs(G) < eps()) ? (result += 1) : ()

    return result
end

end
