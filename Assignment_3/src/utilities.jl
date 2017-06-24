#
#   utilities.jl
#
#   Created by Petter Taule on 25.04.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

typealias RealNum Float64

# This function only integrates one the interval [0,1]
function trapz(vec::Array{RealNum})
    # Calculating f(x_1) + 2f(x_2)...2f(x_n-1) + f(x_n)
    result = sum(vec) + sum(vec[2:length(vec)-1])

    result /= 2* (length(vec)-1)
    return result
end

function innerProduct{T<:Number}(a::Array{T, 1}, b::Array{T, 1})
    return trapz(real(conj(a).*b))
end

