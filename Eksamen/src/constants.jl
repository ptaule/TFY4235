#
#   constants.jl
#
#   Created by Petter Taule on 09.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

# Constants
const H          = 20
const λ          = 1
const a_over_λ   = 3.5
const ζ_0_over_λ = 0.3
const ρ_0        = NaN
const ρ_b        = NaN
const ρ_t        = NaN

const θ          = pi/4
const ϕ          = 0

const ζ_0 = ζ_0_over_λ * λ
const a   = a_over_λ   * λ


# Number of theta values
const N_n = 20

const basename = "../doubleSin/" * "_$(a)_$(N_n)"

# Write parameters to file
f = open(basename * ".log","w")
println(f, "H   = $H")
println(f, "λ   = $λ")
println(f, "ζ_0 = $ζ_0")
println(f, "a   = $a")
println(f, "θ   = $θ")
println(f, "ϕ   = $ϕ")
println(f, "ρ_0 = $ρ_0")
println(f, "ρ_b = $ρ_b")
println(f, "ρ_t = $ρ_t")
println(f, "N_n = $N_n")
close(f)
