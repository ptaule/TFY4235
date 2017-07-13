#
#   hopf_burgers.jl
#
#   Created by Petter Taule on 06.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#


function lax_Wendroff_hopf(U_in::Vector, δt::Real, δx::Real, timeSteps::Integer)
    const S = size(U_in,1)
    const α = δt/δx

    U = U_in
    V = zeros(S)
    for i=1:timeSteps
        for j=2:S-1
            V[j] =( + U[j]
                    - α/4   * (U[j+1]^2 - U[j-1]^2)
                    + α^2/8 * (U[j+1] + U[j]  ) * (U[j+1]^2 - U[j]^2  )
                    - α^2/8 * (U[j]   + U[j-1]) * (U[j]^2   - U[j-1]^2)
            )
        end
        # Periodic boundary conditions
        V[1] = (+ U[1]
                - α/4   * (U[2]^2 - U[S]^2)
                + α^2/8 * (U[2] + U[1]) * (U[2]^2 - U[1]^2)
                - α^2/8 * (U[1] + U[S]) * (U[1]^2 - U[S]^2)
        )
        V[S] = (+ U[S]
                - α/4   * (U[1]^2 - U[S-1]^2)
                + α^2/8 * (U[1] + U[S]  ) * (U[1]^2 - U[S]^2  )
                - α^2/8 * (U[S] + U[S-1]) * (U[S]^2 - U[S-1]^2)
        )
        U = V
    end
    return U
end


function iterative_hopf(X::Vector, t::Real, u_0::Function, epsilon::Real = 1e-3)
    const S = size(X,1)
    U = zeros(S)

    for i=1:S
        u = u_0(X[i])
        while (true)
            u_ = u_0(X[i] - u*t)
            if (abs(u_ - u) < epsilon)
                u = u_
                break
            end
            u = u_
        end
        U[i] = u
    end
    return U
end
