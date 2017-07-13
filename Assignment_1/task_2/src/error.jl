#
#   error.jl
#
#   Created by Petter Taule on 01.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

include("spectral.jl")

using PyPlot

# const N      = 10  # Number of steps
# const h      = 1/N # Step size
# const n_max  = 100   # Number of eigenvalues/vectors computed

# const X      = linspace(0,1,N)
# const Y      = linspace(0,1,N)

# Analytic solution

const N_max = 7

const N = [10*i for i=1:N_max]

error = zeros(N_max,5)

i = 0
for n in N
    i += 1
    const X      = linspace(0,1,n)
    const Y      = linspace(0,1,n)
    (_, eigvecs) = @time helmholtz2D(n-2,5)

    # 1. eigenvector
    begin
        u_num = reindexArrayTo2D(eigvecs[:,1])
        u_num = normalizeAndAddBoundary!(u_num)
        u_ana = analyticHelmholtz(:+,1,1,n)
        error[i,1] = integrateDiff2D(u_num, u_ana)
    end
    # 2. eigenvector
    begin
        u_num = reindexArrayTo2D(eigvecs[:,2])
        u_num = normalizeAndAddBoundary!(u_num)
        u_ana = analyticHelmholtz(:+,2,1,n)
        error[i,2] = integrateDiff2D(u_num, u_ana)
    end
    # 3. eigenvector
    begin
        u_num = reindexArrayTo2D(eigvecs[:,3])
        u_num = normalizeAndAddBoundary!(u_num)
        u_ana = analyticHelmholtz(:-,2,1,n)
        error[i,3] = integrateDiff2D(u_num, u_ana)
    end
    # 4. eigenvector
    begin
        u_num = reindexArrayTo2D(eigvecs[:,4])
        u_num = normalizeAndAddBoundary!(u_num)
        u_ana = - analyticHelmholtz(:+,2,2,n)
        error[i,4] = integrateDiff2D(u_num, u_ana)
    end
    # 5. eigenvector
    begin
        u_num = reindexArrayTo2D(eigvecs[:,5])
        u_num = normalizeAndAddBoundary!(u_num)
        u_ana = analyticHelmholtz(:+,2,3,n)
        error[i,5] = integrateDiff2D(u_num, u_ana)
    end

end
writedlm("./data/error.txt",error)

# surf(X,Y,Z, cmap=ColorMap("winter"))
# xlabel(L"x")
# ylabel(L"y")
# savefig("test.svg")
