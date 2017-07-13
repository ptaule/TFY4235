#
#   utilities.jl
#
#   Created by Petter Taule on 05.05.2017
#   Copyright (c) 2017 Petter Taule. All rights reserved.
#

function trapz2D(matrix::Matrix, L_1::Integer = 1, L_2::Integer = 1)
    const n_rows = size(matrix,1)
    const n_cols = size(matrix,2)

    result = 0.0
    # Corners
    result += 1/4 * (matrix[1,1] + matrix[1,n_cols] + matrix[n_rows,1] + matrix[n_rows,n_cols])
    # Edges
    for i=2:(n_rows-1)
        result += 1/2 *(matrix[i,1] + matrix[i,n_cols])
    end
    for j=2:(n_cols-1)
        result += 1/2 *(matrix[1,j] + matrix[n_rows,j])
    end
    # Inner elements
    for j=2:(n_cols - 1), i=2:(n_rows - 1)
        result += matrix[i,j]
    end

    # Multiply by "tile" size
    result *= (L_1/n_rows) * (L_2/n_cols)
    return result
end


function innerProduct2D(A::Matrix, B::Matrix = A)
    # Assuming L_x = L_y = 1
    return trapz2D(real(conj(A).*B))
end

