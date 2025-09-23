# definition of the laplacian in spherical coordinates (assuming spherical symmetry)

function laplacian_radial(f, n, dr, r, Ndim)
    return (
        + second_derivative(f, n, dr)
        + (Ndim-1)/r * first_derivative(f, n, dr)
    )
end


function laplacian_radial_even_at_left_boundary(f, n, dr, r, Ndim, boundary_idx)
    return (
        + second_derivative_even_at_left_boundary(f, n, dr, boundary_idx)
        + (Ndim-1)/r * first_derivative_even_at_left_boundary(f, n, dr, boundary_idx)
    )
end