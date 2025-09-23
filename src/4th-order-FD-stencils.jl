#####################################################
# INTERIOR
##################################################### 

# First derivative (4th order)
function first_derivative(f, n, dx)
    return (
        - f[n + 2]
        + 8 * f[n + 1]
        - 8 * f[n - 1]
        + f[n - 2]
    ) / (12 * dx)
end

# Second derivative (4th order)
function second_derivative(f, n, dx)
    return (
        - f[n + 2]
        + 16 * f[n + 1]
        - 30 * f[n]
        + 16 * f[n - 1]
        - f[n - 2]
    ) / (12 * dx^2)
end

# Third derivative (4th order)
function third_derivative(f, n, dx)
    return (
        - f[n + 3]
        + 8 * f[n + 2]
        - 13 * f[n + 1]
        + 13 * f[n - 1]
        - 8 * f[n - 2]
        + f[n - 3]
    ) / (8 * dx^3)
end


#####################################################
# LEFT BOUNDARY
#####################################################  

# First derivative (4th order); left boundary version
function first_derivative_even_at_left_boundary(f, n, dx, boundary_idx)
    result = (
        - f[idx_even_at_left_boundary(n + 2, boundary_idx)]
        + 8 * f[idx_even_at_left_boundary(n + 1, boundary_idx)]
        - 8 * f[idx_even_at_left_boundary(n - 1, boundary_idx)]
        + f[idx_even_at_left_boundary(n - 2, boundary_idx)]
    ) / (12 * dx)
    return result
end

# Second derivative (4th order); left boundary version
function second_derivative_even_at_left_boundary(f, n, dx, boundary_idx)
    return (
        - f[idx_even_at_left_boundary(n + 2, boundary_idx)]
        + 16 * f[idx_even_at_left_boundary(n + 1, boundary_idx)]
        - 30 * f[idx_even_at_left_boundary(n, boundary_idx)]
        + 16 * f[idx_even_at_left_boundary(n - 1, boundary_idx)]
        - f[idx_even_at_left_boundary(n - 2, boundary_idx)]
    ) / (12 * dx^2)
end

# Helper function to reflect indices at the left boundary
function idx_even_at_left_boundary(n, boundary_idx)
    # (note that julia starts counting indices at 1)
    return abs(n-boundary_idx) + boundary_idx
end