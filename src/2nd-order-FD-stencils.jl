# First derivative (2nd order)
function first_derivative(f, n, dx)
    return (
        + f[n + 1]
        - f[n - 1]
    ) / (2 * dx)
end

# Second derivative (2nd order)
function second_derivative(f, n, dx)
    return (
        + f[n + 1]
        - 2 * f[n]
        + f[n - 1]
    ) / (dx^2)
end

# Third derivative (2nd order)
function third_derivative(f, n, dx)
    return (
        - f[n + 2]
        + 2 * f[n + 1]
        - 2 * f[n - 1]
        + f[n - 2]
    ) / (2 * dx^3)
end
