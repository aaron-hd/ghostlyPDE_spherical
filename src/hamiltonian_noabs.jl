
using NumericalIntegration
using SpecialFunctions


function nintegrate_simps(y::Vector, h::Number)
    n = length(y)-1
    n % 2 == 0 || error("`y` length (number of intervals) must be odd")
    s = - sum(y[1:2:n] + 4*y[2:2:n] + y[3:2:n+1])
    return h/3 * s
end;



#   Now define a function that determines the Hamiltonian density.

function hamiltonian_density(fields_at_t, pHyper)
    
    # unpack parameters
    (pGrid, pModel) = pHyper
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # unpack fields
    N = size(fields_at_t, 1)    
    y0, y1, y2, y3 = [
            fields_at_t[1:N], 
            fields_at_t[N+1:2N], 
            fields_at_t[2N+1:3N], 
            fields_at_t[3N+1:4N]
        ]
    
    # surface area of (n-1)-dimensional unit sphere
    # A_n = 2 * π^(Ndim/2) / gamma(Ndim/2)
    A_n = 1
    
    # prepare and calculate the hamiltonian density    
    # NOTE: to achieve periodicity for the integrator, 
    # we need to include one of the boundary points
    ham = zeros(Nx - 1)  # needs to be odd for Simpson rule
    for n = NboundaryPadding + 1:Nx + NboundaryPadding - 1
        
        # current index
        idx = n
        # assign the left boundary index        
        boundary_idx = NboundaryPadding + 1 
        
        # reconstruct the current value of r
        rn = (idx - NboundaryPadding - 1) * dx
        
        # calculate Hamiltonian density
        ham[n - NboundaryPadding] = rn^(Ndim-1)* A_n *(
            # phi
            + 1//2 * y1[idx]^2
            + 1//2 * first_derivative_even_at_left_boundary(y0, idx, dx, boundary_idx)^2
            # chi
            + sigma * 1//2 * y3[idx]^2
            + sigma * 1//2 * first_derivative_even_at_left_boundary(y2, idx, dx, boundary_idx)^2
            # masses
            + mphi2 * y0[idx]^2 / 2 
            + sigma * mchi2 * y2[idx]^2 / 2 
            # interactions            
            + lambdaCross * y0[idx]^mV * y2[idx]^nV
            + lambdaSelf * y0[idx]^lV
            + sigma * lambdaSelf * y2[idx]^lV
        )
    end
    ham
end




function hamiltonian_phi(fields_at_t, pHyper)
    
    # unpack parameters
    (pGrid, pModel) = pHyper
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # unpack fields
    N = size(fields_at_t, 1)    
    y0, y1 = [
            fields_at_t[1:N], 
            fields_at_t[N+1:2N]
        ]
    
    # prepare and calculate the hamiltonian density
    hamphi = zeros(Nx - 1)  # needs to be odd for Simpson rule
    for n = NboundaryPadding + 1:Nx + NboundaryPadding - 1
        
        # current index
        idx = n
        
        # reconstruct the current value of r
        rn = (idx - NboundaryPadding - 1) * dx
        
        # calculate Hamiltonian density
        hamphi[n - NboundaryPadding] = rn^(Ndim-1)*(
            # phi
            + 1//2 * y1[idx]^2
            + 1//2 * first_derivative(y0, idx, dx)^2
            # mass phi
            + mphi2 * y0[idx]^2 / 2 
            + lambdaSelf * y0[idx]^lV
        )
    end
    hamphi
end





function hamiltonian_chi(fields_at_t, pHyper)
    
    # unpack parameters
    (pGrid, pModel) = pHyper
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # unpack fields
    N = size(fields_at_t, 1)    
    y2, y3 = [
            fields_at_t[2N+1:3N], 
            fields_at_t[3N+1:4N]
        ]
    
    # prepare and calculate the hamiltonian density
    hamchi = zeros(Nx - 1)  # needs to be odd for Simpson rule
    for n = NboundaryPadding + 1:Nx + NboundaryPadding - 1
        
        # current index
        idx = n
        
        # reconstruct the current value of r
        rn = (idx - NboundaryPadding - 1) * dx
        
        # calculate Hamiltonian density
        hamchi[n-NboundaryPadding] = rn^(Ndim-1)*(
            # kinetic chi
            + sigma * 1//2 * y3[idx]^2
            + sigma * 1//2 * first_derivative(y2, idx, dx)^2
            # mass phi
            + sigma * mchi2 * y2[idx]^2 / 2 
            + sigma * lambdaSelf * y2[idx]^lV
        )
    end
    hamchi
end




function power_radiated_phi(fields_at_t, pHyper)
    
    # unpack parameters
    (pGrid, pModel) = pHyper
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # unpack fields
    N = size(fields_at_t, 1)    
    y0, y1 = [
            fields_at_t[1:N], 
            fields_at_t[N+1:2N]
        ]
    
    # index and radius at the boundary
    idx = Nx + NboundaryPadding
    rn = (idx - NboundaryPadding - 1) * dx
    
    # energy flux density at boundary
    F_boundary_phi = y1[idx] * first_derivative(y0, idx, dx)

    # surface area of (n-1)-dimensional unit sphere
    # A_n = 2 * π^(Ndim/2) / gamma(Ndim/2)
    A_n = 1

    # power through boundary
    P_out_phi = A_n * rn^(Ndim-1) * F_boundary_phi
    
    return P_out_phi
    
end

function power_radiated_chi(fields_at_t, pHyper)
    
    # unpack parameters
    (pGrid, pModel) = pHyper
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # unpack fields
    N = size(fields_at_t, 1)    
    y2, y3 = [
            fields_at_t[2N+1:3N], 
            fields_at_t[3N+1:4N]
        ]
    
    # index and radius at the boundary
    idx = Nx + NboundaryPadding
    rn = (idx - NboundaryPadding - 1) * dx
    
    # energy flux density at boundary
    F_boundary_chi = sigma * y3[idx] * first_derivative(y2, idx, dx)

    # surface area of (n-1)-dimensional unit sphere
    # A_n = 2 * π^(Ndim/2) / gamma(Ndim/2)
    A_n = 1

    # power through boundary
    P_out_chi = A_n * rn^(Ndim-1) * F_boundary_chi
    
    return P_out_chi
    
end





function nintegrate_simps(y::Vector, h::Number)
    n = length(y)-1
    n % 2 == 0 || error("`y` length (number of intervals) must be odd")
    s = - sum(y[1:2:n] + 4*y[2:2:n] + y[3:2:n+1])
    return h/3 * s
end;











