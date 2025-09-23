include("./laplacian_radial.jl")

function finite_differenced_pde_with_bc!(du, u, pHyper, t)
    
    # unpack u and du into separate arrays
    y0, y1, y2, y3 = [u[:, i] for i in 1:4]
    dy0dt, dy1dt, dy2dt, dy3dt = [du[:, i] for i in 1:4]
    
    # unpack parameters
    (pGrid, pModel) = pHyper
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # notation
    dr = dx
    
    # assign the left boundary index        
    boundary_idx = NboundaryPadding + 1 
    
    #####################################################
    # INTERIOR
    #####################################################
    # update non-boundary points
    #(note that NboundaryPadding + 1 is r=0 and will be dealt with below)
    Threads.@threads for n = NboundaryPadding + 2:Nx + NboundaryPadding
        
        # reconstruct the current value of r
        rn = (n - boundary_idx) * dr
        
        # calculate derivatives (which may then be used several times)
        laplacian_phi = laplacian_radial_even_at_left_boundary(y0, n, dr, rn, Ndim, boundary_idx)
        laplacian_chi = laplacian_radial_even_at_left_boundary(y2, n, dr, rn, Ndim, boundary_idx)
        
        dy0dt[n] = (
            # kinetic term
            + y1[n]
        )

        dy1dt[n] = (
            # kinetic term
            + laplacian_phi
            # mass term
            - mphi2 * y0[n] 
            # interaction terms
            - mV * lambdaCross * y0[n] * abs(y0[n])^(mV-2) * abs(y2[n])^nV
            - lV * lambdaSelf * y0[n] * abs(y0[n])^(lV-2)
        )

        dy2dt[n] = (
            # kinetic term
            + y3[n]
        )

        dy3dt[n] = (
            # kinetic term
            + laplacian_chi 
            # mass term
            - mchi2 * y2[n] 
            # interaction terms
            - nV * sigma * lambdaCross * y2[n] * abs(y2[n])^(nV-2) * y0[n]^mV
            - lV * lambdaSelf * y2[n] * abs(y2[n])^(lV-2)
        )
    end
    
    #####################################################
    # INNER BOUNDARY
    #####################################################
    # treat r=0 boundary (enforce regularity and even functions)
    
    # no regularisation necessary
    dy0dt[boundary_idx] = (
        # definition
        + y1[boundary_idx]
    ) 
    # enforce even functions and regularise
    dy1dt[boundary_idx] = (
        # regularised radial laplacian
        + Ndim*second_derivative_even_at_left_boundary(y0,boundary_idx,dr,boundary_idx)
        # mass term
        - mphi2 * y0[boundary_idx] 
        # interaction terms
        - mV * lambdaCross * y0[boundary_idx] * abs(y0[boundary_idx])^(mV-2) * abs(y2[boundary_idx])^nV
        - lV * lambdaSelf * y0[boundary_idx] * abs(y0[boundary_idx])^(lV-2)
    )      
    # no regularisation necessary
    dy2dt[boundary_idx] = (
        # definition
        + y3[boundary_idx]
    )
    # enforce even functions and regularise
    dy3dt[boundary_idx] = (         
        # regularised radial laplacian
        + Ndim*second_derivative_even_at_left_boundary(y2,boundary_idx,dr,boundary_idx)
        # mass term
        - mchi2 * y2[boundary_idx] 
        # interaction terms
        - nV * sigma * lambdaCross * y2[boundary_idx]^(nV-1) * y0[boundary_idx]^mV
        - lV * lambdaSelf * y2[boundary_idx]^(lV-1)
        # interaction terms
        - nV * sigma * lambdaCross * y2[boundary_idx] * abs(y2[boundary_idx])^(nV-2) * y0[boundary_idx]^mV
        - lV * lambdaSelf * y2[boundary_idx] * abs(y2[boundary_idx])^(lV-2)
    )
    
    #####################################################
    # OUTER BOUNDARY
    #####################################################    
    # treat r=\infty boundary points (open boundary conditions)
    for n in 1:NboundaryPadding
        
        # set wave speed and index
        c = 1                                # wave speed
        idx = Nx + NboundaryPadding + n      # current ghost zone index
        
        # reconstruct the current value of r
        rn = (idx - NboundaryPadding - 1) * dr
        
        # or a 4th-order accurate version
        dy0dt[idx] = - c * y0[idx] * (Ndim-1)/(2*rn) + c * (
            -25*y0[idx] 
            + 48*y0[idx-1] 
            - 36*y0[idx-2] 
            + 16*y0[idx-3] 
            - 3*y0[idx-4] 
        ) / (12*dr)
        dy1dt[idx] = - c * y1[idx] * (Ndim-1)/(2*rn) + c * (
            -25*y1[idx] 
            + 48*y1[idx-1] 
            - 36*y1[idx-2] 
            + 16*y1[idx-3] 
            - 3*y1[idx-4] 
        ) / (12*dr)
        dy2dt[idx] = - c * y2[idx] * (Ndim-1)/(2*rn) + c * (
            -25*y2[idx] 
            + 48*y2[idx-1] 
            - 36*y2[idx-2] 
            + 16*y2[idx-3] 
            - 3*y2[idx-4] 
        ) / (12*dr)
        dy3dt[idx] = - c * y3[idx] * (Ndim-1)/(2*rn) + c * ( 
            -25*y3[idx] 
            + 48*y3[idx-1] 
            - 36*y3[idx-2] 
            + 16*y3[idx-3] 
            - 3*y3[idx-4] 
        ) / (12*dr)
        
    end
    
    # Concatenate arrays again and update
    du .= hcat(dy0dt, dy1dt, dy2dt, dy3dt)

end



# Define a callback that checks if any field in 'u' grows too large
function check_field_too_large(u, t, integrator)
    max_value_allowed = 10^12  # Define your threshold here
    return any(x -> x > max_value_allowed, u)  # True if any field exceeds threshold
end

# Define a terminating function
function terminate_if_large!(integrator)
    println("Terminating because one of the fields grew too large.")
    terminate!(integrator)
end

# Create a DiscreteCallback that triggers if the check_field_too_large returns true
field_size_callback = DiscreteCallback(
    check_field_too_large, # Condition to check field size
    terminate_if_large!    # Action to terminate the solver
);
