#   Setup: Prepare initial conditions
#   –––––––––––––––––––––––––––––––––



function add_massless_scattering_ID!(x, y0, y1, y2, y3, pGrid, pModel, pInit)
    # unpack parameters
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    (
        x0phi, x0chi, # location parameters
        a0phi, a0chi, # amplitude parameters
        p0phi, p0chi, # some other parameter (here used as width)
        offsetphi, offsetchi, # parameter to control potential offset (TODO) 
        aStochastic, mink, maxk, stableRandomSeed # preparation for stochastic parameters (TODO)
    ) = pInit
    
    # set widths
    w0phi = p0phi
    w0chi = p0chi
        
    for n = NboundaryPadding + 2:Nx + NboundaryPadding
        # note that we ommit adding to the (singular) central point
        # add the Gaussian wave packet
        y0[n] += a0phi * abs(x[n])^(-(Ndim-1)/2) * exp(- (x[n] - x0phi)^2 / (2 * w0phi^2))
        y1[n] += - a0phi * abs(x[n])^(-(Ndim-1)/2) * (x[n] - x0phi)/w0phi^2 * exp(- (x[n] - x0phi)^2 / (2 * w0phi^2))
        y2[n] += a0chi * abs(x[n])^(-(Ndim-1)/2) * exp(- (x[n] - x0chi)^2 / (2 * w0chi^2))
        y3[n] += - a0chi * abs(x[n])^(-(Ndim-1)/2) * (x[n] - x0chi)/w0chi^2 * exp(- (x[n] - x0chi)^2 / (2 * w0chi^2))
    end
    
    return hcat(y0, y1, y2, y3)

end




function add_offset!(x, y0, y1, y2, y3, pGrid, pModel, pInit)
    # unpack parameters
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    (
        x0phi, x0chi, # location parameters
        a0phi, a0chi, # amplitude parameters
        p0phi, p0chi, # some other parameter (here used as width)
        offsetphi, offsetchi, # parameter to control potential offset (TODO) 
        aStochastic, mink, maxk, stableRandomSeed # preparation for stochastic parameters (TODO)
    ) = pInit
        
    for n = 1:Nx + 2*NboundaryPadding
        
        # add the offset
        
        y0[n] += offsetphi
        
        y2[n] += offsetchi
    
    end
    
    return hcat(y0, y1, y2, y3)

end
