#   Setup: Prepare initial conditions
#   –––––––––––––––––––––––––––––––––


include("./add_initial_data.jl")


function initial_data(pGrid, pModel, pInit)
    
    # unpack parameters
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # define a grid for the x variable
    # (note that this includes boundary points)
    x = range(- NboundaryPadding*dx, step=dx, length=(Nx + 2 * NboundaryPadding))
    
    # initialise the field variables
    y0 = zeros((Nx + 2 * NboundaryPadding))
    y1 = zeros((Nx + 2 * NboundaryPadding))
    y2 = zeros((Nx + 2 * NboundaryPadding))
    y3 = zeros((Nx + 2 * NboundaryPadding))
    
    # add the desired ID (other types of ID can be defined in add_initial_data.jl)
    add_offset!(x, y0, y1, y2, y3, pGrid, pModel, pInit)
    add_massless_scattering_ID!(x, y0, y1, y2, y3, pGrid, pModel, pInit)
    
    # update boundary points (TODO)
    
    # return
    fields = hcat(y0, y1, y2, y3)    
    return hcat(y0, y1, y2, y3)

end

