using DifferentialEquations # for the actual time evolution
using OrdinaryDiffEq # for ODEs
using Plots # for plotting
using Base.Threads # for parallelization
using StaticArrays # somehow needed to use multiple variables in DifferentialEquations.jl

using Plots, LaTeXStrings, Colors
using Plots.PlotMeasures
using LinearAlgebra

using Random, Distributions

using FFTW # discrete Fourier transform

using JLD2 # for file saving

level = "../../../"

include(joinpath(level, "src/4th-order-FD-stencils.jl"));
include(joinpath(level, "src/evolution_noabs.jl"));
include(joinpath(level, "src/hamiltonian_noabs.jl"));
include(joinpath(level, "src/initial_data.jl"));
include(joinpath(level, "src/visualisation_L=4.jl"));


include("../visualisation_for_paper.jl");

#   evolution
#   –––––––––

function artisan_evolution_at_resolution(Nx, stableRandomSeed, pModel, pInit)
    # unpack model parameters
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # define the domain
    Lx = 4.7  # Physical length of the domain in x direction
    dx = Lx/Nx  # Physical grid spacing
    NboundaryPadding = 4 # Number of boundary points
    
    # set the grid parameters at the given resolution
    pGrid = (dx, Nx, NboundaryPadding);

    # set the time span
    tspan = (0, 3.6);

    # set the evolution method
    time_integration_method = RK4();
    
    # generate initial conditions
    u0 = initial_data(pGrid, pModel, pInit);
    
    # TODO: implement export of a plot of the initial data
        
    # set the problem
    pHyper = (pGrid, pModel)
    prob = ODEProblem(finite_differenced_pde_with_bc!, u0, tspan, pHyper);

    sol = solve(
        prob, time_integration_method, 
        saveat = tspan[end]/(3*10^2),
        dt=dx/4, 
        adaptive = false, 
        dense=false, 
        maxiters=typemax(Int),
        #callback=field_size_callback
    );
        
    # obtain the hamiltonian
    hamiltonian = zeros(length(sol.u))
    hamphi = zeros(length(sol.u))
    hamchi = zeros(length(sol.u))
    current_radiationphi = 0
    current_radiationchi = 0
    radiationphi = zeros(length(sol.u))
    radiationchi = zeros(length(sol.u))
    for i = 1:length(sol.u)
        # radiation
        dtMonitor = i <= 1 ? sol.t[i] : sol.t[i] - sol.t[i-1]
        current_radiationphi += power_radiated_phi(sol.u[i], pHyper) * dtMonitor
        radiationphi[i] = current_radiationphi 
        current_radiationchi += power_radiated_chi(sol.u[i], pHyper) * dtMonitor
        radiationchi[i] = current_radiationchi 
        # bulk energies        
        hamiltonian[i] = nintegrate_simps(hamiltonian_density(sol.u[i], pHyper), dx)
        hamphi[i] = nintegrate_simps(hamiltonian_phi(sol.u[i], pHyper), dx)
        hamchi[i] = nintegrate_simps(hamiltonian_chi(sol.u[i], pHyper), dx)
        # add radiated energy to the total hamiltonian (for convergence)
        hamiltonian[i] = hamiltonian[i] + current_radiationphi + current_radiationchi
    end
    
    return (pGrid, sol, hamiltonian, hamphi, hamchi, radiationphi, radiationchi)
end

function evolution(stableRandomSeed, amp)
    
    Random.seed!(stableRandomSeed) 
    print("persistent random seed: ", stableRandomSeed, "\n")    
    print("current characteristic amplitude: ", amp, "\n")
    
    # parameters of the model
    mphi2, mchi2 = 0, 0   # masses
    lambdaCross = 1
    lambdaSelf = 1
    sigma = - 1   # ghostly boolean parameter (+1: no ghost; -1: ghost) 
    Ndim = 3   # number of spatial dimensions
    lV, mV, nV = 6, 2, 2   # positive integer exponents in the potentials
    
    # parameters of the initial data
#     random_sign_phi = 1 - 2*rand(0:1)
#     random_sign_chi = 1 - 2*rand(0:1)
    x0phi, x0chi = 1, 1   # location parameters
    a0phi = + amp
    a0chi = + amp
#     a0phi = random_sign_phi * amp * rand()  # amplitude parameters
#     a0chi = random_sign_chi * amp * rand()  # amplitude parameters
    p0phi, p0chi = 0.1, 0.1   # some other parameter (here used as width)
    offsetphi, offsetchi = 0, 0   # parameter to control potential offset (TODO) 
    aStochastic = 0 
    mink, maxk, = 0, 0  

    # set the combined set of parameters 
    pModel = (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV);
    pInit = (
        x0phi, x0chi, # location parameters
        a0phi, a0chi, # amplitude parameters
        p0phi, p0chi, # some other parameter (here used as width)
        offsetphi, offsetchi, # parameter to control potential offset (TODO) 
        aStochastic, mink, maxk, stableRandomSeed # preparation for stochastic parameters (TODO)
    );
     
    # set some tables to store output
    resTab = [2^i for i in 11:14]
    pGridTab = []
    solTab = []
    hamiltonianTab = []
    hamPhiTab = []
    hamChiTab = []
    radiationPhiTab = []
    radiationChiTab = []
    
    #############################
    # evolution
    #############################

    # run evolution
    for res in resTab
        print("current resolution: ", res, " ... \n")
        # run the evolution
        @time (pGrid, sol, hamiltonian, hamPhi, hamChi, radiationPhi, radiationChi) = artisan_evolution_at_resolution(
            res, stableRandomSeed, pModel, pInit
        )
        print("... terminated", "\n")
        
        # append the results
        push!(pGridTab, pGrid)
        push!(solTab, sol)
        push!(hamiltonianTab, hamiltonian)
        push!(hamPhiTab, hamPhi)
        push!(hamChiTab, hamChi)
        push!(radiationPhiTab, radiationPhi)
        push!(radiationChiTab, radiationChi)
    end
    

    
    #############################
    # SAVE OUTPUT DATA
    #############################
    
#     dir_path = "dat"
#     if !isdir(dir_path)
#         mkpath(dir_path)
#     end
    
#     timesteps = solTab[end].t
    
#     @save joinpath(pwd(), dir_path, string(amp,".jld2")) timesteps hamiltonianTab, hamPhiTab hamChiTab radiationPhiTab radiationChiTab
    
#     print("Finished output.", "\n")
    
    
    #############################
    # EXTRACTED ENERGY RATIOS
    #############################
    
    energyPhi = hamPhiTab[end] + radiationPhiTab[end]
    energyChi = hamChiTab[end] + radiationChiTab[end]
    initialHamiltonian = hamiltonianTab[end][1]
    
    ratioPhi = abs(energyPhi[1] - energyPhi[end])/(abs(energyPhi[1]) + abs(energyChi[1]))
    ratioChi = abs(energyChi[1] - energyChi[end])/(abs(energyPhi[1]) + abs(energyChi[1]))
    
    print("Extracted kinetic energy ratio for phi: ", ratioPhi, "\n")
    print("Extracted kinetic energy ratio for chi: ", ratioChi, "\n")
    
    
    #############################
    # PLOTTING
    #############################
    
    
    dir_path = string("plots/",stableRandomSeed,"/",amp)

    # create the directory if it does not yet exist
    if !isdir(dir_path)
        print("Output plot directory does not exist. Creating it ...\n")
        mkpath(dir_path)
    else
        print("Output plot directory already exists.\n")
    end
    
    # plot and determine convergence 
    loss_of_convergence_time = save_convergence_plots(
        resTab, pGridTab, solTab, hamiltonianTab, 
        dir_path
    )
    if loss_of_convergence_time >= solTab[end].t[end]
        print("Convergence kept at all times.\n")
    else
        print("Convergence lost at time t=",loss_of_convergence_time,"\n")
    end
    loss_of_convergence_time = solTab[end].t[end]
    
    # call further plotting routines
    
    # plot energy components
    save_energies_plot(
        resTab, pGridTab, solTab, 
        hamiltonianTab, hamPhiTab, hamChiTab, radiationPhiTab, radiationChiTab,
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )
#     save_difference_in_energies_plot(
#         resTab, pGridTab, solTab, 
#         hamiltonianTab, hamPhiTab, hamChiTab,
#         dir_path,
#         loss_of_convergence_time=loss_of_convergence_time
#     )
    
    # plot density plots of the fields
    save_density_plots_for_paper(
        solTab[end], pGridTab[end], pModel, pInit,
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )
    
    # determine the index of convergence loss
    loss_of_convergence_index = findfirst(t -> t > loss_of_convergence_time, solTab[end].t)
    if loss_of_convergence_index === nothing
        loss_of_convergence_index = length(solTab[end].t)
    end
#     # and then animate the fields
#     save_animation_abslog(
#         solTab[end][1:max(1,div(loss_of_convergence_index,1*10^2)):loss_of_convergence_index], 
#         pGridTab[end], pModel, pInit,
#         dir_path
#     );
    save_animation(
        solTab[end][1:max(1,div(loss_of_convergence_index,6*10^2)):loss_of_convergence_index], 
        pGridTab[end], pModel, pInit,
        dir_path
    );
    
    print("Finished plotting.", "\n")
    
    return (ratioPhi, ratioChi)
end

#   main()
#   ––––––

function main()
    
    stableRandomSeed = 0;#rand(1:10^7)
    
#     for log2amp in -5:1:5
        
#         amp = 2.0^log2amp
    
    for amp in 0.2:0.2:16.
    
        #############################
        # RUN EVOLUTION
        #############################
        
        (ratioPhi, ratioChi) = evolution(stableRandomSeed, amp)
        
        #############################
        # SAVE OUTPUT DATA
        #############################

        dir_path = "dat"
        if !isdir(dir_path)
            mkpath(dir_path)
        end

        filename = joinpath(
            pwd(),
            dir_path,
            "results.jld2"
        )
        case_id = string("seed=", stableRandomSeed, "A=", amp)

        jldopen(filename, "a+") do file
            # delete if current case already exists
            # TODO: implement checks before (for instance on whether precision of existing case is higher than current)
            if haskey(file, case_id)
                delete!(file, case_id)  # delete existing dataset
            end
            # then write case into jld2
            file[case_id] = (
                stableRandomSeed = stableRandomSeed,
                amp = amp,
                ratioPhi = ratioPhi,
                ratioChi = ratioChi
            )
        end
    end
    
end

main()

#   output
#   ––––––

using NBInclude
nbexport("main.jl", "main.ipynb")