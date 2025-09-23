phicolor = colorant"#c51b8a";
chicolor = colorant"#00FFFF";




# Function to plot the solution and save the GIF
function save_animation(sol, pGrid, pModel, pInit, dir_path; fps=16)
    
    # unpack parameters
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # some definitions
    timesteps = sol.t
    N = size(sol[1], 1)
    x = range(0, length=N)*dx
    timeIter = size(sol)[end]
    
    rangePrev = 2.;

    animation = @animate for i in 1:timeIter
        
        y0, y1, y2, y3 = [
            sol[:, i][1:N], 
            sol[:, i][N+1:2N], 
            sol[:, i][2N+1:3N], 
            sol[:, i][3N+1:4N]
        ]
        
        rangePhi = maximum(abs.(y0));
        rangeChi = maximum(abs.(y2));
        range = maximum([rangePhi,rangeChi,rangePrev]);
#         rangePrev = range;
        
        t = round(timesteps[i],digits=1)
        
        plot(
            layout=(1), framestyle=:box, dpi = 300, 
            legend=:topright,
            guidefont = "Computer Modern", tickfont = "Computer Modern",
            xguidefontsize = 20, yguidefontsize = 20, 
            legendfontsize = 20, titlefontsize = 20, tickfontsize = 20,
            xlims=(0, N*dx-NboundaryPadding*dx), 
            xticks = 0:1:4, 
            ylims=(-range*1.1, range*1.1), 
            title=latexstring("\\textrm{physical\\;\\;time:\\;\\;} t =", t),
            xlabel=L"$\textrm{space}\;\;r/R$",
            ylabel=L"$\textrm{field\;\;value}$",
            right_margin = 20px,
        )
        
        plot!(
            x[NboundaryPadding+1:N-NboundaryPadding], y0[NboundaryPadding+1:N-NboundaryPadding], 
            label=L"$\phi(x)$",
            linewidth=2, linecolor=phicolor
        )
        plot!(
            x[NboundaryPadding+1:N-NboundaryPadding], y2[NboundaryPadding+1:N-NboundaryPadding], 
            label=L"$\chi(x)$",
            linewidth=2, linecolor=chicolor
        )
    end
    
#     filename = join(["animation_Nx=", string(Nx), ".gif"])
#     gif(animation, joinpath(pwd(), dir_path, filename), fps=fps)
    filename = join(["animation_Nx=", string(Nx), ".mp4"])
    mp4(animation, joinpath(pwd(), dir_path, filename), fps=fps)
end





# Function to plot the solution and save the GIF
function save_animation_abslog(sol, pGrid, pModel, pInit, dir_path; fps=16)
    
    # unpack parameters
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # some definitions
    timesteps = sol.t
    N = size(sol[1], 1)
    x = range(0, length=N)*dx
    timeIter = size(sol)[end]
    
    rangePrev = 2.;

    animation = @animate for i in 1:timeIter
        
        y0, y1, y2, y3 = [
            sol[:, i][1:N], 
            sol[:, i][N+1:2N], 
            sol[:, i][2N+1:3N], 
            sol[:, i][3N+1:4N]
        ]
        
        rangePhi = maximum(abs.(y0));
        rangeChi = maximum(abs.(y2));
        range = maximum([rangePhi,rangeChi,rangePrev]);
        rangePrev = range;
        
        t = round(timesteps[i],digits=1)
        
        plot(
            layout=(1), framestyle=:box, dpi = 300, 
            legend=:topright,
            guidefont = "Computer Modern", tickfont = "Computer Modern",
            xguidefontsize = 20, yguidefontsize = 20, 
            legendfontsize = 20, titlefontsize = 20, tickfontsize = 20,
            xlims=(0, N*dx-NboundaryPadding*dx), 
            xticks = 0:1:4, 
            ylims=(1e-3, range*1.1),
            yscale=:log10,
            title=latexstring("\\textrm{physical\\;\\;time:\\;\\;} t =", t),
            xlabel=L"$\textrm{space}\;\;r/R$",
            ylabel=L"$\textrm{field\;\;value}$",
        )
        
        plot!(
            x[NboundaryPadding+1:N-NboundaryPadding], abs.(y0[NboundaryPadding+1:N-NboundaryPadding]), 
            label=L"$\phi(x)$",
            linewidth=2, linecolor=phicolor
        )
        plot!(
            x[NboundaryPadding+1:N-NboundaryPadding], abs.(y2[NboundaryPadding+1:N-NboundaryPadding]), 
            label=L"$\chi(x)$",
            linewidth=2, linecolor=chicolor
        )
    end
    
#     filename = join(["animation_abslog_Nx=", string(Nx), ".gif"])
#     gif(animation, joinpath(pwd(), dir_path, filename), fps=fps)
    filename = join(["animation_abslog_Nx=", string(Nx), ".mp4"])
    mp4(animation, joinpath(pwd(), dir_path, filename), fps=fps)
end





# Function to plot the solution and save the GIF
function save_snaps(sol, pGrid, pModel, pInit; snap_intervals=50, yrangeVal=1, dir_path="plots")
    
    # unpack parameters
    (dx, Nx, NboundaryPadding) = pGrid
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Ndim, lV, mV, nV) = pModel
    
    # some definitions
    timesteps = sol.t
    N = size(sol[1], 1)
    x = range(0, length=N)*dx
    timeIter = size(sol)[end]
    
    rangePrev = 2.;

    snaps_plot = plot(
        layout=(1), framestyle=:box, dpi = 400, size = (400,300), 
        legend=:bottom,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 12, yguidefontsize = 12, 
        legendfontsize = 6, titlefontsize = 12, tickfontsize = 12,
        xlims=(-NboundaryPadding*dx, N*dx-NboundaryPadding*dx), 
        xticks = 0:1:4, xmirror = true,
        #ylims=(-yrangeVal, yrangeVal), 
        #title=latexstring("\\textrm{physical\\;\\;time: }", t),
        xlabel=L"$\textrm{space}\;\;r/R$", 
        ylabel=L"$\textrm{field\;\;values}$",
    )
    
    for i in 1:snap_intervals:timeIter
        
        y0, y1, y2, y3 = [
            sol[:, i][1:N], 
            sol[:, i][N+1:2N], 
            sol[:, i][2N+1:3N], 
            sol[:, i][3N+1:4N]
        ]
                
        rangePhi = maximum(abs.(y0));
        rangeChi = maximum(abs.(y2));
        range = maximum([rangePhi,rangeChi,rangePrev]);
        rangePrev = range;
        
        t = round(timesteps[i],digits=1)
           
        plot!(
            x .- (NboundaryPadding*dx), y0, 
            label=latexstring("\\phi\\;\\textrm{at}\\;t/R=",t),
            linewidth=2, 
            linecolor=RGBA(phicolor.r, phicolor.g, phicolor.b, (1-i/timeIter))
        )
        plot!(
            x .- (NboundaryPadding*dx), y2, 
            label=latexstring("\\chi\\;\\textrm{at}\\;t/R=",t),
            linewidth=2, 
            linecolor=RGBA(chicolor.r, chicolor.g, chicolor.b, (1-i/timeIter)),
            linestyle=:dot
        )
    end
    
    savefig(snaps_plot, joinpath(pwd(), dir_path, "snaps.pdf"));
end










# saves the convergence plot in dir_path
# also returns the timestep at which convergence is lost
# at any point during evolution
function save_convergence_plots(resTab, pTab, solTab, hamiltonianTab, dir_path)
    
    order = 0;
    convergence_maintained = true;
    loss_of_convergence_time = solTab[end].t[end]
    loss_of_convergence_step = length(solTab[end].t)
    
    # CONVERGENCE
    
    #set up the plot environment
    convergenceplot = plot(
        layout=(1), framestyle=:box, dpi = 400, size = (400,300), 
        legend=:bottom,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 12, yguidefontsize = 12, 
        legendfontsize = 8, titlefontsize = 12, tickfontsize = 12,
        #xlims=(0,2),
        yscale=:log, 
        ylims=(2e-15, 2e2), 
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlabel=L"$\textrm{t}$", 
        ylabel=if order==0
            latexstring("|H-H_0|")
        else
            latexstring("|H-H_0|/dx^",order)
        end
    )

    # convergence plot 
    for i = 1:length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        hamiltonian = hamiltonianTab[i]
        #print(hamiltonian[1],"\n")
        #dH = abs.((hamiltonian .- hamiltonian[1])/hamiltonian[1])
        dH = abs.(hamiltonian .- hamiltonian[1])
        #print(dH[end],"\n")
        (dx, Nx, NboundaryPadding) = p
        # append to plot
        plot!(
            convergenceplot, 
            sol.t, dH/dx^order,
            label=latexstring("N_x=", res),
            linewidth=2#, linecolor=colorant"#f03b20"
        )
    end
    
    # Save the plot to a file (e.g., PNG format)
    savefig(convergenceplot, joinpath(pwd(), dir_path,"convergence-hamiltonian.png"));# SELF-CONVERGENCE
    
    #set up the plot environment
    selfconvergenceplot = plot(
        layout=(1), framestyle=:box, dpi = 400, size = (400,300), 
        legend=:bottom,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 12, yguidefontsize = 12, 
        legendfontsize = 8, titlefontsize = 12, tickfontsize = 12,
        #xlims=(0,2),
        #yscale=:log, 
        ylims=(0, 6), 
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlabel=L"$\textrm{t}$", 
        ylabel=L"$\textrm{log}\frac{||(u_i-u_{i+1})||}{||(u_{i+1}-u_{i+2})||}$"
    )

    # convergence plot 
    for i = 1:length(resTab)-2
        # unpack parameters
        p = pTab[i]
        (dx, Nx, NboundaryPadding) = p
        res = resTab[i]
        sol0 = solTab[i]
        sol1 = solTab[i+1]
        sol2 = solTab[i+2]

        selfconv = []
        for j = 1:min(length(sol0),length(sol1),length(sol2))
            u0 = sol0.u[j][1+NboundaryPadding:2^(1-1):end-NboundaryPadding,:]
            u1 = sol1.u[j][1+NboundaryPadding:2^(2-1):end-NboundaryPadding,:]
            u2 = sol2.u[j][1+NboundaryPadding:2^(3-1):end-NboundaryPadding,:]
                        
            currentconv = log2(norm(u0 - u1)/norm(u1 - u2))
            push!(selfconv, currentconv)
            # for the highest set of resolutions
            if convergence_maintained == true && i==length(resTab)-2 && j>min(length(sol0),length(sol1),length(sol2))/2^3
                if currentconv < 3 || currentconv > 5
                    convergence_maintained = false
                    loss_of_convergence_time = sol2.t[j]
                    loss_of_convergence_step = j
                end
            end
        end

        #print(selfconv,"\n")

        #append to plot
        plot!(
            selfconvergenceplot, 
            sol0.t[1:min(length(sol0),length(sol1),length(sol2))], selfconv,
            label=latexstring("N_{x,i}=", res),
            linewidth=2#, linecolor=colorant"#f03b20"
        )
        
        # for the highest set of resolutions
        if i==length(resTab)-2
            # check that convergence rate remains between 1 and 3 throughout
            if minimum(selfconv) < 1 || maximum(selfconv) > 3
                convergence_maintained = false;
            end
        end
    end

    savefig(selfconvergenceplot, joinpath(pwd(), dir_path, "self-convergence.png"));

    return min(loss_of_convergence_time, solTab[end].t[end])
    
end







function save_energies_plot(resTab, pTab, solTab, hamiltonianTab, hamPhiTab, hamChiTab, radiationPhiTab, radiationChiTab, dir_path; loss_of_convergence_time=Inf)
    
    # determine limit values for plots
    if loss_of_convergence_time isa Number
        tlim = (solTab[end].t[2], loss_of_convergence_time)
    else
        tlim = (solTab[end].t[2], solTab[end].t[end])
    end
    loss_of_convergence_time = min(solTab[end].t[end],loss_of_convergence_time)
    dt = solTab[end].t[2] - solTab[end].t[1];
    Nt = Int(min(length(solTab[end].t), div(loss_of_convergence_time,dt)));
    ylimSetAuto = (
        maximum([minimum([
            minimum(abs.(hamiltonianTab[end][1:Nt])),
            minimum(abs.(hamPhiTab[end][1:Nt])),
            minimum(abs.(hamChiTab[end][1:Nt])),
            minimum(abs.(hamiltonianTab[end][1:Nt] - hamPhiTab[end][1:Nt] - hamChiTab[end][1:Nt]))     
        ]),10^-11]),
        minimum([maximum([
            maximum(abs.(hamiltonianTab[end][1:Nt])),
            maximum(abs.(hamPhiTab[end][1:Nt])),
            maximum(abs.(hamChiTab[end][1:Nt])),
            maximum(abs.(hamiltonianTab[end][1:Nt] - hamPhiTab[end][1:Nt] - hamChiTab[end][1:Nt]))
        ]),10^11])
    );
    
    ## PLOT ONCE WITHOUT ABSOLUTE VALUE
    
    #set up the plot environment
    energiesPlot = plot(
        layout=(1), framestyle=:box, dpi = 400, size = (400,400), 
        legend=:bottom,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 12, yguidefontsize = 12, 
        legendfontsize = 8, titlefontsize = 12, tickfontsize = 12,
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlim=tlim,
        #ylims=ylimSetAuto,
        xlabel=L"$t/R$", 
        ylabel=L"$\textrm{component}\;\;\textrm{energies}\;\;\textrm{per}\;\;\textrm{unit}\;\;\textrm{sphere}$"
        #ylabel=L"$\textrm{kinetic}\;\;\textrm{energies:}\;\;T_\phi\;\;\textrm{and}\;\;T_\chi$"
    )

    for i = length(resTab):length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        radiationPhi = radiationPhiTab[i]
        radiationChi = radiationChiTab[i]
        radiation = radiationPhi + radiationChi
        #print(dH[end],"\n")
        (dx, Nx, NboundaryPadding) = p
        # append to plot
        plot!(
            energiesPlot, 
            sol.t, hamiltonian,
            label=latexstring("H^{N=",res,"} + E_\\mathrm{rad}^{N=",res,"}"),
            linewidth=2, linecolor="black"
        )
        plot!(
            energiesPlot, 
            sol.t, hamPhi + radiationPhi,
            label=latexstring("H_\\phi^{N=",res,"} + E_{\\mathrm{rad},\\phi}^{N=",res,"}"),
            linewidth=2, linecolor=phicolor
        )
        plot!(
            energiesPlot, 
            sol.t, hamChi + radiationChi,
            label=latexstring("H_\\chi^{N=",res,"} + E_{\\mathrm{rad},\\chi}^{N=",res,"}"),
            linewidth=2, linecolor=chicolor
        )
        plot!(
            energiesPlot, 
            sol.t, hamiltonian-hamPhi-hamChi-radiation,
            label=latexstring("V_\\mathrm{int}^{N=",res,"}"),
            linewidth=2, linecolor="orange"
        )
        plot!(
            energiesPlot, 
            sol.t, radiation,
            label=latexstring("E_\\mathrm{rad}^{N=",res,"}"),
            linewidth=2, linecolor="orange", linestyle=:dot
        )
    end

    for i = 1:(length(resTab)-1)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        radiationPhi = radiationPhiTab[i]
        radiationChi = radiationChiTab[i]
        radiation = radiationPhi + radiationChi
        #print(dH[end],"\n")
        (dx, Nx, NboundaryPadding) = p
        # append to plot
        plot!(
            energiesPlot, 
            sol.t, hamiltonian,
            label=false,
            linewidth=0.1/(length(resTab) - i), linecolor=colorant"#000000",
        )
        plot!(
            energiesPlot, 
            sol.t, hamPhi + radiationPhi,
            label=false,
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, hamChi + radiationChi,
            label=false,
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, hamiltonian-hamPhi-hamChi-radiation,
            label=false,
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, radiation,
            label=false,
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dot
        )
    end

    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "energies.png"));
    
    ## RE-PLOT WITH ABSOLUTE VALUE FOR LOG-PLOTS
    
    #set up the plot environment
    energiesPlot = plot(
        layout=(1), framestyle=:box, dpi = 400, size = (400,400), 
        legend=:bottom,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 12, yguidefontsize = 12, 
        legendfontsize = 8, titlefontsize = 12, tickfontsize = 12,
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlim=tlim,
        ylims=ylimSetAuto,
        xlabel=L"$t/R$", 
        ylabel=L"$\textrm{component}\;\;\textrm{energies}$"
        #ylabel=L"$\textrm{kinetic}\;\;\textrm{energies:}\;\;T_\phi\;\;\textrm{and}\;\;T_\chi$"
    )

    for i = length(resTab):length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        radiationPhi = radiationPhiTab[i]
        radiationChi = radiationChiTab[i]
        radiation = radiationPhi + radiationChi
        #print(dH[end],"\n")
        (dx, Nx, NboundaryPadding) = p
        # append to plot
        plot!(
            energiesPlot, 
            sol.t, abs.(hamiltonian),
            label=latexstring("|H^{N=",res,"} + E_\\mathrm{rad}^{N=",res,"}|"),
            linewidth=2, linecolor="black"
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamPhi + radiationPhi),
            label=latexstring("|H_\\phi^{N=",res,"} + E_{\\mathrm{rad},\\phi}^{N=",res,"}|"),
            linewidth=2, linecolor=phicolor
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamChi + radiationChi),
            label=latexstring("|H_\\chi^{N=",res,"} + E_{\\mathrm{rad},\\chi}^{N=",res,"}|"),
            linewidth=2, linecolor=chicolor
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamiltonian-hamPhi-hamChi-radiation),
            label=latexstring("|V_\\mathrm{int}^{N=",res,"}|"),
            linewidth=2, linecolor="orange"
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(radiation),
            label=latexstring("|E_\\mathrm{rad}^{N=",res,"}|"),
            linewidth=2, linecolor="orange", linestyle=:dot
        )
    end

    for i = 1:(length(resTab)-1)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        radiationPhi = radiationPhiTab[i]
        radiationChi = radiationChiTab[i]
        radiation = radiationPhi + radiationChi
        #print(dH[end],"\n")
        (dx, Nx, NboundaryPadding) = p
        # append to plot
        plot!(
            energiesPlot, 
            sol.t, abs.(hamiltonian),
            label=false,
            linewidth=0.1/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamPhi + radiationPhi),
            label=false,
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamChi + radiationChi),
            label=false,
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamiltonian-hamPhi-hamChi-radiation),
            label=false,
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(radiation),
            label=false,
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dot
        )
    end
    
    energiesPlot = plot!(
        energiesPlot,
        yscale=:log, 
        ylims=ylimSetAuto,
        xlim=tlim
    )
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "energies_log.png"));
    
    energiesPlot = plot!(
        energiesPlot,
        yscale=:log, 
        ylims=ylimSetAuto,
        xscale=:log,
        xlim=tlim
    )
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "energies_loglog.png"));
    
end






function save_difference_in_energies_plot(resTab, pTab, solTab, hamiltonianTab, hamPhiTab, hamChiTab, dir_path; loss_of_convergence_time=Inf)
    
    # determine limit values for plots
    if loss_of_convergence_time isa Number
        tlim = (solTab[end].t[2], loss_of_convergence_time)
    else
        tlim = (solTab[end].t[2], solTab[end].t[end])
    end
    loss_of_convergence_time = min(solTab[end].t[end],loss_of_convergence_time)  
    dt = solTab[end].t[2] - solTab[end].t[1];
    Nt = Int(min(length(solTab[end].t), div(loss_of_convergence_time,dt)));
    ylimSetAuto = (
        maximum([minimum([
            minimum(abs.(hamiltonianTab[end][2:Nt] .- hamiltonianTab[end][1])),
            minimum(abs.(hamPhiTab[end][2:Nt] .- hamPhiTab[end][1])),
            minimum(abs.(hamChiTab[end][2:Nt] .- hamChiTab[end][1]))    
        ]),10^-11]),
        minimum([maximum([
            maximum(abs.(hamiltonianTab[end][2:Nt] .- hamiltonianTab[end][1])),
            maximum(abs.(hamPhiTab[end][2:Nt] .- hamPhiTab[end][1])),
            maximum(abs.(hamChiTab[end][2:Nt] .- hamChiTab[end][1]))  
        ]),10^11])
    );
   
    ## PLOT ONCE WITHOUT ABSOLUTE VALUE
    
    #set up the plot environment
    energiesPlot = plot(
        layout=(1), framestyle=:box, dpi = 400, size = (400,400), 
        legend=:bottom,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 12, yguidefontsize = 12, 
        legendfontsize = 8, titlefontsize = 12, tickfontsize = 12,
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlim=tlim,
        ylims=ylimSetAuto,
        xlabel=L"$t/R$", 
        ylabel=L"$\textrm{component}\;\;\textrm{energies}$"
        #ylabel=L"$\textrm{kinetic}\;\;\textrm{energies:}\;\;T_\phi\;\;\textrm{and}\;\;T_\chi$"
    )

    for i = length(resTab):length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        (dx, Nx, NboundaryPadding) = p
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        #build difference
        hamiltonian = hamiltonian.-hamiltonian[1]
        hamPhi = hamPhi.-hamPhi[1]
        hamChi = hamChi.-hamChi[1]
        # append to plot
        plot!(
            energiesPlot, 
            sol.t, hamiltonian,
            label=latexstring("H^{N=",res,"} - H^{0}"),
            linewidth=2, linecolor="black", linestyle=:dash
        )
        plot!(
            energiesPlot, 
            sol.t, hamPhi,
            label=latexstring("H_\\phi^{N=",res,"} - H_\\phi^{0}"),
            linewidth=2, linecolor=phicolor
        )
        plot!(
            energiesPlot, 
            sol.t, hamChi,
            label=latexstring("H_\\chi^{N=",res,"} - H_\\chi^{0}"),
            linewidth=2, linecolor=chicolor
        )
        plot!(
            energiesPlot, 
            sol.t, hamiltonian-hamPhi-hamChi,
            label=latexstring("V_\\mathrm{int}^{N=",res,"} - V_\\mathrm{int}^{0}"),
            linewidth=2, linecolor="orange", linestyle=:dot
        )
    end

    for i = 1:(length(resTab)-1)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        (dx, Nx, NboundaryPadding) = p
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        #build difference
        hamiltonian = hamiltonian.-hamiltonian[1]
        hamPhi = hamPhi.-hamPhi[1]
        hamChi = hamChi.-hamChi[1]
        # append to plot
        plot!(
            energiesPlot, 
            sol.t, hamiltonian,
            label=false,#latexstring("H^{N=",res,"} - H^{0}"),
            linewidth=0.1/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dash
        )
        plot!(
            energiesPlot, 
            sol.t, hamPhi,
            label=false,#latexstring("H_\\phi^{N=",res,"} - H_\\phi^{0}"),
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, hamChi,
            label=false,#latexstring("H_\\chi^{N=",res,"} - H_\\chi^{0}"),
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, hamiltonian-hamPhi-hamChi,
            label=false,#latexstring("V_\\mathrm{int}^{N=",res,"} - V_\\mathrm{int}^{0}"),
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dot
        )
    end
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "energies_difference.png"));
    
    ## RE-PLOT WITH ABSOLUTE VALUE FOR LOG-PLOTS
   
    #set up the plot environment
    energiesPlot = plot(
        layout=(1), framestyle=:box, dpi = 400, size = (400,400), 
        legend=:bottom,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 12, yguidefontsize = 12, 
        legendfontsize = 8, titlefontsize = 12, tickfontsize = 12,
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlim=tlim,
        ylims=ylimSetAuto,
        xlabel=L"$t/R$", 
        ylabel=L"$\textrm{component}\;\;\textrm{energies}$"
        #ylabel=L"$\textrm{kinetic}\;\;\textrm{energies:}\;\;T_\phi\;\;\textrm{and}\;\;T_\chi$"
    )

    for i = length(resTab):length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        (dx, Nx, NboundaryPadding) = p
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        #build difference
        hamiltonian = hamiltonian.-hamiltonian[1]
        hamPhi = hamPhi.-hamPhi[1]
        hamChi = hamChi.-hamChi[1]
        # append to plot
        plot!(
            energiesPlot, 
            sol.t, abs.(hamiltonian),
            label=latexstring("|H^{N=",res,"} - H^{0}|"),
            linewidth=2, linecolor="black", linestyle=:dash
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamPhi),
            label=latexstring("|H_\\phi^{N=",res,"} - H_\\phi^{0}|"),
            linewidth=2, linecolor=phicolor
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamChi),
            label=latexstring("|H_\\chi^{N=",res,"} - H_\\chi^{0}|"),
            linewidth=2, linecolor=chicolor
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamiltonian-hamPhi-hamChi),
            label=latexstring("|V_\\mathrm{int}^{N=",res,"} - V_\\mathrm{int}^{0}|"),
            linewidth=2, linecolor="orange", linestyle=:dot
        )
    end

    for i = 1:(length(resTab)-1)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        (dx, Nx, NboundaryPadding) = p
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        #build difference
        hamiltonian = hamiltonian.-hamiltonian[1]
        hamPhi = hamPhi.-hamPhi[1]
        hamChi = hamChi.-hamChi[1]
        # append to plot
        plot!(
            energiesPlot, 
            sol.t, abs.(hamiltonian),
            label=false,#latexstring("|H^{N=",res,"} - H^{0}|"),
            linewidth=0.1/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dash
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamPhi),
            label=false,#latexstring("|H_\\phi^{N=",res,"} - H_\\phi^{0}|"),
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamChi),
            label=false,#latexstring("|H_\\chi^{N=",res,"} - H_\\chi^{0}|"),
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamiltonian-hamPhi-hamChi),
            label=false,#latexstring("|V_\\mathrm{int}^{N=",res,"} - V_\\mathrm{int}^{0}|"),
            linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dot
        )
    end

    energiesPlot = plot!(
        energiesPlot,
        yscale=:log, 
        ylims=ylimSetAuto
    )
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "energies_difference_log.png"));
    
    energiesPlot = plot!(
        energiesPlot,
        yscale=:log, 
        ylims=ylimSetAuto,
        xscale=:log,
        xlim=tlim
    )
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "energies_difference_loglog.png"));
end





function save_density_plots(sol, pGrid, pModel, pInit, dir_path; loss_of_convergence_time=Inf)
    # determine the desired end time
    if !(loss_of_convergence_time isa Number)
        loss_of_convergence_time = sol.t[end]
    end
    loss_of_convergence_time = min(loss_of_convergence_time, sol.t[end])
    
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
    
    # define further parameters
    N = Nx + 2*NboundaryPadding;
    dt = sol.t[2] - sol.t[1];
    Nt = Int(min(length(sol.t), div(loss_of_convergence_time,dt)));

    time = range(0, stop=Nt*dt, length=Nt);
    space = range(0, stop=N*dx, length=N);
    phiField = [sol[:, i][j] for i = 1:Nt, j=1:N]; 
    chiField = [sol[:, i][j] for i = 1:Nt, j=2N+1:3N];
    phiField = [sol[:, i][j] - offsetphi for i = 1:Nt, j=1:N]; 
    chiField = [sol[:, i][j] - offsetchi for i = 1:Nt, j=2N+1:3N];
    
    # set density ledgend range
    rangePhi = maximum((maximum(abs.(phiField)),2))/2;
    rangeChi = maximum((maximum(abs.(chiField)),2))/2;
    #rangePhi = minimum((maximum(abs.(phiField[1:div(Nt, 2),:])),3.5));
    #rangeChi = minimum((maximum(abs.(chiField[1:div(Nt, 2),:])),3.5));
    
    spacingPhi = 10^(floor(log10(rangePhi)))
    rangePhiFloor = floor(rangePhi/spacingPhi)*spacingPhi
    spacingChi = 10^(floor(log10(rangeChi)))
    rangeChiFloor = floor(rangeChi/spacingChi)*spacingChi
    
    # define colors for gradients
    black = RGBA(0, 0, 0, 1)
    white = RGBA(1, 1, 1, 1)
    target_color_phi = RGBA(phicolor.r, phicolor.g, phicolor.b, 1)
    mid_color_phi_1 = RGBA(0.5 * phicolor.r, 0.5 * phicolor.g, 0.5 * phicolor.b, 1)
    mid_color_phi_2 = RGBA(0.5 + 0.5 * phicolor.r, 0.5 + 0.5 * phicolor.g, 0.5 + 0.5 * phicolor.b, 1)
    target_color_chi = RGBA(chicolor.r, chicolor.g, chicolor.b, 1)
    mid_color_chi_1 = RGBA(0.5 * chicolor.r, 0.5 * chicolor.g, 0.5 * chicolor.b, 1)
    mid_color_chi_2 = RGBA(0.5 + 0.5 * chicolor.r, 0.5 + 0.5 * chicolor.g, 0.5 + 0.5 * chicolor.b, 1)
    
    # define color gradients with alpha transparency
    cgradphiAlpha = cgrad([
        RGBA(0, 0, 0, 1),
        RGBA(phicolor.r, phicolor.g, phicolor.b, 0), 
        RGBA(phicolor.r, phicolor.g, phicolor.b, 1)
    ], -1:0.1:1)
    cgradchiAlpha = cgrad([
        RGBA(0, 0, 0, 1),
        RGBA(chicolor.r, chicolor.g, chicolor.b, 0), 
        RGBA(chicolor.r, chicolor.g, chicolor.b, 1)
    ], -1:0.1:1)
    
    # define non-transparent color gradients 
    # (workaround to avoid bad legend behaviour in the single field plots)
    # (for some reason julia's Plots.jl legends to not account for the alpha value)
    cgradphi = cgrad([
            black,
            mid_color_phi_1,
            white,
            mid_color_phi_2,
            target_color_phi
        ], -1:0.1:1)
    cgradchi = cgrad([
            black,
            mid_color_chi_1,
            white,
            mid_color_chi_2,
            target_color_chi
        ], -1:0.1:1)
    
    # generate plot with phi only
    
    # reset density ledend range
    #rangePhi = maximum((maximum(abs.(phiField)),0));
    #rangeChi = maximum((maximum(abs.(chiField)),0));
    
    densityplotphi = plot(
        layout=(1), size = (400,600), framestyle=:box, dpi = 400, 
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{field}\\;\\;\\phi"),
        xlabel=L"$\textrm{space}\;\;r/R$",
        xticks = 0:1:4, 
        xflip=true,
        ylabel=L"$\textrm{time}\;\;t/R$",
        left_margin = 60px,
        right_margin = -6px,
        bottom_margin = 40px,
        top_margin = 10px
    )
    heatmap!(
        densityplotphi,
        space, time, phiField,
        color=cgradphi,
        cbar=false,
        clim=(-rangePhi,rangePhi),
        colorbar_ticks=-rangePhiFloor:spacingPhi:rangePhiFloor
    )
    contour!(
        densityplotphi,
        space, time, phiField;
        levels=10,#Int(rangePhiFloor/spacingPhi),
        linewidth=0.1, color=:black
    )

    # generate plot with chi only
    
    densityplotchi = plot(
        layout=(1), size = (400,600), framestyle=:box, dpi = 400, 
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{field}\\;\\;\\chi"),
        xlabel=L"$\textrm{space}\;\;r/R$",
        xticks = 0.5:0.5:2, 
#         ylabel=L"$\textrm{time}\;\;t/R$",
        ylabel="",
        yticks=false,
        left_margin = -6px,
        right_margin = 60px,
        bottom_margin = 40px,
        top_margin = 10px
    )
    heatmap!(
        densityplotchi,
        space, time, chiField,
        color=cgradchi,
        cbar=false,
        clim=(-rangeChi,rangeChi),
        colorbar_ticks=-rangeChiFloor:spacingChi:rangeChiFloor
    )
    contour!(
        densityplotchi,
        space, time, chiField;
        levels=10,#Int(rangeChiFloor/spacingChi), 
        linewidth=0.1, color=:black
    )
        
    # generate legends
               
    densityplotphilegend = plot(
        densityplotphi,
        colorbar=true,
        framestyle=:none,
        xticks=false,
        yticks=false,
        xlims=(-2, -1),   # Move x out of visible range
        ylims=(-2, -1),   # Move y out of visible range
        size=(80, 600),
        top_margin = 10px
    )
    
    densityplotchilegend = plot(
        densityplotchi,
        colorbar=true,
        framestyle=:none,
        xticks=false,
        yticks=false,
        xlims=(-2, -1),   # Move x out of visible range
        ylims=(-2, -1),   # Move y out of visible range
        size=(80, 600),
        top_margin = 10px
    )
    
    # Combine plots
    final_plot_1 = plot(
        densityplotphilegend, densityplotphi, densityplotchi, densityplotchilegend,
        layout=@layout[a{0.1w} b{0.4w} c{0.4w} d{0.1w}],
        spacing=0px,
        size = (1300,600)
    )

    # Save combined figures    
    savefig(
        final_plot_1, 
        joinpath(pwd(), dir_path, "evolution_density.png")
    )
    
end



