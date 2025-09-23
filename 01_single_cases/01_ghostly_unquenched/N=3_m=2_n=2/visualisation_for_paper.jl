phicolor = colorant"#c51b8a";
chicolor = colorant"#00FFFF";




function save_density_plots_for_paper(sol, pGrid, pModel, pInit, dir_path; loss_of_convergence_time=Inf)
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
    
    # Clip space
    mask = space .<= 2.1
    space = space[mask]
    phiField = phiField[:, mask]   # keep all time, only restricted space
    chiField = chiField[:, mask]   # keep all time, only restricted space
    
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
    
    # generate plot with both fields (and without color bars)
    
    densityplot = plot(
        layout=(1), size = (600,600), framestyle=:box, dpi = 400, 
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{fields}\\;\\;\\phi\\;\\;\\textrm{and}\\;\\;\\chi"),
        xlabel=L"$\textrm{space}\;\;r/R$",
        xticks = 0:0.5:2, 
        ylabel=L"$\textrm{time}\;\;t/R$",
        left_margin = 40px,
        right_margin = 40px,
        bottom_margin = 0px,
        top_margin = 10px
    )
    heatmap!(
        densityplot,
        space, time, phiField,
        color=cgradphiAlpha,
        cbar=false,
        clim=(-rangePhi,rangePhi)
    )
    heatmap!(
        densityplot,
        space, time, chiField,
        color=cgradchiAlpha,
        cbar=false,
        clim=(-rangeChi,rangeChi)
    )
    
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
        xticks = 0:0.5:2, 
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
    
    final_plot_2 = plot(
        densityplotphilegend, densityplot, densityplotchilegend,
        layout=@layout[a{0.2w} b{0.6w} c{0.2w}],
        #spacing=0px,
        size = (800,600)
    )

    # Save combined figures    
    savefig(
        final_plot_1, 
        joinpath(pwd(), dir_path, "evolution_density_for_paper.png")
    )
    savefig(
        final_plot_2, 
        joinpath(pwd(), dir_path, "evolution_density_overlayed_for_paper.png")
    )
    
end



