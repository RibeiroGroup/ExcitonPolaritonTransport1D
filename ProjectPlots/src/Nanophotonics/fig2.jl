using Makie
using HDF5
using LaTeXStrings

function fig2(;σxvals=[120, 480], ΩR=0.1, σMvals = [0.02, 0.04, 0.1, 0.2], show_std=false)
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 800))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1], ylabel=L"Escape Probability ($\chi$)", title=L"$\sigma_x = %$(σxvals[1])$ nm")
    ax2 = Axis(gd[2,1], ylabel=L"RMSD ($\mu$m)", xticks=0:5)
    ax3 = Axis(gd[1,2], title=L"$\sigma_x = %$(σxvals[2])$ nm")
    ax4 = Axis(gd[2,2], xticks=0:5)
    axs = [ax1 ax3; ax2 ax4]

    # Link axes so they have the same plotting range
    linkxaxes!(axs...)
    linkyaxes!(ax3, ax1)
    linkyaxes!(ax4, ax2)

    xlims!(ax4, -0.1, 5)
    ylims!(ax3, 0, 0.7)
    ylims!(ax4, 0, 5)

    # Hide redundant axis info
    hideydecorations!(ax3, grid=false)
    hideydecorations!(ax4, grid=false)

    hidexdecorations!(ax1, grid=false)
    hidexdecorations!(ax3, grid=false)

    rel_dis = ["$(round(σM/ΩR, digits=2))" for σM in σMvals]

    clrs = Makie.wong_colors()
    for (j, σM) in enumerate(σxvals)
        σx = σxvals[j]
        for i = eachindex(σMvals) 
            c = clrs[i]
            escp_over_time!(axs[1,j], ΩR=ΩR, σx=σx, σM=σMvals[i],color=c,label=rel_dis[i])
            rmsd_propagation!(axs[2,j], ΩR=ΩR, σx=σx, σM=σMvals[i], color=c, show_std=show_std)
        end
    end

    linkxaxes!(ax1, ax2)

    Label(gd[3,1:2], "Time (ps)")

    Legend(gd[4,1:2], ax1, L"\sigma_M/\Omega_R", orientation=:horizontal, merge=true, titleposition=:left)
    fig
end


function wvp_fig2(;σxvals=[120, 480], ΩR=0.1, σMvals = [0.02, 0.04, 0.1, 0.2], show_std=false)

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 800))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    axs = [Axis(gd[i,j]) for i = 1:length(σMvals), j = 1:2]
    linkxaxes!(axs...)
    for i = eachindex(σMvals)
        linkyaxes!(axs[i,2], axs[i,1])
    end
    xlims!(axs[end,2], 15, 35)

    axs[1,1].title = L"$\sigma_x$ = %$(σxvals[1]) nm"
    axs[1,2].title = L"$\sigma_x$ = %$(σxvals[2]) nm"

    rel_dis = ["$(round(σM/ΩR, digits=2))" for σM in σMvals]

    clrs = Makie.wong_colors()
    for j in 1:2
        σx = σxvals[j]
        for i = eachindex(σMvals) 
            c = clrs[i]
            σM = σMvals[i]
            get_wvp!(axs[i,j], ΩR=ΩR, σx=σx, σM=σM, t=5, normalize=true, color=c)
            text!(axs[i,j], 0.9, 0.9, text=L"\sigma_M / \Omega_R = %$(rel_dis[i])",align=(:right,:top), fontsize=20, space=:relative)
        end
    end

    Label(gd[length(σMvals)+1,1:2], "Distance μm")

    fig
end

function altfig2()
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 800))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    ax = Axis(gd[1,1], xlabel="Disorder (eV)", ylabel=L"v_0\;(\text{ps}^{-1})")
    ax2 = Axis(gd[2,1], xlabel="Rabi Splitting (eV)", ylabel=L"v_0\;(\text{ps}^{-1})")
    #axs = [Axis(gd[i,j]) for i = 1:3, j = 1:3]

    clrs = Makie.wong_colors()
    mkers = [:circle, :rect, :cross, :utriangle, :star5]
    Rvals = [0.05, 0.1, 0.2, 0.3]
    σxvals = [60.0, 120.0, 240.0, 360.0, 480.0]
    for i = 1:5
        c = clrs[i]
        m = mkers[i]
        σx = σxvals[i]
        v0_vs_dis!(ax, σx=σx, ΩR=0.1 ,color=c, marker = m, label="$σx")
        v0_vs_rabi!(ax2, σx=σx, marker=m, color=c, σM=0.04)
    end


    axislegend(ax)
    fig
end

function fig3()
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 800))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    ax = Axis(gd[1,1:2], xlabel=L"\text{Relative disorder}\;(\sigma_M/\Omega_R)", ylabel=L"v_0\;(\text{ps}^{-1})")
    v0_vs_reldis!(ax)
    xlims!(ax, 0, 2.05)

    axislegend(ax, [[MarkerElement(marker=m, color=:black) for m = (:circle, :rect, :cross, :utriangle, :star5)],
    [PolyElement(color=c, strokecolor=:transparent) for c in Makie.wong_colors()[1:4]]],
    [string.([60,120,240,360,480]), string.([0.05, 0.1, 0.2, 0.3])],
    [L"\sigma_x\;(\text{nm})", L"\Omega_R\;(\text{eV})"], orientation=:vertical, titleposition=:left, nbanks=2)

    vlines!(ax, [0.2, 0.4],  color=:gray, linestyle=:dash)

    ax2 = Axis(gd[2,1], ylabel=L"v_0\;(\text{ps}^{-1})")
    v0_vs_σx_fixed_reldis!(ax2; reldis=0.2)
    text!(ax2, 0.5, 0.95, text=L"\sigma_M/\Omega_R = 0.2", align=(:center, :top), space=:relative)

    ax3 = Axis(gd[2,2])
    v0_vs_σx_fixed_reldis!(ax3; reldis=0.4)
    text!(ax3, 0.5, 0.95, text=L"\sigma_M/\Omega_R = 0.4", align=(:center, :top), space=:relative)

    linkyaxes!(ax2, ax3) 

    Label(gd[3,1:2], L"Initial wave packet spread $\sigma_M$ (nm)")
    fig
end
