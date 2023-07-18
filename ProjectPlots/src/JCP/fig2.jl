using Makie
using HDF5
using LaTeXStrings

function fig2()
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
