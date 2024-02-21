function fig8()

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 400))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()

    ax1 = Axis(gd[1,1], xticks=[0, 0.005, 0.010], ylabel=L"$v_{\alpha q}^\text{eff}$ ($\mu$m$\cdot$ps$^{-1}$)")
    #hidexdecorations!(ax1, grid=false, ticks=false)

    ax2 = Axis(gd[1,2], xticks=[0, 0.005, 0.010])
    #hidexdecorations!(ax2, grid=false, ticks=false)
    hideydecorations!(ax2, grid=false, ticks=false)

    linkaxes!(ax1, ax2)

    ylims!(ax2, 0, 20)
    xlims!(ax2, -0.002, 0.015)

    wvps0 = [(60, 0.0, 5/6), (120, 0.0, 4/6), (240, 0.0, 3/6), (360, 0.0, 2/6), (480, 0.0, 1/6)]
    wvps = [(60, 0.0055, 5/6), (120, 0.0055, 4/6), (240, 0.0055, 3/6), (360, 0.0055, 2/6), (480, 0.0055, 1/6)]

    resonant = Makie.wong_colors()[3]
    blueshifted = Makie.wong_colors()[1]
    redshifted = Makie.wong_colors()[4]

    # Rabi analysis
    plot_kwvp!(ax1, 0.1, 2.0, wvps=wvps0, ymax=20)
    plot_effective_vg!(ax1, 0.1, 1.8, c=blueshifted)
    plot_effective_vg!(ax1, 0.1, 2.0, c=resonant)
    plot_effective_vg!(ax1, 0.1, 2.2, c=redshifted)

    plot_kwvp!(ax2, 0.1, 2.0, wvps=wvps, ymax=20)
    plot_effective_vg!(ax2, 0.1, 1.8, c=blueshifted)
    plot_effective_vg!(ax2, 0.1, 2.0, c=resonant)
    plot_effective_vg!(ax2, 0.1, 2.2, c=redshifted)

    Label(gd[2,1:2], L"$q$ (nm$^{-1}$)")

    text!(ax1, 0.02, 1, align=(:left, :top), space=:relative, text="(a)", font=:bold, fontsize=20)
    text!(ax2, 0.02, 1, align=(:left, :top), space=:relative, text="(b)", font=:bold, fontsize=20)

    group_line = [LineElement(color=:black, linestyle=:solid), LineElement(color=:black, linestyle=:dash)]
    group_colors = [PolyElement(color=c, strokecolor=:transparent) for c in [redshifted, resonant, blueshifted]]
    Legend(gd[1,3], [group_line, group_colors], [["LP", "UP"], 
    [L"-0.2",  L"0.0", L"+0.2"]], 
    ["", L"$\Omega_R$ (eV)"], orientation=:vertical, titleposition=:top)

    rowgap!(gd, 1, 5)
    colgap!(gd, 1, 7)
    colgap!(gd, 2, 7)

    fig
end