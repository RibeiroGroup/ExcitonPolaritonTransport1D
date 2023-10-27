function fig3()
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 400))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1], ylabel=L"v_0\;(\mu\text{m}\cdot \text{ps}^{-1})")
    ax2 = Axis(gd[1,2])
    hideydecorations!(ax2, ticks=false, grid=false)
    v0_vs_reldis!(ax1, σx=120)
    v0_vs_reldis!(ax2, σx=480)
    xlims!(ax1, 0, 2.05)
    ylims!(ax1, 0, 20)

    xlims!(ax2, 0, 2.05)
    ylims!(ax2, 0, 20)

    text!(ax1, 0.01, 0.95, text=L"(a) $\sigma_x = 120$ nm", space=:relative, align=(:left, :top), font=:bold)
    text!(ax2, 0.01, 0.95, text=L"(b) $\sigma_x = 480$ nm", space=:relative, align=(:left, :top), font=:bold)
    Label(gd[2,1:2], L"\text{Relative disorder}\;(\sigma_M/\Omega_R)")
    rowgap!(gd, 1, 3)
    axislegend(ax2, ax2, L"$\Omega_R$ (eV)", merge=true)

    fig
end