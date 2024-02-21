function fig5()

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

    ylims!(ax2, 0, 12)
    xlims!(ax2, -0.002, 0.015)

    wvps0 = [(60, 0.0, 5/6), (120, 0.0, 4/6), (240, 0.0, 3/6), (360, 0.0, 2/6), (480, 0.0, 1/6)]
    wvps = [(60, 0.0055, 5/6), (120, 0.0055, 4/6), (240, 0.0055, 3/6), (360, 0.0055, 2/6), (480, 0.0055, 1/6)]

    c1 = Makie.wong_colors()[3]
    c2 = Makie.wong_colors()[2]
    #c3 = Makie.wong_colors()[1]
    #c4 = Makie.wong_colors()[4]

    # Rabi analysis
    plot_kwvp!(ax1, 0.1, 2.0, wvps=wvps0, ymax=12)
    plot_effective_vg!(ax1, 0.1, 2.0, c=c1)
    plot_effective_vg!(ax1, 0.2, 2.0, c=c2)

    plot_kwvp!(ax2, 0.1, 2.0, wvps=wvps, ymax=12)
    plot_effective_vg!(ax2, 0.1, 2.0, c=c1)
    plot_effective_vg!(ax2, 0.2, 2.0, c=c2)

    Label(gd[2,1:2], L"$q$ (nm$^{-1}$)")

    text!(ax1, 0.02, 1, align=(:left, :top), space=:relative, text="(a)", font=:bold, fontsize=20)
    text!(ax2, 0.02, 1, align=(:left, :top), space=:relative, text="(b)", font=:bold, fontsize=20)

    group_line = [LineElement(color=:black, linestyle=:solid), LineElement(color=:black, linestyle=:dash)]
    group_colors = [PolyElement(color=c, strokecolor=:transparent) for c in [c1, c2]]
    Legend(gd[1,3], [group_line, group_colors], [["LP", "UP"], 
    [L"0.1",  L"0.2"]], 
    ["", L"$\Omega_R$ (eV)"], orientation=:vertical, titleposition=:top)

    rowgap!(gd, 1, 5)
    colgap!(gd, 1, 7)
    colgap!(gd, 2, 7)

    fig
end

function plot_kwvp!(ax, ΩR, Em; wvps, ymax)
    path = joinpath(@__DIR__, "../../../dispersion_curves/effective_group_velocity/group_vel.h5")
    Rstr = replace(string(ΩR), "." => "p") 
    Emstr = replace(string(Em), "." => "p") 

    qvals = h5read(path, "qvals_R$(Rstr)_Em$(Emstr)")

    # Get k-space wave packet
    for w in wvps
        σx = w[1]
        q0 = w[2]
        h = w[3]

        qstr = replace(string(q0), "." => "p") 
        kwvp = h5read(path, "kwvp_q$(qstr)_sm$(σx)_R$(Rstr)_Em$(Emstr)")
        kwvp = kwvp ./ maximum(kwvp)
        if q0 > 0 
            cmap = [(:salmon3, t) for t in kwvp]
        else
            cmap = [(:steelblue2, t) for t in kwvp]
        end
        lines!(ax, qvals, repeat([ymax*h], length(qvals)), color=cmap, linewidth=12)
        text!(ax, 0.01, h+0.02, align=(:left, :bottom), space=:relative, text=L"$\sigma_x = %$σx$ nm", fontsize=15)
    end
end

function plot_effective_vg!(ax, ΩR, Em; c)

    path = joinpath(@__DIR__, "../../../dispersion_curves/effective_group_velocity/group_vel.h5")
    Rstr = replace(string(ΩR), "." => "p") 
    Emstr = replace(string(Em), "." => "p") 

    ∂UP = h5read(path, "effvg_UP_R$(Rstr)_Em$(Emstr)")
    ∂LP = h5read(path, "effvg_LP_R$(Rstr)_Em$(Emstr)")
    qvals = h5read(path, "qvals_R$(Rstr)_Em$(Emstr)")

    #scatter!(ax, LP.qvals, ∂LP, label="LP")
    lines!(ax, qvals, ∂LP, label="LP", linewidth=3, color=c)

    #scatter!(ax, UP.qvals, ∂UP, label="UP")
    lines!(ax, qvals, ∂UP, label="UP", linewidth=3, linestyle=:dash, color=c)
end