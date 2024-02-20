using HDF5
using LaTeXStrings

function dispersion_plot()
    path = joinpath(@__DIR__, "out.h5")

    LP_k = h5read("dispersion_curves/out.h5", "2.0_LP_k")
    LP = h5read("dispersion_curves/out.h5", "2.0_LP")

    UP_k = h5read("dispersion_curves/out.h5", "2.0_UP_k")
    UP = h5read("dispersion_curves/out.h5", "2.0_UP")

    qvals = h5read("dispersion_curves/out.h5", "qvals")
    ephot = h5read("dispersion_curves/out.h5", "phot2")

    fontsize_theme = Theme(fontsize = 30)
    set_theme!(fontsize_theme)
    fig = Figure(backgroundcolor = :transparent)

    ax = Axis(fig[1,1], xlabel=L"q \; (\text{nm}^{-1})", ylabel=L"\text{Energy (eV)}", backgroundcolor = :transparent)

    hlines!(ax, [2.0], linestyle=:dash, color=:gray)
    lines!(ax, qvals, ephot, linestyle=:dash, color=:gray)

    ylims!(ax, 1.85, 2.3)
    xlims!(ax, -0.01, +0.01)

    lines!(ax, LP_k, LP, linewidth=5)
    lines!(ax, UP_k, UP, linewidth=5)

    text!(ax, 0.0, 1.9, align=(:center, :bottom), text="LP", color=Makie.wong_colors()[1])
    text!(ax, 0.0, 2.1, align=(:center, :bottom), text="UP", color=Makie.wong_colors()[2])
    text!(ax, 0.0075, 2.0, align=(:center, :bottom), text=L"E_M", color=:gray)
    text!(ax, 0.008, 2.15, align=(:center, :bottom), text=L"\hbar \omega_q", color=:gray)

    fig
end