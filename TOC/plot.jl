using Makie
using WGLMakie
using HDF5

function get_wvp(t) 

    mr = 402:5401
    pr = 1:401

    fig = Figure(backgroundcolor = :transparent)
    ax1 = Axis(fig[1,1], backgroundcolor=:transparent)
    xlims!(1500, 3500)
    ylims!(0,0.009)

    wvp = h5read("out.h5", "wvp")

    exc = wvp[mr, t]
    println("EXC: $(sum(exc))")

    hidedecorations!(ax1)
    hidespines!(ax1)

    band!(ax1, 1:5000, 0, exc, color=:royalblue2)
    lines!(ax1, 1:5000, exc, color=:black, linewidth=2)

    fig
end

function get_phot(t) 

    mr = 402:5401
    pr = 1:401

    fig = Figure(backgroundcolor = :transparent)
    ax2 = Axis(fig[1,1], backgroundcolor=:transparent)
    ylims!(ax2, 0,0.009)
    xlims!(ax2, -0.012, 0.012)

    wvp = h5read("out.h5", "wvp")
    cav_w = h5read("cavity.h5", "w")

    s = sortperm(cav_w)

    phot = wvp[pr, t]

    hidedecorations!(ax2)
    hidespines!(ax2)

    band!(ax2, cav_w[s], 0, phot[s], color=:lightgoldenrod1)
    lines!(ax2, cav_w[s], phot[s], color=:black, linewidth=2)
    fig
end
