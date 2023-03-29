using Makie
using WGLMakie
using HDF5

function plot_E2(ΩR, σM, σx)

    Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 

    zpm = h5read(joinpath(@__DIR__, "Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "$(σx)_avg_zpm")
    ppm = h5read(joinpath(@__DIR__, "Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5"), "$(σx)_avg_ppm")


    fig2 = Figure()
    ax2 = Axis(fig2[1,1], xlabel="Distance", ylabel=L"|E|^2")
    ylims!(ax2, ymin, ymax)
    xlims!(ax2, xmin, xmax)

    sl_t = Makie.Slider(fig2[2, 1], range = 1:110, horizontal = true, startvalue = 1)

    wv = lift(sl_t.value) do t
       	wvp[:,t]
    end

    σwv = lift(sl_t.value) do t
       	wvp_std[:,t]
    end

    tstr = lift(sl_t.value) do t
    	L"t = %$(wvp_t[t]) \; \mathrm{fs}"
    end

    barplot!(ax2, wv)
    errorbars!(ax2, 1:100, wv, σwv, color = :red4, whiskerwidth = 10)
    text!(ax2, tstr, space = :relative, position = Point2f(0.05, 0.65), fontsize=30)

    fig2
end