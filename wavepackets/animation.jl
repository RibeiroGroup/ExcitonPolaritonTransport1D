using WGLMakie
using Makie
using HDF5


function animate_wvp(;ΩR=0.1, σx=120, smvals=[0.02, 0.04, 0.1])

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

	Rstr = replace(string(ΩR), "." => "p") 
    smstr = [replace(string(s), "." => "p") for s in smvals]
    path1 = joinpath(@__DIR__, "../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$(smstr[1])/out.h5")
    path2 = joinpath(@__DIR__, "../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$(smstr[2])/out.h5")
    path3 = joinpath(@__DIR__, "../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$(smstr[3])/out.h5")

    r1 = 0:0.005:0.5

    time = Observable(0.0)

    y1 = @lift begin
        get_wvp(path1, σx, $time)
    end

    y2 = @lift begin
        wvp_idx = findfirst(x->x==$time, r1)
        h5read(path2, "$(Int(σx))_avg_wvp", (:, wvp_idx))
    end

    y3 = @lift begin
        wvp_idx = findfirst(x->x==$time, r1)
        h5read(path3, "$(Int(σx))_avg_wvp", (:, wvp_idx))
    end

    fig = Figure()
    ax1 = Axis(fig[1,1], title = @lift("t = $(round($time, digits = 2))"))
    ax2 = Axis(fig[2,1])
    ax3 = Axis(fig[3,1], xlabel="Distance (μm)", xticks=0:5:50)
    axs = [ax1, ax2, ax3]
    xlims!(ax1, 15,35)
    xlims!(ax2, 15,35)
    xlims!(ax3, 15,35)

    for (ax,sm) in zip(axs, smvals)
        reldis = 10*sm
        text!(ax, 0.01, 0.9, text=L"\sigma_M /\Omega_R = %$reldis", space=:relative, align=(:left, :top), font=:bold)
    end

    Label(fig[1:3, 0], L"|\Psi|^2", rotation=π/2)

    hidexdecorations!(ax1)
    hidexdecorations!(ax2)
    hideydecorations!.(axs)
    hidexdecorations!(ax3, label=false, ticklabels=false, ticks=false)

    barplot!(ax1, 0.25:0.5:49.755, y1)
    barplot!(ax2, 0.25:0.5:49.755, y2)
    barplot!(ax3, 0.25:0.5:49.755, y3)

    framerate = 20
    timestamps = range(0, 0.5, step=0.005)
    record(fig, "test.gif", timestamps;
            framerate = framerate) do t
        time[] = t
    end
end

function get_wvp(path,σx, t)
    r1 = 0:0.005:0.5
    idx = findfirst(x->x==t, r1)
    return h5read(path, "$(Int(σx))_avg_wvp", (:, idx))
end