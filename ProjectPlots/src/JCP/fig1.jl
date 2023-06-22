using Makie
using HDF5

function fig1(; ΩR=0.1, σx=120, σMvals=[0.005, 0.02, 0.1], tvals=[0.1, 0.5, 1], h=1)
    fig = Figure(resolution=(800, 400))
    axs = [Axis(fig[1,1]),
           Axis(fig[1,2]),
           Axis(fig[1,3])]

    colors = [Makie.wong_colors()[i] for i = [1,3,4]]
    for i = 1:3
        σM = σMvals[i]
        c = colors[i]
        ax = axs[i]
        yshift = 0.0
        for t in tvals
            get_wvp!(ax, ΩR=ΩR, σx=σx, σM=σM, t=t, normalize=true, yshift=yshift, color=c)
            yshift += h
        end
    end

    w = 300
    for ax in axs
        xlims!(ax, 500-w,500+w)
        ylims!(ax, -0.1, 2.6*h)
        hideydecorations!(ax)
        hidexdecorations!(ax, label=false, ticks=false, ticklabels=false)
        ax.xticks = [300, 500, 700]
    end

    return fig
end

function get_wvp(; ΩR=0.1, σx=120, σM=0.005, t=0, yshift=0.0, normalize=false)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Distance (μm)", ylabel="Probability")

    get_wvp!(ax, ΩR=ΩR, σx=σx, σM=σM, t=t, yshift=yshift, normalize=normalize)

    return fig
end

function get_wvp!(ax::Axis; ΩR=0.1, σx=120, σM=0.005, t=0, normalize=false, yshift=0.0, color=Makie.wong_colors()[1])

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5")

    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5
    tvals = vcat(r1, r2)
    idx = findfirst(x->x==t, tvals)

    wvp = h5read(path, "$(Int(σx))_avg_wvp", (:, idx))
    P = sum(wvp)
    if normalize
        wvp = wvp ./ P
    end
    wvp .+= yshift
    println("Pmol = $P")

    barplot!(ax, 0:10:990, wvp, fillto=yshift, color=color)
end