using Makie
using HDF5
using LaTeXStrings

function fig1(; ΩR=0.1, σx=120, σMvals=[0.005, 0.02, 0.1], tvals=[0.1, 0.5, 1])
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 400))

    # This code only works for a 3x3 grid (3 time and 3 disorder values)
    @assert length(tvals) == length(σMvals) == 3

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    axs = [Axis(gd[i,j]) for i = 1:3, j = 1:3]

    # Array of colors for each disorder value
    colors = [Makie.wong_colors()[i] for i = [1,3,4]]

    # w controls the width of the x axis relative to the center: xmin = 500-w, xmax = 500+w
    w = 260 

    # Loop through grid
    for i = 1:3
        for j = 1:3
            ax = axs[i,j]
            c = colors[j]
            σM = σMvals[j]
            t = tvals[i]

            # Create a plot in each axis
            get_wvp!(ax, ΩR=ΩR, σx=σx, σM=σM, t=t, normalize=true, color=c)

            # Set axis range
            xlims!(ax, 500-w,500+w)
            ylims!(ax, -0.12, 0.46)

            # Set x ticks
            ax.xticks = [300, 500, 700]

            # For the top row, create a label with relative disorder
            if i == 1
                rd = Int(round(100*σM/ΩR, digits=0))
                text!(ax, 0.15, 1, text="($('`'+j))", align=(:right, :top), space=:relative, font=:bold, fontsize=20)
                text!(ax, 0.18, 1, text="$rd%", align=(:left, :top), space=:relative, fontsize=20)
            end

            # Hide y-axes decorations, but keep labels on the first column
            if j == 1
                ax.ylabel = L"t = %$(Int(t*1000))\;\text{fs}"
                hideydecorations!(ax, label=false)
            else
                hideydecorations!(ax)
            end

            # Hide x-axes decorations, but keep ticks on the bottom row
            if i == 3
                hidexdecorations!(ax, ticks=false, ticklabels=false, label=false)
            else
                hidexdecorations!(ax)
            end

            # Hide spines, only keep a grid between each column
            spines = []
            if j != 1
                push!(spines, :l)
            end
            if i != 3
                push!(spines, :b)
            end
            if i != 1
                push!(spines, :t)
            end
            hidespines!(ax, spines...)
        end
    end

    # Label x-axis using the middle column
    axs[3,2].xlabel="Distance (nm)"

    # Remove gaps between columns and rows
    colgap!(gd, 0)
    rowgap!(gd, 0)

    # Done!
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
    d = h5read(path, "$(Int(σx))_avg_d")[idx]
    P = sum(wvp)
    if normalize
        wvp = wvp ./ P
    end
    wvp .+= yshift

    text!(ax, 500, -0.01, text=L"d \approx %$(Int(round(d, digits=0)))\;\; P_M = %$(round(P, digits=2))", align=(:center,:top), fontsize=18)

    barplot!(ax, 0:10:990, wvp, fillto=yshift, color=color)
end