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

    ax2 = Axis(gd[2,1], ylabel=L"v_0\;(\text{ps}^{-1})")
    v0_vs_σx_fixed_reldis!(ax2; reldis=0.2)

    ax3 = Axis(gd[2,2])
    v0_vs_σx_fixed_reldis!(ax3; reldis=0.4)

    linkyaxes!(ax2, ax3) 
    fig
end

function d_vs_reldis(rel_dis=0.2; σx::Int=120, show_std=false)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel=L"d = \sqrt{\langle x^2 \rangle} / a")

    colors = Makie.wong_colors()
    ΩRvals = [0.05, 0.1, 0.2, 0.3]
    smvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]
    found = false
    for (ΩR,c) in zip(ΩRvals, colors)
        sm = round(ΩR * rel_dis / 100, digits=3)
        if sm in smvals
            println("Relative disorder of $rel_dis found for ΩR = $ΩR (σM = $sm)")
            d_propagation!(ax, ΩR=ΩR, σx=σx, σM=sm, show_std=show_std, color=c, label=L"%$ΩR")
            found = true
        end
    end

    if found
        xlims!(ax, 0, 0.5)

        ax.title = L"\sigma_M / \Omega_R = %$rel_dis"
        axislegend(ax, ax, L"\Omega_R\;\text{(eV)}")
    else
        @warn "No data found for relative disorder $rel_dis"
    end

    fig
end

function d_vs_rabi(;ΩRvals::Vector=[0.05, 0.1, 0.2, 0.3], σx::Int=120, σM::Float64=0.005, show_std=false)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel=L"d = \sqrt{\langle x^2 \rangle} / a")

    colors = Makie.wong_colors()
    for (ΩR,c) in zip(ΩRvals, colors)
        d_propagation!(ax, ΩR=ΩR, σx=σx, σM=σM, show_std=show_std, color=c, label=L"%$ΩR", fit=true)
    end

    xlims!(ax, 0, 0.5)

    axislegend(ax, ax, L"\Omega_R\;\text{(eV)}")

    fig
end

function d_vs_disorder(;ΩR::Float64=0.05, σx::Int=120, σMvals::Vector=[0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5], show_std=false)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel=L"d = \sqrt{\langle x^2 \rangle} / a")

    colors = Makie.wong_colors()
    for (σM,c) in zip(σMvals, colors)
        d_propagation!(ax, ΩR=ΩR, σx=σx, σM=σM, show_std=show_std, color=c, label=L"%$σM", fit=true)
    end

    xlims!(ax, 0, 0.5)

    axislegend(ax, ax, L"\sigma_M\;\text{(eV)}", position=:lt)

    fig
end

function d_vs_σx(;ΩR::Float64=0.05, σxvals::Vector=[60, 120, 240, 360, 480], σM::Float64=0.005, show_std=false)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel=L"d = \sqrt{\langle x^2 \rangle} / a")

    colors = Makie.wong_colors()
    for (σx,c) in zip(σxvals, colors)
        d_propagation!(ax, ΩR=ΩR, σx=σx, σM=σM, show_std=show_std, color=c, label=L"%$σx")
    end

    xlims!(ax, 0, 0.5)

    axislegend(ax, ax, L"\sigma_x\;\text{(nm)}", position=:lt)

    fig
end

function d_propagation(; ΩR::Float64=0.1, σx::Int=120, σM::Float64=0.005, show_std=false, color=Makie.wong_colors()[1], fit=false)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel=L"d = \sqrt{\langle x^2 \rangle} / a")
    d_propagation!(ax, ΩR=ΩR, σx=σx, σM=σM, show_std=show_std, color=color, fit=fit)
    fig
end

function d_propagation!(ax::Axis; ΩR::Float64=0.1, σx::Int=120, σM::Float64=0.005, show_std=false, color=Makie.wong_colors()[1], label="", fit=false)

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5")

    # Get time values
    r1 = 0:0.005:0.5
    r2 = 0.51:0.01:5
    tvals = vcat(r1, r2)

    d = h5read(path, "$(Int(σx))_avg_d")
    σd = h5read(path, "$(Int(σx))_std_d")

    #scatter!(ax, tvals, d, color=color, label=label)
    lines!(ax, tvals, d, color=color, label=label)

    if fit
        a,b = get_linear_fit(r1, d[1:length(r1)])
        println(b)
        lfit = [a+b*t for t in tvals]
        lines!(ax, tvals, lfit, color=:black)
    end
end

function get_linear_fit(x, y)
    X = zeros(length(x), 2)
    X[:,1] .= 1
    X[:,2] .= x
    return X \ y
end

function v0_vs_rabi!(ax::Axis; σx=120, σM=0.005, color=Makie.wong_colors()[1], marker=:circle)


    r1 = 0:0.005:0.5

    ΩRvals = [0.05, 0.1, 0.2, 0.3]
    v0vals = zeros(4)

    for i in 1:4
	    Rstr = replace(string(ΩRvals[i]), "." => "p") 
        sm = replace(string(σM), "." => "p") 
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5")

        d = h5read(path, "$(Int(σx))_avg_d")
        a,b = get_linear_fit(r1, d[1:length(r1)])

        v0vals[i] = b
    end

    scatter!(ax, ΩRvals, v0vals, color=color, marker=marker)
    lines!(ax, ΩRvals, v0vals, color=color)
end

function v0_vs_dis!(ax::Axis; σx=120, ΩR=0.05, color=Makie.wong_colors()[1], marker=:circle, label="")


    r1 = 0:0.005:0.5
	Rstr = replace(string(ΩR), "." => "p") 
    smvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]

    v0vals = zeros(10)

    for i in 1:10
        sm = replace(string(smvals[i]), "." => "p") 
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5")

        d = h5read(path, "$(Int(σx))_avg_d")
        a,b = get_linear_fit(r1, d[1:length(r1)])

        v0vals[i] = b
    end

    scatter!(ax, smvals, v0vals, color=color, marker=marker, label=label)
    lines!(ax, smvals, v0vals, color=color)
end

function v0_vs_reldis()
    fig = Figure()
    ax = Axis(fig[1,1], xlabel=L"\text{Relative disorder}\;(\sigma_M/\Omega_R)", ylabel=L"v_0\;(\text{ps}^{-1})")
    v0_vs_reldis!(ax)
    xlims!(ax, 0, 2.05)

    axislegend(ax, [[MarkerElement(marker=m, color=:black) for m = (:circle, :rect, :cross, :utriangle, :star5)],
    [PolyElement(color=c, strokecolor=:transparent) for c in Makie.wong_colors()[1:4]]],
    [string.([60,120,240,360,480]), string.([0.05, 0.1, 0.2, 0.3])],
    [L"\sigma_x\;(\text{nm})", L"\Omega_R\;(\text{eV})"], orientation=:vertical, titleposition=:left, nbanks=2)

    #vlines!(ax, [0.2, 0.4],  color=:gray, linestyle=:dash)
    fig
end

function v0_vs_reldis!(ax::Axis)

    r1 = 0:0.005:0.5

    ΩRvals = [0.05, 0.1, 0.2, 0.3]
    σMvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]
    σxvals = [60, 120, 240, 360, 480]
    v0vals = zeros(length(ΩRvals), length(σMvals), length(σxvals))

    for i in eachindex(ΩRvals)
        for j in eachindex(σMvals)
            for k in eachindex(σxvals)
	            Rstr = replace(string(ΩRvals[i]), "." => "p") 
                sm = replace(string(σMvals[j]), "." => "p") 
                path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5")
                d = h5read(path, "$(σxvals[k])_avg_d")
                a,b = get_linear_fit(r1, d[1:length(r1)])
                v0vals[i,j,k] = b
            end
        end
    end

    mkers = [:circle, :rect, :cross, :utriangle, :star5]
    clrs = Makie.wong_colors()[1:4]
    for k in eachindex(σxvals)
        for i in eachindex(ΩRvals)
            reldis = σMvals / ΩRvals[i]
            v0 = v0vals[i,:,k]
            scatter!(ax, reldis, v0, marker=mkers[k], color=clrs[i])
            #if k == 3
            #    lines!(ax, reldis, v0, marker=mkers[k], color=clrs[i])
            #end
        end
    end
end

function v0_vs_σx_fixed_reldis(;reldis=0.2)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel=L"\text{Initial wave packet spread}\;\sigma_x\;()", ylabel=L"v_0\;(\text{ps}^{-1})")
    v0_vs_σx_fixed_reldis!(ax; reldis=reldis)
    fig
end

function v0_vs_σx_fixed_reldis!(ax::Axis; reldis=0.2)

    r1 = 0:0.005:0.5
    mkers = [:circle, :rect, :cross, :utriangle, :star5]

    ΩRvals = [0.05, 0.1, 0.2, 0.3]
    σMvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]
    σxvals = [60, 120, 240, 360, 480]
    v0vals = zeros(5)

    for i in eachindex(ΩRvals)
        σM = round(reldis * ΩRvals[i], digits=3)
        println(σM)
        if !(σM in σMvals)
            continue
        end
	    Rstr = replace(string(ΩRvals[i]), "." => "p") 
        sm = replace(string(σM), "." => "p") 
        v0vals .= 0
        for k in eachindex(σxvals)
            path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5")
            d = h5read(path, "$(σxvals[k])_avg_d")
            a,b = get_linear_fit(r1, d[1:length(r1)])
            v0vals[k] = b
        end
        scatter!(ax, σxvals, v0vals)
        lines!(ax, σxvals, v0vals)
    end
end