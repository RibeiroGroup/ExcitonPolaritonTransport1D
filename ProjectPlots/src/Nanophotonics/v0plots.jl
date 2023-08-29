using Makie
using HDF5
using LaTeXStrings

function v0_vs_rabi!(ax::Axis; σx=120, σM=0.005, color=Makie.wong_colors()[1], marker=:circle)


    r1 = 0:0.005:0.5

    ΩRvals = [0.05, 0.1, 0.2, 0.3]
    v0vals = zeros(4)

    for i in 1:4
	    Rstr = replace(string(ΩRvals[i]), "." => "p") 
        sm = replace(string(σM), "." => "p") 
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

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
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

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
                path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")
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
            path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")
            d = h5read(path, "$(σxvals[k])_avg_d")
            a,b = get_linear_fit(r1, d[1:length(r1)])
            v0vals[k] = b
        end
        scatter!(ax, σxvals, v0vals, marker=mkers)
        lines!(ax, σxvals, v0vals)
    end
end

function linearity_ratio(;ΩRvals::Vector=[0.05, 0.1, 0.2, 0.3], σx::Int=120, σM::Float64=0.005)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel=L"Fit residual ($R$)")

    colors = Makie.wong_colors()
    for (ΩR,c) in zip(ΩRvals, colors)
        linearity_ratio!(ax, ΩR=ΩR, σx=σx, σM=σM, color=c, label=L"%$ΩR")
    end

    xlims!(ax, 0, 2)
    ylims!(ax, 0.7, 1.05)
    axislegend(ax, ax, L"\Omega_R\;\text{(eV)}", position=:rt)
    hlines!(ax, [0.95], linestyle=:dot, color=:gray)
    fig
end

function linearity_ratio!(ax::Axis; ΩR::Float64=0.1, σx::Int=120, σM::Float64=0.005, color=Makie.wong_colors()[1], label="", polyorder=8)

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    # Get time values
    r1 = 0:0.005:0.5
    r2 = 0.51:0.01:5
    tvals = vcat(r1, r2)

    d = h5read(path, "$(Int(σx))_avg_d")
    σd = h5read(path, "$(Int(σx))_std_d")

    coef = polyfit(tvals, d, polyorder)
    dsmooth = [eval_poly(t, coef) for t in tvals]

    idx, tfinal = ProjectPlots.Nanophotonics.complete_span(ΩR, σM, σx, 0.01)

    t1 = 0.1:0.005:0.5
    tfits = vcat(t1, r2)
    Rvals = zeros(length(tfits))
    f = true
    for i in eachindex(tfits)
        tt = tfits[i]
        idx = 20 + i
        cs = polyfit(tvals[1:idx], dsmooth[1:idx], 1)
        #lines!(ax, tvals, [eval_poly(t, cs) for t in tvals])
        Rvals[i] = residual(tvals[1:idx], dsmooth[1:idx], cs)
    end

    m,mi = findmax(Rvals)
    println(mi)
    critical = findfirst(x->x<0.95, Rvals[mi:end]) + mi
    println(tfits[critical])
    vlines!(ax, [tfits[critical]], linestyle=:dash, color=color)

    #cs = complete_span(ΩR, σM, σx, 0.01)
    #vlines!(ax, [cs[2]], linestyle=:solid)

    lines!(ax, tfits, Rvals, label=label, color=color)
end