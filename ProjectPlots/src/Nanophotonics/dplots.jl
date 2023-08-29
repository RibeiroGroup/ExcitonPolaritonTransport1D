using Makie
using HDF5
using LaTeXStrings


function rmsd_propagation!(ax::Axis; ΩR::Float64=0.1, σx::Int=120, σM::Float64=0.005, show_std=false, color=Makie.wong_colors()[1], label="", fit=false)

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    # Get time values
    r1 = 0:0.005:0.5
    r2 = 0.51:0.01:5
    tvals = vcat(r1, r2)

    # Note that a = 10 here!
    # with a = 10 and converting to μm we get a 10/1000 = 1/100 = 0.01

    rmsd = 0.01*h5read(path, "$(Int(σx))_avg_d")
    std = 0.01*h5read(path, "$(Int(σx))_std_d")

    if show_std
        band!(ax, tvals, rmsd .- std, rmsd .+ std, color=(color, 0.5))
    end
    lines!(ax, tvals, rmsd, color=color, label=label)

    if fit
        a,b = get_linear_fit(r1, d[1:length(r1)])
        println(b)
        lfit = [a+b*t for t in tvals]
        lines!(ax, tvals, lfit, color=color, linewidth=2, linestyle=:dot)
    end
end


# Below are unused
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

    xlims!(ax, 0, 2)
    ylims!(ax, 0, 2000)

    axislegend(ax, ax, L"\Omega_R\;\text{(eV)}", position=:lt)

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
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

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
        lines!(ax, tvals, lfit, color=color, linewidth=2, linestyle=:dot)
    end
end

function d_polyfit!(ax::Axis; ΩR::Float64=0.1, σx::Int=120, σM::Float64=0.005, color=Makie.wong_colors()[1], label="", order=1)

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    # Get time values
    r1 = 0:0.005:0.5
    r2 = 0.51:0.01:5
    tvals = vcat(r1, r2)

    d = h5read(path, "$(Int(σx))_avg_d")
    σd = h5read(path, "$(Int(σx))_std_d")

    #scatter!(ax, tvals, d, color=color, label=label)
    lines!(ax, tvals, d, color=color, label=label)
    a,b = get_linear_fit(r1, d[1:length(r1)])
    lfit = [a+b*t for t in tvals]
    #lines!(ax, tvals, lfit, color=color)

    coef = polyfit(tvals, d, order)
    println("Polynomial fit for ΩR = $ΩR")
    println(coef)
    pfit = [eval_poly(t, coef) for t in tvals]
    lines!(ax, tvals, pfit, color=color, linewidth=2, label=label)
end

function d_polyfit(;ΩRvals::Vector=[0.05, 0.1, 0.2, 0.3], σx::Int=120, σM::Float64=0.005, order=3)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel=L"d = \sqrt{\langle x^2 \rangle} / a")

    colors = Makie.wong_colors()
    for (ΩR,c) in zip(ΩRvals, colors)
        d_polyfit!(ax, ΩR=ΩR, σx=σx, σM=σM, color=c, label=L"%$ΩR", order=order)
    end

    xlims!(ax, 0, 5)
    ylims!(ax, 0, 2000)

    axislegend(ax, ax, L"\Omega_R\;\text{(eV)}", position=:lt, merge=true)

    fig
end

function breakingpoint(coef, tvals)
    lin = [coef[1] + coef[2]*t for t in tvals]
    nonlin = [eval_poly(t, coef) for t in tvals] .- lin

    diff = abs.(lin .- nonlin)
    val, idx = findmin(diff)
    return tvals[idx]
end

function complete_span(ΩR, σM, σx, p=0.001)
	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")
    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5
    tvals = vcat(r1, r2)
    wvp = h5read(path, "$(Int(σx))_avg_wvp")

    for i in eachindex(tvals)
        if (wvp[1,i]+wvp[end,i]) ≥ p
            return (i, tvals[i])
        end
    end

    return (Inf, Inf)
end

function d_log!(ax::Axis; ΩR::Float64=0.1, σx::Int=120, σM::Float64=0.005, show_std=false, color=Makie.wong_colors()[1], label="", fit=false)

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    # Get time values
    r1 = 0:0.005:0.5
    r2 = 0.51:0.01:5
    tvals = vcat(r1, r2)

    d = h5read(path, "$(Int(σx))_avg_d")

    logd = [log(abs((10*d[i])^2 - (10*d[1])^2)) for i = 2:length(d)]
    τ = log.(tvals[2:end])

    der = [(logd[i] - logd[i-1])/(τ[i] - τ[i-1]) for i = 2:length(logd)]


    #scatter!(ax, tvals, d, color=color, label=label)
    #lines!(ax, τ[2:end], der, color=color, label=label)
    lines!(ax, τ, logd, color=color, label=label)
end