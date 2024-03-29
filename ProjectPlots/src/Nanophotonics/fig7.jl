using Makie
using HDF5
using LaTeXStrings
using Statistics
using Measurements

function fig7(;σx=240)
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 400))
    gd = fig[1,1] = GridLayout()

    ax1 = Axis(gd[1,1], ylabel=L"$v_0$ ($\mu$m$\cdot$ps$^{-1}$)", xticks=[-0.6, -0.4,-0.2, 0.0, 0.2])
    ax2 = Axis(gd[1,2], ylabel=L"Maximum RMSD ($\mu$m)", xticks=[-0.6, -0.4,-0.2, 0.0, 0.2])
    linkxaxes!(ax1, ax2)

    for (i,σM) in enumerate([0.005, 0.01, 0.02, 0.04, 0.1, 0.2])
    # Plot detuning effect on v0
        v0_vs_detuning!(ax1, σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i))
        #nzq_v0_vs_detuning!(ax1, σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i))

    # Plot detuning effect on max RMSD
        maxd_vs_detuning!(ax2, σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i))
    end

    Label(gd[2,1:2], L"Detuning $\delta = \hbar\omega_0 - E_M$ (eV)")
    Legend(gd[3,1:2], ax1, L"\sigma_M/\Omega_R", orientation=:horizontal, merge=true, titleposition=:left)

    text!(ax1, 0.01, 0.98, text="(a)", space=:relative, align=(:left, :top), font=:bold)#, fontsize=22)
    text!(ax2, 0.01, 0.98, text="(b)", space=:relative, align=(:left, :top), font=:bold)#, fontsize=22)

    rowgap!(gd, 1, Relative(0.008))
    rowgap!(gd, 2, Relative(0.008))

    fig
end

function v0_vs_detuning!(ax; ΩR=0.1, σx=480, σM=0.005, Emvals=[1.8, 1.9, 2.0, 2.1, 2.2, 2.5], color=Makie.wong_colors()[1], label="")

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 

    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5

    v0 = zeros(length(Emvals))

    println("----σM = $σM")
    for i in eachindex(v0)
        δstr = replace(string(Emvals[i]), "." => "p")
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em$δstr/R$Rstr/sm$sm/out.h5")
        d = h5read(path, "$(Int(σx))_avg_d")
        std = h5read(path, "$(Int(σx))_std_d")
        println("Em = $(Emvals[i])")
        a,b = get_linear_fit(r1, d[1:length(r1)])
        v0[i] = b / 100 # Convert from 1/ps to μm/ps using a = 10 nm
    end

    lines!(ax, 2.0 .- Emvals, v0, color=color, label=label)
    scatter!(ax, 2.0 .- Emvals, v0, color=color, label=label)
end

function maxd_vs_detuning!(ax; ΩR=0.1, σx=480, σM=0.005, Emvals=[1.8, 1.9, 2.0, 2.1, 2.2, 2.5], color=Makie.wong_colors()[1], label="")

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 

    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5

    maxd = zeros(length(Emvals))

    for i in eachindex(maxd)
        δstr = replace(string(Emvals[i]), "." => "p")
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em$δstr/R$Rstr/sm$sm/out.h5")
        d = h5read(path, "$(Int(σx))_avg_d")
        maxd[i] = maximum(d) / 100 # Convert from unitless to μm using a = 10 nm = 1/100 μm
        #maxd[i] = mean(d[451:end]) / 100 # Convert from unitless to μm using a = 10 nm = 1/100 μm
    end

    lines!(ax, 2.0 .- Emvals, maxd, color=color, label=label)
    scatter!(ax, 2.0 .- Emvals, maxd, color=color, label=label)
end

function nzq_v0_vs_detuning!(ax; σx=480, σM=0.005, Emvals=[1.8, 2.0, 2.2], color=Makie.wong_colors()[1], label="")

    sm = replace(string(σM), "." => "p") 

    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5

    v0 = zeros(length(Emvals))

    println("----σM = $σM")
    for i in eachindex(v0)
        δstr = replace(string(Emvals[i]), "." => "p")
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em$δstr/R0p1/non_zero_q/sm$sm/out.h5")
        d = h5read(path, "$(Int(σx))_avg_d")
        std = h5read(path, "$(Int(σx))_std_d")
        println("Em = $(Emvals[i])")
        a,b = get_linear_fit(r1, d[1:length(r1)])
        v0[i] = b / 100 # Convert from 1/ps to μm/ps using a = 10 nm
    end

    lines!(ax, 2.0 .- Emvals, v0, color=color, label=label)
    scatter!(ax, 2.0 .- Emvals, v0, color=color, label=label)
end

function propagation_and_v0(;ΩR, σx, σM, Em)

    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()

    ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel="RMSD")

    rmsd_propagation!(ax, ΩR=ΩR, σx=σx, σM=σM, Em=Em, fit=true)

    xlims!(ax, 0, 0.7)

    fig
end