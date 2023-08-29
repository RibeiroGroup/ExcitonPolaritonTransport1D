using Makie
using HDF5
using LaTeXStrings
using Statistics
using Measurements

function fig4()

    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(600, 700))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    axs = [Axis(gd[i,j]) for i = 1:2, j = 1:2]
    for (k,σx) in enumerate([60, 240, 120, 480])
        for (i,σM) in enumerate([0.005, 0.01, 0.02, 0.04, 0.1, 0.2])
            v0_vs_detuning!(axs[k], σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i))
        end

        ylims!(axs[k], 0, 13)

        if k == 1 || k == 3
            axs[k].xticklabelsvisible = false
        end

        if k > 2
            axs[k].yticklabelsvisible = false
        end
        lett = ['a', 'c', 'b', 'd'][k]
        text!(axs[k], 0.05, 0.99, text=L"(%$lett) $\sigma_x = %$(σx)$ nm", align=(:left, :top), space=:relative, fontsize=20, font=:bold)
    end

    Label(gd[3,1:2], L"$\delta = E_\text{min} - E_M$ (eV)")
    Label(gd[1:2,0], L"$v_0$ ($\mu$m ps$^{-1}$)", rotation=π/2, fontsize=25)
    Legend(gd[4,1:2], axs[1], L"\sigma_M/\Omega_R", orientation=:horizontal, merge=true, titleposition=:left)

    rowgap!(gd, 2, 10)
    rowgap!(gd, 3, 10)
    colgap!(gd, 1, 10)
    fig
end

function v0_vs_detuning!(ax; ΩR=0.1, σx=480, σM=0.005, Emvals=[1.8, 1.9, 2.0, 2.1, 2.2, 2.5], color=Makie.wong_colors()[1], label="")

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 

    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5

    v0 = zeros(length(Emvals))

    for i in eachindex(v0)
        δstr = replace(string(Emvals[i]), "." => "p")
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em$δstr/R$Rstr/sm$sm/out.h5")
        d = h5read(path, "$(Int(σx))_avg_d")
        std = h5read(path, "$(Int(σx))_std_d")
        a,b = get_linear_fit(r1, d[1:length(r1)])
        v0[i] = b / 100 # Convert from 1/ps to μm/ps using a = 10 nm
    end

    lines!(ax, 2.0 .- Emvals, v0, color=color, label=label)
    scatter!(ax, 2.0 .- Emvals, v0, color=color, label=label)
end

function fig5()

    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(600, 700))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    axs = [Axis(gd[i,j]) for i = 1:2, j = 1:2]
    for (k,σx) in enumerate([60, 240, 120, 480])
        for (i,σM) in enumerate([0.005, 0.01, 0.02, 0.04, 0.1, 0.3])
            maxd_vs_detuning!(axs[k], σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i))
        end

        ylims!(axs[k], 0, 25)

        if k == 1 || k == 3
            axs[k].xticklabelsvisible = false
        end

        if k > 2
            axs[k].yticklabelsvisible = false
        end

        lett = ['a', 'c', 'b', 'd'][k]
        text!(axs[k], 0.05, 0.99, text=L"(%$(lett)) $\sigma_x = %$(σx)$ nm", align=(:left, :top), space=:relative, fontsize=20, font=:bold)
    end

    Label(gd[3,1:2], L"Cavity Detuning $E_\text{min} - E_M$ (eV)")
    Label(gd[1:2,0], L"Maximum $d$ value ($\mu$m)", rotation=π/2, fontsize=25)
    Legend(gd[4,1:2], axs[1], L"\sigma_M/\Omega_R", orientation=:horizontal, merge=true, titleposition=:left)

    rowgap!(gd, 2, 10)
    rowgap!(gd, 3, 10)
    colgap!(gd, 1, 10)
    fig
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
        #maxd[i] = maximum(d) / 100 # Convert from unitless to μm using a = 10 nm = 1/100 μm
        maxd[i] = mean(d[451:end]) / 100 # Convert from unitless to μm using a = 10 nm = 1/100 μm
    end

    lines!(ax, 2.0 .- Emvals, maxd, color=color, label=label)
    scatter!(ax, 2.0 .- Emvals, maxd, color=color, label=label)
end