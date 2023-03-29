using Makie
using WGLMakie
using Statistics
using HDF5
using LinearRegression
using LaTeXStrings
using Measurements

function plot_mw2(σx)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    Rvals = ["R0p05", "R0p1", "R0p2", "R0p3"]
    Emvals = ["Em1p9", "Em2p0", "Em2p1", "Em2p2"]
    NRvals = [0.05, 0.1, 0.2, 0.3]

    fig = Figure()
    #ax1p8 = Axis(fig[1,1], xlabel="Mode Energy (eV)", title=L"\delta = 0.2\; \mathrm{eV}", xticks=[2, 2.2, 2.4])
    ax1p9 = Axis(fig[1,1], xlabel="Mode Energy (eV)", ylabel="Relative mean prob.", title=L"\delta = 0.1\; \mathrm{eV}", xticks=[2, 2.2, 2.4],xminorticksvisible=true)
    ax2p0 = Axis(fig[1,2], xlabel="Mode Energy (eV)", ylabel="Relative mean prob.", title=L"\delta = 0.0\; \mathrm{eV}", xticks=[2, 2.2, 2.4],xminorticksvisible=true)
    ax2p1 = Axis(fig[2,1], xlabel="Mode Energy (eV)", ylabel="Relative mean prob.", title=L"\delta = -0.1\;\mathrm{eV}", xticks=[2, 2.2, 2.4],xminorticksvisible=true)
    ax2p2 = Axis(fig[2,2], xlabel="Mode Energy (eV)", ylabel="Relative mean prob.", title=L"\delta = -0.2\;\mathrm{eV}", xticks=[2, 2.2, 2.4],xminorticksvisible=true)
    axs = [ax1p9, ax2p0, ax2p1, ax2p2]

    (x-> xlims!(x, 2.0, 2.5)).(axs)
    e = h5read(joinpath(@__DIR__, "cavity.h5"), "e")

    i = 0
    for R in Rvals
        i += 1
        NR = NRvals[i]
        for k in eachindex(axs)
            E = Emvals[k]
            dist = h5read(joinpath(@__DIR__,"zero_exc_momentum/$E/$R/out.h5"), "$(σx)_phot_weight")
            lines!(axs[k], e, dist, label=L"\Omega_R = %$(NRvals[i]) \; \mathrm{eV}", linewidth=2.5)
        end
    end

    fig[3,1:2] = Legend(fig, ax2p0, orientation=:horizontal)

    fig
end

function plot_mw_NZ(σx)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    Rvals = ["R0p05", "R0p1", "R0p2", "R0p3"]
    Emvals = ["Em1p9", "Em2p0", "Em2p1", "Em2p2"]
    NRvals = [0.05, 0.1, 0.2, 0.3]

    fig = Figure()
    #ax1p8 = Axis(fig[1,1], xlabel="Mode Energy (eV)", title=L"\delta = 0.2\; \mathrm{eV}", xticks=[2, 2.2, 2.4])
    ax1p9 = Axis(fig[1,1], xlabel="Mode Energy (eV)", ylabel="Relative mean prob.", title=L"\delta = 0.1\; \mathrm{eV}", xticks=[2, 2.2, 2.4],xminorticksvisible=true)
    ax2p0 = Axis(fig[1,2], xlabel="Mode Energy (eV)", ylabel="Relative mean prob.", title=L"\delta = 0.0\; \mathrm{eV}", xticks=[2, 2.2, 2.4],xminorticksvisible=true)
    ax2p1 = Axis(fig[2,1], xlabel="Mode Energy (eV)", ylabel="Relative mean prob.", title=L"\delta = -0.1\;\mathrm{eV}", xticks=[2, 2.2, 2.4],xminorticksvisible=true)
    ax2p2 = Axis(fig[2,2], xlabel="Mode Energy (eV)", ylabel="Relative mean prob.", title=L"\delta = -0.2\;\mathrm{eV}", xticks=[2, 2.2, 2.4],xminorticksvisible=true)
    axs = [ax1p9, ax2p0, ax2p1, ax2p2]

    (x-> xlims!(x, 2.0, 3.0)).(axs)
    e = h5read(joinpath(@__DIR__, "cavity.h5"), "e")

    even = [i%2 == 0 for i = eachindex(e)]
    odd = even .== false
    clrs = Makie.wong_colors()[1:4]

    i = 0
    for R in Rvals
        i += 1
        NR = NRvals[i]
        for k in eachindex(axs)
            E = Emvals[k]
            dist = h5read(joinpath(@__DIR__,"nonzero_exc_momentum/$E/$R/out.h5"), "$(σx)_phot_weight")
            lines!(axs[k], e[even], dist[even], label=L"\Omega_R = %$(NRvals[i]) \; \mathrm{eV}", linewidth=2.5, color=clrs[i])
            lines!(axs[k], e[odd], dist[odd], linewidth=2.5, linestyle=:dash, color=clrs[i])
        end
    end

    fig[3,1:2] = Legend(fig, ax2p0, orientation=:horizontal)

    fig
end


function plot_dis_mw(;σx, δ, ΩR, early=false)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    σxvals = [0.005, 0.02, 0.05, 0.1]
    # Create figure and axes
    titles = ["(a)" "(e)"; "(b)" "(f)"; "(c)" "(g)"; "(d)" "(h)"]
    fig = Figure()
    axs = [Axis(fig[i,j], xticklabelsvisible= i == 4 ? true : false, 
    yticklabelsvisible= j == 1 ? true : false,  xticks = 2:0.1:2.5) for i = 1:4, j = 1:2]

    for ax in axs
        xlims!(ax, 2, 2.55)
    end

    for i = 1:4
        linkyaxes!(axs[i,1], axs[i,2])
    end

    Label(fig[0,1], L"\bar{q}_0 = 0", justification=:center, tellwidth=false, lineheight=0.8)
    Label(fig[0,2], L"\bar{q}_0 \neq 0", justification=:center, tellwidth=false, lineheight=0.8)
    Label(fig[1:end,0], "Relative average weight", justification=:center, rotation=pi/2)
    Label(fig[5,1:2], "Mode energy (eV)", justification=:center)

    # Translate numerical δ and ΩR to a string, for file handling 
    Emstr = "Em" * replace(string(2 - δ), '.'=>'p')
    Rstr = "R" * replace(string(ΩR), '.'=>'p')

    cav_e = h5read(joinpath(@__DIR__, "cavity.h5"), "e")

    ER = early ? "EARLY_" : ""

    even = [i%2 == 0 for i = eachindex(cav_e)]
    odd = even .== false

    upbounds = 2.02:0.02:3.5

    mw_even = zeros(Measurement, length(upbounds))
    mw_odd = zeros(Measurement, length(upbounds))

    for (i,σM) in enumerate([0.005, 0.02,  0.05, 0.1])
        σMstr = "sm" * replace(string(σM), '.'=>'p')

        for (j, q) in enumerate(["", "nzq_"])

            mw_even .= 0
            mw_odd .= 0

            avg_mw = h5read(joinpath(@__DIR__, "disorder/$Emstr/$Rstr/$σMstr/out.h5"), "$(ER)$(q)$(σx)_avg_mode_weight")
            std_mw = h5read(joinpath(@__DIR__, "disorder/$Emstr/$Rstr/$σMstr/out.h5"), "$(ER)$(q)$(σx)_std_mode_weight")

            mw = [measurement(avg_mw[i], 2 .* std_mw[i]) for i = eachindex(avg_mw)]

            k = 1
            for u in eachindex(upbounds)
                while cav_e[even][k] < upbounds[u]
                    mw_even[u] += mw[even][k]
                    k += 1
                    if k > length(cav_e[even])
                        break
                    end
                end
            end
            k = 1
            for u in eachindex(upbounds)
                while cav_e[odd][k] < upbounds[u]
                    mw_odd[u] += mw[odd][k]
                    k += 1
                    if k > length(cav_e[even])
                        break
                    end
                end
            end

            mw_even = mw_even ./ maximum(mw_even)
            mw_odd = mw_odd ./ maximum(mw_odd)

            vals_even = [x.val for x = mw_even]
            vals_odd  = [x.val for x = mw_odd]
            err_even = [x.err for x = mw_even]
            err_odd  = [x.err for x = mw_odd]

            barplot!(axs[i,j], upbounds, vals_even)
            barplot!(axs[i,j], upbounds, vals_odd)
	        errorbars!(axs[i,j], upbounds, vals_even, err_even, color = :red4, whiskerwidth = 10)

            vlines!(axs[i,j], [2 - δ], color = :red, linestyle=:dash)
            if q == "nzq_"
                vlines!(axs[i,j], [2.1], color = :blue, linestyle=:dash)
            end

            text!(axs[i,j], 0.95, 0.95, text=titles[i,j], align=(:right, :top), space=:relative, font=:bold)
            text!(axs[i,j], 0.99, 0.65, text=L"{\sigma_M}/{\Omega_R} = %$(10*σM)", align=(:right, :top), space=:relative, fontsize=18)
        end
    end

    fig
end

function plot_mw(;δ)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    # Create figure and axes
    titles = ["(a)" "(e)"; "(b)" "(f)"; "(c)" "(g)"; "(d)" "(h)"]
    fig = Figure()
    axs = [Axis(fig[i,j], xticklabelsvisible= i == 3 ? true : false, 
    yticklabelsvisible= j == 1 ? true : false,  xticks = 2:0.1:2.5) for i = 1:3, j = 1:2]

    for ax in axs
        xlims!(ax, 2, 2.55)
    end

    Label(fig[0,1], L"\bar{q}_0 = 0", justification=:center, tellwidth=false, lineheight=0.8)
    Label(fig[0,2], L"\bar{q}_0 \neq 0", justification=:center, tellwidth=false, lineheight=0.8)
    Label(fig[1:end,0], "Relative average weight", justification=:center, rotation=pi/2)
    Label(fig[4,1:2], "Mode energy (eV)", justification=:center)

    # Translate numerical δ and ΩR to a string, for file handling 
    Emstr = "Em" * replace(string(2.0 - δ), '.'=>'p')

    cav_e = h5read(joinpath(@__DIR__, "cavity.h5"), "e")

    even = [i%2 == 0 for i = eachindex(cav_e)]
    odd = even .== false

    for (i,σx) in enumerate([60, 120, 180])

        for (j, q) in enumerate(["", "non"])

            for (k,ΩR) in enumerate([0.05, 0.1, 0.2, 0.3])

                Rstr = "R" * replace(string(ΩR), '.'=>'p')

                mw = h5read(joinpath(@__DIR__,"$(q)zero_exc_momentum/$Emstr/$Rstr/out.h5"), "$(σx)_phot_weight")
                #lines!(axs[i,j], cav_e[even], mw[even], label=L"\Omega_R = %$(ΩR) \; \mathrm{eV}", linewidth=2.5)

                if q == ""
                    lines!(axs[i,j], cav_e, mw, label=L"\Omega_R = %$(ΩR) \; \mathrm{eV}", linewidth=2.5)
                elseif q == "non"
                    lines!(axs[i,j], cav_e[even], mw[even], label=L"\Omega_R = %$(ΩR) \; \mathrm{eV}", linewidth=2.5)
                    lines!(axs[i,j], cav_e[odd], mw[odd], linewidth=2.5, linestyle=:dash, color=Cycled(k))
                end

                vlines!(axs[i,j], [2 - δ], color = :red, linestyle=:dot)
                if q == "non"
                    vlines!(axs[i,j], [2.1], color = :blue, linestyle=:dot)
                end
            end

            text!(axs[i,j], 0.95, 0.95, text=titles[i,j], align=(:right, :top), space=:relative, font=:bold)
            text!(axs[i,j], 0.99, 0.7, text=L"\sigma_x / a = %$(Int(σx/10))", align=(:right, :top), space=:relative, fontsize=18)
        end
    end

    fig[5,1:2] = Legend(fig, axs[1,1], orientation=:horizontal)
    fig
end
