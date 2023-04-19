using Makie
using WGLMakie
using Statistics
using HDF5
using LaTeXStrings
using Measurements

"""
    Produces a graph with two columns and several plots of the mode weight at different detunings
"""
function plot_ideal_mw(;δvals=[0.1, 0.0, -0.1, -0.2], σx=120, Em=2.0, zeroq=true)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    δvals = reshape(δvals, length(δvals)÷2, 2)

    fig = Figure()
    axs = [Axis(fig[i,j], title=L"\delta = %$(δvals[i,j])\; \mathrm{eV}", xticks=2:0.1:2.4, xticklabelsvisible=(i==2),
    yticklabelsvisible=(j==1)) for i = axes(δvals, 1), j = axes(δvals, 2)]

    for (ax, δ) in zip(axs, δvals)
        plot_ideal_mw!(ax, σx=σx, Em=Em-δ, zeroq=zeroq)
        xlims!(ax, 2, 2.5)
    end

    Label(fig[1:end,0], "Relative Weight", rotation=π/2)
    Label(fig[end+1,1:end], "Mode energy (eV)")

    #fig[end+1, 1:end] = Legend(fig[end+1, 1:end], axs[1], L"\Omega_R\;\mathrm{(eV)}", orientation=:horizontal)
    axislegend(axs[1], axs[1], L"\Omega_R\;\mathrm{(eV)}")

    fig
end

function plot_ideal_mw!(ax::Axis; σx=120, Em=2.0, zeroq=true)

    e = h5read(joinpath(@__DIR__, "../../mode_weight/cavity.h5"), "e")

    even = [i%2 == 0 for i = eachindex(e)]
    odd = even .== false

    mstr = zeroq ? "zero" : "nonzero"

    i = 1
    for ΩR in [0.05, 0.1, 0.2, 0.3]
        Rstr = "R" * replace(string(ΩR), "."=>"p")
        Estr = "Em" * replace(string(Em), "."=>"p")
        mw = h5read(joinpath(@__DIR__,"../../mode_weight/$(mstr)_exc_momentum/$Estr/$Rstr/out.h5"), "$(σx)_phot_weight")
        lines!(ax, e[even], mw[even], label=L"%$(ΩR)", linewidth=2.5, color=Cycled(i))
        lines!(ax, e[odd], mw[odd], linestyle=:dash, linewidth=2.5, color=Cycled(i))
        i += 1
    end
end

function plot_ideal_mw_q!(ax::Axis; σx=120, Em=2.0, zeroq=true)

    e = h5read(joinpath(@__DIR__, "../../mode_weight/cavity.h5"), "e")
    w = h5read(joinpath(@__DIR__, "../../mode_weight/cavity.h5"), "w")

    perm = sortperm(w)

    mstr = zeroq ? "zero" : "nonzero"

    for ΩR in [0.05, 0.1, 0.2, 0.3]
        Rstr = "R" * replace(string(ΩR), "."=>"p")
        Estr = "Em" * replace(string(Em), "."=>"p")
        mw = h5read(joinpath(@__DIR__,"../../mode_weight/$(mstr)_exc_momentum/$Estr/$Rstr/out.h5"), "$(σx)_phot_weight")
        lines!(ax, w[perm], mw[perm], label=L"%$(ΩR)", linewidth=2.5)
    end
end

function barplot_dis_mw(;σx=120, σM=0.005, Em=2.0, ΩR=0.1, early=false, zeroq=true, Δω=0.02)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Mode Energy (eV)", ylabel="Average Relative Weight")
    xlims!(ax, 2, 2.55)

    plot_dis_mw!(ax, σx=σx, σM=σM, Em=Em, ΩR=ΩR, early=early, zeroq=zeroq, Δω=Δω)
    axislegend(ax)

    fig
end

function barplot_dis_mw!(ax::Axis; σx=120, σM=0.005, Em=2.0, ΩR=0.1, early=false, zeroq=true, Δω=0.02)

    # Translate numerical δ and ΩR to a string, for file handling 
    Emstr = "Em" * replace(string(Em), '.'=>'p')
    Rstr = "R" * replace(string(ΩR), '.'=>'p')

    cav_e = h5read(joinpath(@__DIR__, "../../mode_weight/cavity.h5"), "e")

    ER = early ? "EARLY_" : ""

    even = [i%2 == 0 for i = eachindex(cav_e)]
    odd = even .== false

    # Upbound hold the maximum (exclusive) mode energy grouped in each bar
    upbounds = (2 + Δω):Δω:3.5

    σMstr = "sm" * replace(string(σM), '.'=>'p')
    path = joinpath(@__DIR__, "../../mode_weight/disorder/$Emstr/$Rstr/$σMstr/out.h5")

    # Arrays to hold the sum of weights for each bar
    mw_even = zeros(Measurement, length(upbounds))
    mw_odd = zeros(Measurement, length(upbounds))           

    q = zeroq ? "" : "nzq_"
    avg_mw = h5read(path, "$(ER)$(q)$(σx)_avg_mode_weight")
    std_mw = h5read(path, "$(ER)$(q)$(σx)_std_mode_weight")

    # Use twice the standard deviation as error bars
    mw = [measurement(avg_mw[i], 2 .* std_mw[i]) for i = eachindex(avg_mw)]

    # Loop thorugh even mode weights
    k = 1
    for u in eachindex(upbounds)
        # If the energy of that mode is below the u-th upbound, add to the corresponding bar
        while cav_e[even][k] < upbounds[u]
            mw_even[u] += mw[even][k]
            k += 1
            # Since we are using an explicit iterator (k)
            # Add a break statement to prevent from going being the array length
            if k > length(cav_e[even])
                break
            end
        end
        # Once we go over the u-th ceiling, iterate to the next bar
    end

    # Repeat process for odd modes
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

    # Normalized using the max bar
    N = max(maximum(mw_even), maximum(mw_odd))

    mw_even = mw_even ./ N
    mw_odd = mw_odd ./ N

    vals_even = [x.val for x = mw_even]
    vals_odd  = [x.val for x = mw_odd]
    err_even = [x.err for x = mw_even]
    err_odd  = [x.err for x = mw_odd]

    barplot!(ax, upbounds, vals_even, label=L"q\;>\;0")
    barplot!(ax, upbounds, vals_odd, label=L"q \leq 0")
    #scatter!(ax, cav_e[even], avg_mw[even], label=L"q\;>\;0")
    #errorbars!(ax, cav_e[even], avg_mw[even], std_mw[even])
	errorbars!(ax, upbounds, vals_even, err_even, color = :red4, whiskerwidth = 10)
	errorbars!(ax, upbounds, vals_odd, err_odd, color = :red4, whiskerwidth = 10)
end

function plot_dis_mw(;σx=120, σM=0.005, Em=2.0, ΩR=0.1, early=false, zeroq=true)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Mode Energy (eV)", ylabel="Average Relative Weight")
    xlims!(ax, 2, 2.55)

    plot_dis_mw!(ax, σx=σx, σM=σM, Em=Em, ΩR=ΩR, early=early, zeroq=zeroq, Δω=Δω)
    axislegend(ax)

    fig
end

function plot_dis_mw!(ax::Axis;σx=120, σM=0.005, Em=2.0, ΩR=0.1, early=false, zeroq=true)

    # Translate numerical δ and ΩR to a string, for file handling 
    Emstr = "Em" * replace(string(Em), '.'=>'p')
    Rstr = "R" * replace(string(ΩR), '.'=>'p')

    cav_e = h5read(joinpath(@__DIR__, "../../mode_weight/cavity.h5"), "e")

    ER = early ? "EARLY_" : ""

    even = [i%2 == 0 for i = eachindex(cav_e)]
    odd = even .== false

    σMstr = "sm" * replace(string(σM), '.'=>'p')
    path = joinpath(@__DIR__, "../../mode_weight/disorder/$Emstr/$Rstr/$σMstr/out.h5")

    maxp = zeroq ? maximum(h5read(path, "q0_$(σx)_phot_cont")) : maximum(h5read(path, "nzq_$(σx)_phot_cont"))
    println("Maximum photon probability for Em=$Em σx=$σx ΩR=$ΩR σM=$σM =>  $maxp")

    q = zeroq ? "" : "nzq_"
    avg_mw = h5read(path, "$(ER)$(q)$(σx)_avg_mode_weight")
    std_mw = h5read(path, "$(ER)$(q)$(σx)_std_mode_weight")

    # Use twice the standard deviation as error bars
    mw = [measurement(avg_mw[i], std_mw[i]) for i = eachindex(avg_mw)]
    mw_even = mw[even]
    mw_odd = mw[odd]

    # Normalized using the maximum weight
    mw_even = mw_even ./ maximum(avg_mw)
    mw_odd = mw_odd ./ maximum(avg_mw)

    vals_even = [x.val for x = mw_even]
    vals_odd  = [x.val for x = mw_odd]
    err_even = [x.err for x = mw_even]
    err_odd  = [x.err for x = mw_odd]

    c1 = :royalblue3
    c2 = :darkgoldenrod

    band!(ax, cav_e[even], vals_even .- err_even, vals_even .+ err_even, transparency=true, color=(c1, 0.5))
    band!(ax, cav_e[odd], vals_odd .- err_odd, vals_odd .+ err_odd, transparency=true, color=(c2, 0.5))
    lines!(ax, cav_e[even], vals_even, linewidth=3, color=c1)
    lines!(ax, cav_e[odd], vals_odd, linewidth=3, color=c2)
end
