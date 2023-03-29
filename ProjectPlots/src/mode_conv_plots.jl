using HDF5
using Makie
using WGLMakie
using Measurements
using Statistics
using PrettyTables

function get_path(Nm, Nc, ΩR, fname="out.h5")
    fpath  = "ideal/Nm$(Nm)_Nc$(Nc)_R$(replace(string(ΩR), '.'=>'p'))_a10_0_Em2_0_sx60/$fname"
    return joinpath(@__DIR__, fpath)
end

"""
    Produces three panels (corresponding to different system sizes) with d(t) plots from 0 - 500 fs for several number of cavity modes.
"""
function plot_propagation()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    # Create figure object
    fig = Figure()

    # List with number of molecules for each panel
    Nmvals = [5000, 10000, 20000]

    # List with number of cavity modes
    Ncvals = [50, 100, 200, 400, 500, 800]

    # Create axis into the figure
    axs = [Axis(fig[i,1]) for i = eachindex(Nmvals)]

    # Axis configuration
    for (ax, Nm) in zip(axs, Nmvals)
        ax.xlabel = "Time (fs)"
        ax.title  = L"N_M = %$Nm"
        ax.xticks = 0:100:500
        ax.yticks = [100, 300, 500]
        xlims!(ax, 0, 500)
        ylims!(ax, 0, 500)
    end

    # Common label to all y-axes
    Label(fig[1:end,0], L"d = \sqrt{\left\langle x^2 \right\rangle} / a", rotation = pi/2)

    # Get time values
    (ti,δt,tf) = h5read(get_path(5000, 800, 0.1), "time_range")

    # Convert pico -> femtoseconds
    t = collect(ti:δt:tf) .* 1000

    for (i,Nm) in enumerate(Nmvals)
        for Nc in Ncvals

            # Fetch mean squared disp
            x2 = h5read(get_path(Nm, Nc, 0.1), "mean_square_disp")

            # Compute d, note that a = 10
            d = (sqrt.(x2) .- sqrt(x2[1])) ./ 10

            # Produce lineplot on the correspoding axis
            lines!(axs[i], t, d, label = L"N_c = %$(2*Nc+1)")
        end
    end

    fig[1:3,2] = Legend(fig, axs[1])
    fig
end

"""
    Plots the error in d (with respect to a computation with 1601 cavity modes) as a function of the number of cavity modes 
    for several system sizes
"""
function plot_error()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    # Create figure object
    fig = Figure()

    # Create axis into the figure
    ax1 = Axis(fig[1,1], xlabel="Number of cavity modes", title="(a)", titlealign = :left)
    ax2 = Axis(fig[2,1], xlabel="Cavity energy cutoff (eV)", title="(b)", titlealign = :left)
    Label(fig[1:end,0], "Error", rotation=pi/2)
    
    # Axis configuration
    xlims!(ax1, 50, 1050)
    xlims!(ax2, (2,3.5))
    ylims!(ax2, (-0.05,0.9))
    linkyaxes!(ax1, ax2)

    # List with number of molecules for each panel
    Nmvals = [5000, 8000, 10000, 15000, 20000]

    # List with number of cavity modes
    Ncvals = [50, 100, 200, 400, 500]

    # Different markers to be used for each Nc value
    mkers = [:circle, :rect, :utriangle, :diamond, :star5]

    # Slice of the simulation used to compute error
    # 501 = 5 ps
    slice = 1:101
    println("AAAAAAAAAA")
    for (Nm, mk) in zip(Nmvals, mkers)

        # Get reference data - Nc = 1601 - Note that we are computing d
        ref = 0.1 .* sqrt.(h5read(get_path(Nm, 800, 0.1), "mean_square_disp")[slice]) 

        # Preallocate arrays for cavity energies and error values
        err_vals = zeros(5)
        cav_e = zeros(5)
        for (j,Nc) in enumerate(Ncvals)

            # Fetch data for given Nc
            d = 0.1 .* sqrt.(h5read(get_path(Nm, Nc, 0.1), "mean_square_disp")[slice])
            
            # Read output files to find the maximum cavity energy (i.e. energy cutoff)
            open(get_path(Nm, Nc, 0.1, "log.out"), "r") do io
                m = match(r"Maximum Cavity energy:\s+?(\d+?\.\d+)", read(io, String))
                if isnothing(m)
                    println("Couldn't match $(get_path(Nm, Nc, 0.1, "log.out"))")
                else
                    cav_e[j] = parse(Float64, m.captures[1])
                end
            end

            # Compute error w.r.t to reference as the mean of the relative absolute error
            err_vals[j] = mean(abs.(d .- ref) ./ ref)
        end

        # Produce scatter and line plots
        scatter!(ax1, (2 .* Ncvals .+ 1), err_vals, label="$(Nm)", marker=mk)
        lines!(ax1, (2 .* Ncvals .+ 1), err_vals)

        scatter!(ax2, cav_e, err_vals, marker=mk)

        # Print out human readable table
        println("Nm = $Nm")
        pretty_table(hcat(Ncvals, cav_e, [round(100*e, digits=2) for e in err_vals]);
        header = ["Nc", "Energy (eV)", "Error"])
    end

    lines!(ax2, 2:0.05:3.5, [0.8*exp(-9.5*(e-2)) for e in 2:0.05:3.5], color=:black, linestyle=:dash, label=L"\exp[-\alpha(E_\mathrm{cutoff}-E_\mathrm{min})]")

    #axislegend(title=L"N_\mathrm{M}", ax = ax1)
    Legend(fig[1:end,2], ax1, L"N_\mathrm{M}") 
    #tellheight = false,
    #tellwidth = false,)
    axislegend(ax2)

    fig
end
 
function plot_dis_propagation(σM)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    # Create figure object
    fig = Figure()

    # Create axis into the figure
    ax = Axis(fig[1,1], xlabel="Time (fs)", ylabel=L"d")#, title = L"N_M = 5000", xticks = 0:100:500, yticks = [100, 300, 500])

    tvals = 0:10:5000 

    for Nc in [50, 75, 100, 200, 400, 500, 800]

        d = h5read(joinpath(@__DIR__, "disorder/R0p1/$(σM)/Nc$Nc/out.h5"), "NR_100_sm60_avg_mode_weight")
        lines!(ax, tvals, d, label = L"N_c= %$(Nc)")
    end

    fig[1,2] = Legend(fig, ax)
    fig
end

function plot_realizations(σM)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (fs)", ylabel=L"d")#, title = L"N_M = 5000", xticks = 0:100:500, yticks = [100, 300, 500])

    tvals = 0:10:5000 

    for NR in [20, 40, 60, 80, 100]

        d = h5read(joinpath(@__DIR__, "disorder/R0p1/$(σM)/Nc800/out.h5"), "NR_$(NR)_sm60_avg_mode_weight")
        σ = h5read(joinpath(@__DIR__, "disorder/R0p1/$(σM)/Nc800/out.h5"), "NR_$(NR)_sm60_std_mode_weight")
        lines!(ax, tvals, d, label = L"N_R= %$(NR)")
        band!(ax, tvals, d .- σ, d .+ σ)
    end

    fig[1,2] = Legend(fig, ax)
    fig
end

function plot_dis_error(;NR = 100, σx = 60)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    # Create figure object
    #fig = Figure(resolution=(800,400))
    fig = Figure()

    # Create axis into the figure
    ax1 = Axis(fig[1,1], title="(a)", titlealign = :left)
    ax2 = Axis(fig[2,1], title="(b)", titlealign = :left, xlabel="Cavity energy cutoff (eV)")
    linkyaxes!(ax1, ax2)
    linkxaxes!(ax1, ax2)
    Label(fig[1:end,0], "Error", rotation=pi/2)

    # List with number of cavity modes
    Ncvals = [0, 1, 5, 10, 20, 35, 50, 75, 100]
    
    # Preallocate arrays for error values, error deviations and cavity energies
    err_vals = zeros(length(Ncvals))
    err_devs = zeros(length(Ncvals))
    cav_e = zeros(length(Ncvals))

    # Slice on time steps - this controls whether we average early or late propagation times
    slice = 1:101

    # Different markers to be used for each Nc value
    mkers = [:circle, :rect, :utriangle, :diamond]

    # Plot results for ΩR = 0.1 onto axis 1
    for (mk, σM) in zip(mkers, [0.005, 0.02, 0.05, 0.1])

        σMstr = "sm" * replace(string(σM), '.'=>'p')
        r = 10 * σM

        # Get data for Nc = 1601V
        ref = h5read(joinpath(@__DIR__, "disorder/Em2p0/R0p1/$(σMstr)/Nc800/out.h5"), "NR_$(NR)_sm$(σx)_avg_mode_weight")[slice]
        std_ref = h5read(joinpath(@__DIR__, "disorder/Em2p0/R0p1/$(σMstr)/Nc800/out.h5"), "NR_$(NR)_sm$(σx)_std_mode_weight")[slice]

        # Create Measurement objects, with the standard deviation representing the error of d
        ref_mmt = [measurement(ref[i], std_ref[i]) for i = eachindex(ref)]

        for (j,Nc) in enumerate(Ncvals)

            # Fetch d values along with their standard deviation
            d = h5read(joinpath(@__DIR__, "disorder/Em2p0/R0p1/$(σMstr)/Nc$Nc/out.h5"), "NR_$(NR)_sm$(σx)_avg_mode_weight")[slice]
            std = h5read(joinpath(@__DIR__, "disorder/Em2p0/R0p1/$(σMstr)/Nc$Nc/out.h5"), "NR_$(NR)_sm$(σx)_std_mode_weight")[slice]

            # Create Measurement objects, with the standard deviation representing the error of d
            d_mmt = [measurement(d[i], 2 .* std[i]) for i = eachindex(d)]

            # Read output files to find the maximum cavity energy (i.e. energy cutoff)
            open(joinpath(@__DIR__, "disorder/Em2p0/R0p1/$(σMstr)/Nc$Nc/log.out"), "r") do io
                m = match(r"Maximum Cavity energy:\s+?(\d+?\.\d+)", read(io, String))
                if isnothing(m)
                    println("Couldnt match $(get_path(Nm, Nc, 0.1, "log.out"))")
                else
                    cav_e[j] = parse(Float64, m.captures[1])
                end
            end

            # Compute error and associated uncertainty 
            err = mean(abs.(d_mmt .- ref_mmt) ./ ref_mmt)
            err_vals[j] = err.val
            err_devs[j] = err.err
        end

        # Produce plots
        scatter!(ax1, cav_e, err_vals, label="$r", marker=mk)
        lines!(ax1, cav_e, err_vals)
        errorbars!(ax1, cav_e, err_vals, err_devs, color=:red)

    end

    # Plot results for ΩR = 0.2 onto axis 2
    for (mk, σM) in zip(mkers, [0.01, 0.04, 0.1, 0.2])

        σMstr = "sm" * replace(string(σM), '.'=>'p')
        r =  5 * σM

        # Get data for Nc = 1601V
        ref = h5read(joinpath(@__DIR__, "disorder/Em2p0/R0p2/$(σMstr)/Nc800/out.h5"), "NR_$(NR)_sm$(σx)_avg_d")[slice]
        std_ref = h5read(joinpath(@__DIR__, "disorder/Em2p0/R0p2/$(σMstr)/Nc800/out.h5"), "NR_$(NR)_sm$(σx)_std_d")[slice]

        # Create Measurement objects, with the standard deviation representing the error of d
        ref_mmt = [measurement(ref[i], std_ref[i]) for i = eachindex(ref)]

        for (j,Nc) in enumerate(Ncvals)

            # Fetch d values along with their standard deviation
            d = h5read(joinpath(@__DIR__, "disorder/Em2p0/R0p2/$(σMstr)/Nc$Nc/out.h5"), "NR_$(NR)_sm$(σx)_avg_d")[slice]
            std = h5read(joinpath(@__DIR__, "disorder/Em2p0/R0p2/$(σMstr)/Nc$Nc/out.h5"), "NR_$(NR)_sm$(σx)_std_d")[slice]

            # Create Measurement objects, with the standard deviation representing the error of d
            d_mmt = [measurement(d[i], 2 .* std[i]) for i = eachindex(d)]

            # Read output files to find the maximum cavity energy (i.e. energy cutoff)
            open(joinpath(@__DIR__, "disorder/Em2p0/R0p2/$(σMstr)/Nc$Nc/log.out"), "r") do io
                m = match(r"Maximum Cavity energy:\s+?(\d+?\.\d+)", read(io, String))
                if isnothing(m)
                    println("Couldnt match $(get_path(Nm, Nc, 0.1, "log.out"))")
                else
                    cav_e[j] = parse(Float64, m.captures[1])
                end
            end

            # Compute error and associated uncertainty 
            err = mean(abs.(d_mmt .- ref_mmt) ./ ref_mmt)
            err_vals[j] = err.val
            err_devs[j] = err.err
        end

        # Produce plots
        scatter!(ax2, cav_e, err_vals, label=r, marker=mk)
        lines!(ax2, cav_e, err_vals)
        errorbars!(ax2, cav_e, err_vals, err_devs, color=:red)

    end

    fig[1:end,2] = Legend(fig, ax1, L"\sigma_M/\Omega_\mathrm{R}")
    fig
end
