using Measurements
using Statistics
using PrettyTables

"""
    Plot the propagation (d x time) for several system sizes using different number of cavity modes.
"""
function ideal_propagation_conv(;Nmvals=[1000, 2000, 5000, 10000], Ncvals=[50, 100, 200, 400, 500, 800], a=10, ΩR=0.1)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    # Create figure object
    fig = Figure()

    # Create axis into the figure
    axs = [Axis(fig[i,1]) for i = eachindex(Nmvals)]

    # Axis configuration
    for (ax, Nm) in zip(axs, Nmvals)
        #ax.title  = L"N_M = %$Nm"
        ax.xticklabelsvisible = false
        ax.xticks = 0:100:500
        ax.yticks = [100, 300, 500]
        xlims!(ax, 0, 500)
        ylims!(ax, 0, 500)
    end
    axs[end].xticklabelsvisible = true
    axs[end].xlabel = "Time (fs)"

    # Common label to all y-axes
    Label(fig[1:end,0], L"d = \sqrt{\left\langle x^2 \right\rangle} / a", rotation = pi/2)

    for (i,Nm) in enumerate(Nmvals)
        for Nc in Ncvals

            # Fetch mean squared disp
            t, x2 = get_ideal_x2(Nm=Nm, Nc=Nc, a=a, ΩR=ΩR)

            # Compute d, note that a = 10
            d = (sqrt.(x2) .- sqrt(x2[1])) ./ 10

            # Produce lineplot on the correspoding axis
            lines!(axs[i], t .* 1000, d, label = L"%$(2*Nc+1)")
            text!(axs[i], 0, 1, text=L"N_M = %$Nm", align=(:left, :top), space=:relative)
        end
    end

    fig[1:3,2] = Legend(fig, axs[1], L"N_c")
    fig
end

"""
    Plots the error in d (with respect to a computation with 1601 cavity modes) as a function of the number of cavity modes 
    for several system sizes. No disorder.
"""
function ideal_error(;Nmvals=[5000, 8000, 10000, 15000, 20000], Ncvals=[50, 100, 200, 400, 500], Ncref=800, ΩR=0.1, a=10, tmax::Int=5000)

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

    ideal_error!(ax1, ax2, Nmvals=Nmvals, Ncvals=Ncvals, Ncref=Ncref, ΩR=ΩR, a=a, tmax=tmax)

    Legend(fig[1:end,2], ax1, L"N_\mathrm{M}") 
    axislegend(ax2)
    fig
end

function ideal_error!(ax1::Axis, ax2::Axis; Nmvals=[5000, 8000, 10000, 15000, 20000], Ncvals=[50, 100, 200, 400, 500], Ncref=800, ΩR=0.1, a=10, tmax::Int=5000)

    # Different markers to be used for each Nc value
    mkers = [:circle, :rect, :utriangle, :diamond, :star5, :cross, :xcros, :pentagon]

    # Slice of the simulation used to compute error
    # Using time step of 10 fs compute the tmax index
    slice = 1:(tmax ÷ 10 + 1)
    for (Nm, mk) in zip(Nmvals, mkers)

        # Get reference data - Nc = 1601 - Note that we are computing d
        _, x2 = get_ideal_x2(Nm=Nm, Nc=Ncref, ΩR=ΩR, a=a)
        ref = sqrt.(x2[slice]) ./ a

        # Preallocate arrays for cavity energies and error values
        err_vals = zeros(length(Ncvals))
        cav_e = zeros(length(Ncvals))
        for (j,Nc) in enumerate(Ncvals)

            # Fetch data for given Nc
            _, x2 = get_ideal_x2(Nm=Nm, Nc=Nc, ΩR=ΩR, a=a)
            d = sqrt.(x2[slice]) ./ a
            
            # Read output files to find the maximum cavity energy (i.e. energy cutoff)
            cav_e[j] = max_cavity_energy(Nm=Nm, Nc=Nc, a=a)

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
    lines!(ax2, 2:0.05:3.5, [0.8*exp(-9.5*(e-2)) for e in 2:0.05:3.5], color=:black, linestyle=:dash, label=L"\exp[-\alpha(E_\mathrm{cutoff}-\hbar\omega_0)]")
end

function ideal_Ecutoff()
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    # Create figure object
    fig = Figure()

    # Create axis into the figure
    #ax1 = Axis(fig[1,1], xlabel="Number of cavity modes", title="(a)", titlealign = :left)
    ax = Axis(fig[1,1], xlabel=L"\Omega_R \; \mathrm{(eV)}", titlealign = :left)
    xlims!(ax, 0, 0.35)

    ideal_Ecutoff!(ax, σx=60, Em=2.0)
    ideal_Ecutoff!(ax, σx=120, Em=2.0)
    ideal_Ecutoff!(ax, σx=180, Em=2.0)
    ideal_Ecutoff!(ax, σx=60, Em=2.2)
    ideal_Ecutoff!(ax, σx=120, Em=2.2)
    ideal_Ecutoff!(ax, σx=180, Em=2.2)

    lines!(ax, 0:0.01:0.4, [2.40 + 2*R for R = 0:0.01:0.4], linestyle=:dash, linewidth=3, color=:black, label=L"2.4 + 2\Omega_R")
    axislegend(ax)

    fig
end

function ideal_Ecutoff!(ax::Axis; σx=120, Em=2.0, maxerror=1, m, c)

    Estr = replace(string(Em), "."=>"p")
    path = joinpath(@__DIR__, "../../../mode_convergence/eband/Em$Estr/")
    Rvals = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
    Etgt = zeros(length(Rvals))

    for (i,ΩR) in enumerate(Rvals)

        R = replace(string(ΩR), "."=>"p")
        ref = h5read(path * "R$R/Nc800/out.h5", "sx$(σx)_d")

        # Preallocate arrays for cavity energies and error values
        Ncvals = [0, 45, 64, 79, 93, 105, 116, 127, 137, 147, 156, 165, 174, 183, 192, 201] 
        err_vals = zeros(length(Ncvals))
        cav_e = zeros(length(Ncvals))
        for (j,Nc) in enumerate(Ncvals)

            # Fetch data for given Nc
            d = h5read(path * "R$R/Nc$Nc/out.h5", "sx$(σx)_d")
            
            # Read output files to find the maximum cavity energy (i.e. energy cutoff)
            cav_e[j] = max_cavity_energy(Nm=5000, Nc=Nc, a=10, ideal=false, path=path*"R$R/Nc$Nc/log.out")

            # Compute error w.r.t to reference as the mean of the relative absolute error
            err_vals[j] = mean(abs.(d .- ref) ./ ref)
        end

        # Produce scatter and line plots
        #scatter!(ax, cav_e, err_vals, label=L"%$(ΩR)")
        idx = findfirst(x -> x < maxerror/100, err_vals)
        Etgt[i] = cav_e[idx]

        # Print out human readable table
        println("ΩR = $ΩR - Ecutoff for target accuracy = $(cav_e[idx])")
        pretty_table(hcat(Ncvals, cav_e, [round(100*e, digits=2) for e in err_vals]);
        header = ["Nc", "Energy (eV)", "Error"])
    end

    scatter!(ax, Rvals, Etgt, marker=m, color=c, markersize=15)
end

function plot_dis_error_and_propagation()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 23, palette=(color=cgrad(:Dark2_7),))
    set_theme!(fontsize_theme)

    # Create figure object
    fig = Figure()

    # Create axis for error plot
    ax1 = Axis(fig[2:4,3:4], xlabel="Cavity energy cutoff (eV)", ylabel="Error", xticks=2:0.1:2.4)
    ylims!(-0.01, 0.85)
    xlims!(1.995, 2.48)

    # Plot error
    plot_dis_error!(ax1)

    # Create axis for propagation plots
    ax2 = Axis(fig[1:2,1:2], xticks=[0,300,600,900], yticks=0:100:300, xticklabelsvisible=false)
    ax3 = Axis(fig[3:4,1:2], xticks=[0,300,600,900], yticks=0:100:300, xlabel="Time (fs)")
    Label(fig[1:4,0], L"d = \sqrt{\left \langle x^2 \right \rangle} / a", rotation=π/2)
    linkaxes!(ax2, ax3)
    xlims!(ax3, 0, 1000)
    ylims!(ax3, 0, 350)

    # Plot propagations
    plot_dis_propagation!(ax2, σM=0.02, Ncvals=[0, 10, 75, 100, 200, 400])#Ncvals=[0, 1, 5, 10, 20, 35, 50, 75, 100, 200, 400, 800])
    plot_dis_propagation!(ax3, σM=0.05, Ncvals=[0, 10, 75, 100, 200, 400])#Ncvals=[0, 1, 5, 10, 20, 35, 50, 75, 100, 200, 400, 800])

    # Legend for the propagation plots
    fig[1,3] = Legend(fig, ax2, L"N_c", nbanks=2)
    fig[1,4] = Legend(fig, ax1, L"\sigma_M/\Omega_R", nbanks=2)

    # Add plot label (a, b, c)
    for (l, ax) in zip(("(a)", "(b)", "(c)"), [ax2, ax3, ax1])
        text!(ax, 0.05, 0.99, text=l, space=:relative, align=(:left, :top), font=:bold)
    end
    text!(ax2, 0.98, 0.03, text=L"\sigma_M/\Omega_R = 0.2", space=:relative, align=(:right, :bottom), font=:bold)
    text!(ax3, 0.98, 0.03, text=L"\sigma_M/\Omega_R = 0.5", space=:relative, align=(:right, :bottom), font=:bold)

    trim!(fig.layout)
    fig
end
 
function plot_dis_propagation(; σM, ΩR=0.1, σx=60, NR=100, Ncvals=[0, 5, 50, 75, 100, 800], tmax=1000)
    # Plot global settings
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    # Create figure object
    fig = Figure()

    # Create axis into the figure
    ax = Axis(fig[1,1], xlabel="Time (fs)", ylabel=L"d = \sqrt{\left \langle x^2 \right \rangle} / a")
    xlims!(ax, 0,tmax+20)

    # Plot onto axis
    plot_dis_propagation!(ax, σM=σM, ΩR=ΩR, σx=σx, NR=NR, Ncvals=Ncvals, tmax=tmax)

    fig[2,1] = Legend(fig, ax, L"N_c", orientation=:horizontal, nbanks=2)
    fig
end

function plot_dis_propagation!(ax::Axis; σM, ΩR=0.1, σx=60, NR=100, ωM=2.0, Ncvals=[0, 5, 50, 75, 100, 800], threshold=0.05, tmax=1000)


    Rstr = "R" * replace(string(ΩR), "."=>"p")
    σMstr = "sm" * replace(string(σM), '.'=>'p')
    ωMstr = "Em" * replace(string(ωM), '.'=>'p')

    tfinal = tmax ÷ 10 + 1
    slice = 1:tfinal
    tvals = (0:10:5000)[slice]

    ref_path = joinpath(@__DIR__, "../../../mode_convergence/$ωMstr/$Rstr/$σMstr/Nc800/out.h5")
    dref = h5read(ref_path, "NR_$(NR)_sm$(σx)_avg_d")[slice]
    σ = h5read(ref_path, "NR_$(NR)_sm$(σx)_std_d")[slice]
    band!(ax, tvals, dref .- σ, dref .+ σ, color=:ivory3)

    for Nc in Ncvals

        path = joinpath(@__DIR__, "../../../mode_convergence/$ωMstr/$Rstr/$σMstr/Nc$Nc/out.h5")
        d = h5read(path, "NR_$(NR)_sm$(σx)_avg_d")[slice]
        max_cav_e = max_cavity_energy(path=replace(path, "out.h5"=>"log.out"), ideal=false, Nm=5000, Nc=Nc, a=10)
        error = sum(abs.(d - dref) ./ dref) / length(1:tfinal)
        println("Nc = $Nc Error = $error   Ecutoff = $max_cav_e")

        lsty = error > threshold ? :dot : :solid

        lines!(ax, tvals, d, label = L"%$(2*Nc+1)", linewidth=3, linestyle=lsty)

    end
end

function plot_dis_error(;NR = 100, σx = 60, Ncvals=[0, 1, 5, 10, 20, 35, 50, 75, 100], tmax=1000)
    # Plot global settings
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    # Create figure object
    fig = Figure()

    # Create axis
    ax1 = Axis(fig[1,1], title="(a)", titlealign = :left)
    ax2 = Axis(fig[2,1], title="(b)", titlealign = :left, xlabel="Cavity energy cutoff (eV)")

    # Label ΩR
    text!(ax1, 0.7, 0.8, text=L"\Omega_R = 0.1\;\mathrm{eV}", align=(:left, :top), space=:relative)
    text!(ax2, 0.7, 0.8, text=L"\Omega_R = 0.2\;\mathrm{eV}", align=(:left, :top), space=:relative)

    # Plot onto axes
    plot_dis_error!(ax1, ΩR=0.1, NR=NR, σx=σx, Ncvals=Ncvals, tmax=tmax)
    plot_dis_error!(ax2, ΩR=0.2, NR=NR, σx=σx, Ncvals=Ncvals, tmax=tmax)

    # Create axis into the figure
    linkyaxes!(ax1, ax2)
    linkxaxes!(ax1, ax2)
    Label(fig[1:end,0], "Error", rotation=pi/2)

    # Create legend
    fig[1:end,2] = Legend(fig, ax1, L"\sigma_M/\Omega_\mathrm{R}")

    fig
end

function plot_dis_error!(ax::Axis; ΩR=0.1, ωM=2.0, NR = 100, σx = 60, Ncvals=[0, 1, 5, 10, 20, 35, 50, 75, 100], tmax=1000)

    # Preallocate arrays for error values, error deviations and cavity energies
    err_vals = zeros(length(Ncvals))
    err_devs = zeros(length(Ncvals))
    cav_e = zeros(length(Ncvals))

    # Slice on time steps - this controls whether we average early or late propagation times
    tfinal = tmax ÷ 10 + 1
    slice = 1:tfinal

    # Different markers to be used for each Nc value
    mkers = [:circle, :rect, :utriangle, :diamond]
    clrs = [:royalblue2, :darkolivegreen, :darkgoldenrod1 ,:red3]
    rats = [0.05, 0.2, 0.5, 1]

    # Loop through Axes and corresponding Rabi splittings 

    Rstr = "R" * replace(string(ΩR), "."=>"p")
    Emstr = "Em" * replace(string(ωM), '.'=>'p')

    for i in eachindex(rats)

        mk = mkers[i]
        ratio = rats[i]

        σM = round(ΩR * ratio, digits=3)
        σMstr = "sm" * replace(string(σM), '.'=>'p')

        # Get data for Nc = 1601V
        ref_path = joinpath(@__DIR__, "../../../mode_convergence/$Emstr/$Rstr/$σMstr/Nc800/out.h5")
        dref = h5read(ref_path, "NR_$(NR)_sm$(σx)_avg_d")[slice]
        std_dref = h5read(ref_path, "NR_$(NR)_sm$(σx)_std_d")[slice]

        # Create Measurement objects, with the standard deviation representing the error of d
        ref_mmt = [measurement(dref[i], 2 .*std_dref[i]) for i = eachindex(dref)]

        for (j,Nc) in enumerate(Ncvals)

            # Path to the output data
            path = joinpath(@__DIR__,"../../../mode_convergence/$Emstr/$Rstr/$σMstr/Nc$Nc/out.h5")

            # Fetch d values along with their standard deviation
            d = h5read(path, "NR_$(NR)_sm$(σx)_avg_d")[slice]
            std = h5read(path, "NR_$(NR)_sm$(σx)_std_d")[slice]

            # Create Measurement objects, with the standard deviation representing the error of d
            d_mmt = [measurement(d[i], 2 .* std[i]) for i = eachindex(d)]

            # Read output files to find the maximum cavity energy (i.e. energy cutoff)
            cav_e[j] = max_cavity_energy(path=replace(path, "out.h5"=>"log.out"), ideal=false, Nm=5000, Nc=Nc, a=10)

            # Compute error and associated uncertainty 
            err = mean(abs.(d_mmt .- ref_mmt) ./ ref_mmt)
            err_vals[j] = err.val
            err_devs[j] = err.err
        end

        # Produce plots
        scatter!(ax, cav_e, err_vals, label=L"%$ratio", marker=mk, color=clrs[i])
        lines!(ax, cav_e, err_vals, color=clrs[i])
        errorbars!(ax, cav_e, err_vals, err_devs, color=clrs[i])
    end

    # Plot no disorder for reference
    ideal_ref = h5read(joinpath(@__DIR__, "../../../mode_convergence/$Emstr/$Rstr/sm0/Nc800/out.h5"), "sx$(σx)_d")
    for (j,Nc) in enumerate(Ncvals)
        d_ideal = h5read(joinpath(@__DIR__, "../../../mode_convergence/$Emstr/$Rstr/sm0/Nc$Nc/out.h5"), "sx$(σx)_d")
        err_vals[j] = mean(abs.(d_ideal .- ideal_ref) ./ ideal_ref)
    end
    lines!(ax, cav_e, err_vals, linestyle=:dot, color=:black, linewidth=3)
end