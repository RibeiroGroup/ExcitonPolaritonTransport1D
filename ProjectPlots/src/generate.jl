using CairoMakie

"""
    Generate all articles figures in the given path.
"""
function generate_figures(;path=".", quality=5, SI=false, maxnum=100)

    # Scan module for all functions named figX, where X ∈ N
    println("Looking for plot functions 👀")
    plot_functions = []
    fstr = SI ? "SI_fig" : "fig"
    i = 1
    while i ≤ maxnum
        if isdefined(ProjectPlots, Symbol(fstr*"$i"))
            push!(plot_functions, eval(Symbol(fstr*"$i")))
        elseif i > 1 # Make sure to check for at least fig2
            break
        end
        i += 1
    end

    nf = length(plot_functions)
    println("$nf functions found!\n")

    # Call plot functions and save figures
    println("📜 Generating article figures...")
    for (i, func) in enumerate(plot_functions)
        print("($i/$nf) $(string(func)) ")
        save(joinpath(path, string(func)*".png"), func(), px_per_unit=quality)
        println("✅")
    end
    println("Done 🎉")
end

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 1601, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
fig2() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=800, ΩR=0.1, a=10)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 5 ps of simulation.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
fig3() = ideal_error()

"""
(a) Propagation under disorder for σM = 0.02 eV and several Nc vals.
(b) Propagation under disorder for σM = 0.05 eV and several Nc vals.
Band plot is the reference trajectory (Nc=1601)
(c) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
for several values of disorder.
Error computed over 1 ps of simulation.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
"""
function fig4()

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

"""
For a given detuning value, produces a graph with two columns. The first column being q0 = 0 (that is, average exciton)
and the second q0 != 0. Each row then shows the mode weight for a different initial state spread (σx) at different 
Rabi splitting values. Dashed lines are used for modes with q ≤ 0, whereas solid lines are for q > 0.
Red dotted lines indicates where the photon is energy resonant with the molecules, whereas the blue dotted lineheight
indicates the photon mode that matches the average momentum of the initial state.
No disorder. Weights computed over 5 ps.
Nm = 5000, a = 10 nm,
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
"""
function fig5(;δ=-0.2)

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

    axs[1,1].title = L"\bar{q}_0 = 0"
    axs[1,2].title = L"\bar{q}_0 \neq 0"
    Label(fig[1:end,0], "Relative weight", justification=:center, rotation=pi/2)
    Label(fig[4,1:2], "Mode energy (eV)", justification=:center)

    for (j,q0) = enumerate([true, false])
        for (i,σx) in enumerate([60, 120, 180])
            plot_ideal_mw!(axs[i,j], σx=σx, Em=2-δ, zeroq=q0)
            text!(axs[i,j], 0.95, 0.95, text=titles[i,j], align=(:right, :top), space=:relative, font=:bold)
            text!(axs[i,j], 0.99, 0.7, text=L"\sigma_x / a = %$(Int(σx/10))", align=(:right, :top), space=:relative, fontsize=20)
            vlines!(axs[i,j], [2 - δ], color = :red, linestyle=:dot)
            if !q0
                vlines!(axs[i,j], [2.1], color = :blue, linestyle=:dot)
            end
        end
    end

    # Manualy creates a legend... a bit convoluted
    Ωmarkers = [LineElement(color=Makie.wong_colors()[i]) for i = 1:4]
    Ωvals = [L"0.05", L"0.1", L"0.2", L"0.3"]

    qmarkers = [LineElement(color=:black, linestyle=s) for s = [:solid, :dash]]
    qvals = [L"q \; > \; 0", L"q \leq 0"]

    vmarkers = [LineElement(color =:red, linestyle=:dot), LineElement(color=:blue, linestyle=:dot)]
    vvals = [L"\hbar\omega_q = E_M", L"q = \bar{q}_0"]
    fig[5,1:2] = Legend(fig, [Ωmarkers, qmarkers, vmarkers], [Ωvals, qvals, vvals], [L"\mathbf{\Omega_R\:\mathrm{(eV)}}", "Momentum", "Resonance"], orientation=:horizontal)
    fig
end

"""
For Em = 2.2 (δ = -0.2), produces plots for the mode weight of a disordered system (σM = 0.02)
Band plots cover twice the standard deviation in each case. Dotted lines show energy (red) and momentum (blue)
resonances. First and second column show results colected over a short and long period of times, respectively.  
Nm = 5000, a = 10 nm, σx = 120 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
50 Realizations.
"""
function fig6()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    ΩRvals = [0.05, 0.1, 0.2]

    fig = Figure()
    axs = [Axis(fig[i,j], yticklabelsvisible=(j==1), xticks=2:0.1:2.4, xticklabelsvisible=(i==length(ΩRvals))) 
    for i=eachindex(ΩRvals), j=1:2]

    for j = 1:2
        for i = eachindex(ΩRvals)
            plot_dis_mw!(axs[i,j], ΩR=ΩRvals[i], early=(j==1), zeroq=false,  Em=2.2, σx=120, σM=0.02)
            xlims!(axs[i,j], 2, 2.45)
            ylims!(axs[i,j], -0.1, 1.2)
            vlines!(axs[i,j], [2.2], color = :red, linestyle=:dot)
            vlines!(axs[i,j], [2.1], color = :blue, linestyle=:dot)
            text!(axs[i,j], 2.36, 1.0, align=(:center, :center), text=L"\Omega_R = %$(ΩRvals[i])\; \mathrm{eV}", fontsize=22)
            text!(axs[i,j], 2.04, 1.0, align=(:center, :center), text="($('`'+i+3*(j-1)))", fontsize=22, font=:bold)
        end
    end
    axs[1,1].title = "Early (500 fs)"
    axs[1,2].title = "Late (5000 fs)"

    Label(fig[end+1, 1:2], "Mode energy (eV)")
    Label(fig[1:end, 0], "Average relative weight", rotation=π/2)

    # Manualy creates a legend... a bit convoluted
    mkers = [:rtriangle, :ltriangle]
    qmarkers = [MarkerElement(color=Makie.wong_colors()[i], marker=mkers[i]) for i=1:2]
    qvals = [L"q \; > \; 0", L"q \leq 0"]

    vmarkers = [LineElement(color =:red, linestyle=:dot), LineElement(color=:blue, linestyle=:dot)]
    vvals = [L"\hbar\omega_q = E_M", L"q = \bar{q}_0"]
    fig[end+1,1:2] = Legend(fig, [qmarkers, vmarkers], [qvals, vvals], ["Momentum", "Resonance"], orientation=:horizontal, titlesize=18)#, titleposition=:left)
    fig
end

"""
Bar plots of the wavepacket propagating in time with and without photon with negative momentum along x.
Each bar groups the position of 50 adjacent molecules. Error bars are taken as one standard deviation.
Nm = 5000, a = 10 nm, σx = 120 nm, Em = 2.0 (δ=0), ΩR = 0.1 eV
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
100 Realizations
"""
function fig7()

    # Plot parameters
    σ=0.04
    times = [0, 300, 600, 900]
    xmin = 45
    xmax = 55

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)
    
    fig = Figure()
    
    # Create axis
    axs = [Axis(fig[i,j], xticks = (1:100, ["[$(i+1)-$(i+50)]" for i = 0:50:4950]), xticklabelrotation = pi/4) for i = 1:4, j = 1:2]
    
    ## Set axis limits
    for ax in axs
        xlims!(ax, xmin, xmax)
        ax.xticklabelsize=12
    end
    
    # Fetch data
    fname = "Nm5000_Nc200_R0p1_a10_1"*"_Em" *replace(string(σ), "."=>"p")
    
    vals = zeros(100,100)
    pvals = zeros(100,100)
    for i = 1:4
    
        # Get data index from time step
        ti = times[i] ÷ 10 + 1
    
        # Loop through realizations and build an array with wave packets at the same time step but different realizations
        for n = 1:100
            vals[:,n]  .= h5read(joinpath(@__DIR__, "../../momentum_study/all/$fname/out$n.h5"), "120_wavepacket_bars", (:,ti))
            pvals[:,n] .= h5read(joinpath(@__DIR__, "../../momentum_study/pos/$fname/out$n.h5"), "120_wavepacket_bars", (:,ti))
        end
    
        # Average over realizations
        avg =  [mean(vals[n,:]) for n = 1:100]
        pavg = [mean(pvals[n,:]) for n = 1:100]
    
        # Standard deviation of the realizations
        dev = [std(vals[n,:]) for n = axes(vals,2)]
        pdev = [std(pvals[n,:]) for n = axes(pvals,2)]
    
        println(times[i])
        println(sum(avg[52:end]))
        println(sum(pavg[52:end]))
        println("\n")
    
        # Plot wavepackets 
        barplot!(axs[i,1], avg)
        barplot!(axs[i,2], pavg, color=Cycled(3))
    
        # Plo error bars as one standard deviation
        errorbars!(axs[i,1], 1:100, avg, dev  , color = :red4, whiskerwidth = 10)
        errorbars!(axs[i,2], 1:100, pavg, pdev, color = :red4, whiskerwidth = 10)
    
        # Add text indicating time step
        text!(axs[i,1], L"t = %$(times[i]) \; \mathrm{fs}", space = :relative, position = Point2f(0.05, 0.65))
    
        linkaxes!(axs[i,1], axs[i,2])
    
        # Hide redundent axis labels for clarity
        axs[i,2].yticklabelsvisible = false
        if i != 4
            axs[i,1].xticklabelsvisible = false
            axs[i,2].xticklabelsvisible = false
        end
    end
    
    # Add global axis label
    Label(fig[5, 1:2], "Molecule index", fontsize = 20)
    Label(fig[1:4, 0], L"|\left \langle n | \Psi \right \rangle|^2", rotation = pi/2, justification=:right, fontsize = 25)
    
    axs[1,1].title = L"N_c = 401"
    axs[1,2].title = L"N_c = 201\;(q\;>\;0\;\mathrm{only})"
    axs[1,1].titlesize = 20
    axs[1,2].titlesize = 20
    
    fig
end

"""
Disordered propagation (d x time) at early and late times for ΩR = 0.05 and 0.1 eV and several
disorder values. Zoom-in highlights the crossover due to disorder enhanced transport.
Nm = 5000, a = 10 nm, σx = 120 nm, Em = 2.0 (δ=0), 
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
100 Realizations
"""
function fig8()

    σx = 120
    σMvals = [0.01, 0.02, 0.04, 0.1, 0.3, 0.5]

	fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)
	
    fig = Figure()

    ax1 = Axis(fig[1,1], xticks=0:600:2000, title = L"\Omega_R = 0.05\;\mathrm{eV}")
	xlims!(ax1, 0, 2500)
	ylims!(ax1, 0, 800)

    for σM in σMvals
        dis_propagation!(ax1, ΩR=0.05, σx=σx, σM=σM)
    end

    ax2 = Axis(fig[1,2], yticklabelsvisible=false, title = L"\Omega_R = 0.1\;\mathrm{eV}", xticks=0:600:2000)
	xlims!(ax2, 0, 2500)
	ylims!(ax2, 0, 800)

    for σM in σMvals
        dis_propagation!(ax2, ΩR=0.1, σx=σx, σM=σM)
    end

    ax3 = Axis(fig[2,1], xticks=0:100:400)
	xlims!(ax3, 0, 500)
	ylims!(ax3, 0, 300)
    for σM in σMvals
        dis_propagation!(ax3, ΩR=0.05, σx=σx, σM=σM)
    end
    vlines!(ax3, [250], linestyle=:dot, linewidth=3, color=:black)

    ax4 = Axis(fig[2,2], yticklabelsvisible=false, xticks=0:100:400)
	xlims!(ax4, 0, 500)
	ylims!(ax4, 0, 300)
    for σM in σMvals
        dis_propagation!(ax4, ΩR=0.1, σx=σx, σM=σM)
    end
    vlines!(ax4, [100], linestyle=:dot, linewidth=3, color=:black)

    # Label axes
    Label(fig[1:end,0], L"d = \sqrt{\left \langle x^2 \right \rangle} / a", rotation=π/2, justification=:center, fontsize=25)
    Label(fig[end+1, 1:2], "Time (fs)", fontsize=25)

    # Create legend
    Legend(fig[1:2, 3], ax1, L"\sigma_M\;\mathrm{(eV)}")

    # Add graph labels (a,b,c.. etc)
    for (i,ax) in enumerate([ax1, ax3, ax2, ax4])
        text!(ax, 0.03, 1.0, align=(:left, :top), text="($('`'+i))", fontsize=22, font=:bold, space=:relative)
    end
    fig
end

"""
Upper panel: Exciton amplitude at the system boundaries (first and final 100 molecules) as a function of time.
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nm = 1000, Nc = 1601, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm.
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
function SI_fig1()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure()

    # Create plot for probability of finding exciton at the systems edge
    ax1 = Axis(fig[1,1], xlabel="Time (fs)", ylabel="Molecular amplitude\nat the boundaries", ylabelsize=20) 
    xlims!(ax1, 0, 405)
    pend = h5read(joinpath(@__DIR__, "../../propagation_study/size_effects/out.h5"), "pend")
    lines!(ax1, 0:1:400, pend[1:401])
    vlines!(ax1, [75], linestyle=:dash, color=:black)
    text!(ax1, 75, 0.03, text="Tails reach\nsystem boundaries", fontsize=15, justification=:center, align=(:center, :top))

    # Create plot for probability of finding photon
    ax2 = Axis(fig[2,1], xlabel="Time (fs)", ylabel=L"P_\mathrm{phot}")
    xlims!(ax2, 0, 205)
    pphot = h5read(joinpath(@__DIR__, "../../propagation_study/size_effects/out.h5"), "pphot")
    vlines!(ax2, [18.5, 59.5], linestyle=:dot, color=:black)
    scatter!(ax2, 0:1:200, pphot[1:201], color=:lightgoldenrod2)
    text!(ax2, 13, 0.18, text=L"T = 41\; \mathrm{fs} \approx \Omega_R^{-1}")

    fig
end

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 1601, ΩR = 0.1 eV, a = 20 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig2() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=800, ΩR=0.1, a=20, title=true)

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 1601, ΩR = 0.2 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig3() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=800, ΩR=0.2, a=10, title=true)

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 1, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig4() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=0,   ΩR=0.1, a=10, title=true)

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 201, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig5() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=100, ΩR=0.1, a=10, title=true)

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 801, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig6() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=400, ΩR=0.1, a=10, title=true)

"""
Propagations using various system sizes (Nm) and number of cavity modes (Nc).
ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig7() = ideal_propagation_conv(ΩR=0.1, a=10)

"""
Propagations using various system sizes (Nm) and number of cavity modes (Nc).
ΩR = 0.2 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig8() = ideal_propagation_conv(ΩR=0.2, a=10)

"""
Propagations using various system sizes (Nm) and number of cavity modes (Nc).
ΩR = 0.1 eV, a = 20 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig9() = ideal_propagation_conv(ΩR=0.1, a=20)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 0.5 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig10() = ideal_error(tmax=500)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 1 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig11() = ideal_error(tmax=1000)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 20 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig12() = ideal_error(tmax=20000)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 5 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.05 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig13() = ideal_error(ΩR=0.05)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 5 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.3 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig14() = ideal_error(ΩR=0.3)