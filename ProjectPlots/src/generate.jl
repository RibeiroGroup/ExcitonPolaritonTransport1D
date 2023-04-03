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