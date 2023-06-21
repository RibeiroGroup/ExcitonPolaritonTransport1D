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
    #ax1 = Axis(fig[1,1], xlabel="Time (fs)", ylabel="Molecular amplitude\nat the boundaries", ylabelsize=20) 
    #xlims!(ax1, 0, 405)
    #pend = h5read(joinpath(@__DIR__, "../../propagation_study/size_effects/out.h5"), "pend")
    #lines!(ax1, 0:1:400, pend[1:401])
    #vlines!(ax1, [75], linestyle=:dash, color=:black)
    #text!(ax1, 75, 0.03, text="Tails reach\nsystem boundaries", fontsize=15, justification=:center, align=(:center, :top))

    # Create plot for probability of finding photon
    ax2 = Axis(fig[1,1], xlabel="Time (fs)", ylabel=L"P_\mathrm{phot}")
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
SI_fig2() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=800, ΩR=0.1, a=20)

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 1601, ΩR = 0.2 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig3() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=800, ΩR=0.2, a=10)

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 1, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig4() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=0,   ΩR=0.1, a=10)

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 201, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig5() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=100, ΩR=0.1, a=10)

"""
Upper panel: Short time (up to 1 ps) wavepacket width (d) for several system sizes (Nm).
Lower panel: Long time (up to 30 ps) wavepacket width (d) for several system sizes (Nm).
Nc = 801, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig6() = ideal_propagation(Nmvals=[1000, 5000, 10000, 15000, 20000], Nc=400, ΩR=0.1, a=10)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 0.5 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig7() = ideal_error(tmax=500)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 1 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig8() = ideal_error(tmax=1000)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 20 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig9() = ideal_error(tmax=20000)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 5 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.05 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig10() = ideal_error(ΩR=0.05)

"""
(a) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Nc
(b) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
Error computed over 5 ps of simulation.
Exponential profile shown is the same as in Fig. 3.
Nm = 5000, ΩR = 0.3 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
SI_fig11() = ideal_error(ΩR=0.3)

"""
(a) Propagation under disorder for σM = 0.02 eV and several Nc vals.
(b) Propagation under disorder for σM = 0.05 eV and several Nc vals.
Band plot is the reference trajectory (Nc=1601)
(c) Error due to cavity modes truncation (w.r.t to Nc = 1601) as a function of Ecutoff
for several values of disorder.
Error computed over 1 ps of simulation.
Nm = 5000, ΩR = 0.2 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
CHANGE FROM FIG 4: ΩR + 0.1
"""
SI_fig12() = fig4(ΩR=0.2, σM1=0.04, σM2=0.1, TrajNcs=[0, 10, 75, 100], dymax=450)
SI_fig13() = fig4(σx=180, dymax=450)
SI_fig14() = fig4(ωM=2.2, dymax=450)

"""
Photon probability over time for various disoder and initial states.
Nm = 5000, ΩR = 0.1 eV, a = 10 nm, ωM = 2.0 eV, σx = 60 nm
Ly = 200 nm, Lz = 400 nm, Lx = Nm*a, ϵ = 3, nz = ny = 1.
No disorder.
"""
function SI_fig15()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 23, palette=(color=cgrad(:Dark2_7),))
    set_theme!(fontsize_theme)

    # Create figure object
    σMvals = [0.005 0.02; 0.05 0.1]
    fig = Figure()
    axs = [Axis(fig[i,j]) for i = 1:2, j=1:2]
    linkaxes!(axs...)
    tvals = 0:5:5000

    for i = 1:2, j=1:2
        σM = σMvals[i,j]
        σMstr = replace(string(σM), "."=>"p")
        ax = axs[i,j]
        xlims!(ax, 0, 300)
        ax.title = L"\sigma_M = %$(σM)\; \mathrm{eV}"
        if j == 2
            ax.yticklabelsvisible = false
        end
        if i == 1
            ax.xticklabelsvisible = false
        end
        for (k,σx) in enumerate([60, 120, 180])
            p0 = h5read(joinpath(@__DIR__, "../../mode_weight/disorder/Em2p0/R0p1/sm$σMstr/out.h5"), "q0_$(σx)_phot_cont")
            px = h5read(joinpath(@__DIR__, "../../mode_weight/disorder/Em2p0/R0p1/sm$σMstr/out.h5"), "nzq_$(σx)_phot_cont")
            lines!(ax, tvals, p0, label=L"\sigma_x = %$σx\;nm\;\bar{q}_0 = 0", linewidth=2, color=Makie.wong_colors()[k])
            scatter!(ax, tvals, p0, label=L"\sigma_x = %$σx\;nm\;\bar{q}_0 = 0", marker=:circle, color=Makie.wong_colors()[k])
            lines!(ax, tvals, px, label=L"\sigma_x = %$σx\;nm\;\bar{q}_0 \neq 0",linewidth=2, color=Makie.wong_colors()[k])
            scatter!(ax, tvals, px, label=L"\sigma_x = %$σx\;nm\;\bar{q}_0 \neq 0", marker=:diamond, color=Makie.wong_colors()[k])
        end
    end

    Label(fig[1:2,0], "Total photon probability", rotation=π/2)
    Label(fig[3,1:2], "Time (fs)")

    Legend(fig[4, 1:2], [[MarkerElement(marker=m, color=:black) for m = (:circle, :diamond)],
    [LineElement(color=c, linewidth=2) for c in Makie.wong_colors()[1:3]]],
    [[L"\bar{q}_0 = 0", L"\bar{q}_0\;>\;0"], string.([60, 120, 180])],
    ["", L"\sigma_x\;\mathrm{(nm)}"], orientation=:horizontal, titleposition=:left)

    fig
end

"""
Average photon probability for various σx and \bar{q}_0 as a function of disorder.
"""
function SI_fig16()

    # Plot global settings
    fontsize_theme = Theme(fontsize = 23, palette=(color=cgrad(:Dark2_7),))
    set_theme!(fontsize_theme)

    # Create figure object
    σxvals = [60, 120, 180]
    σMvals = [0.005, 0.02, 0.05, 0.1]
    ΩRvals = [0.05, 0.1, 0.2, 0.3]
    lbs = ["(a)" "(d)"; "(b)" "(e)"; "(c)" "(f)"]
    fig = Figure()
    axs = [Axis(fig[i,j], xticks=[0.005, 0.05, 0.1]) for i = 1:3, j=1:2]
    linkaxes!(axs...)

    for i = 1:3, j=1:2
        glab = lbs[i,j]
        σx = σxvals[i]
        ax = axs[i,j]
        xlims!(ax, 0, 0.11)
        if j == 2
            ax.yticklabelsvisible = false
        end
        if i == 1
            ax.xticklabelsvisible = false
            if j == 1
                ax.title = L"\bar{q}_0 = 0"
            else
                ax.title = L"\bar{q}_0\;>\;0"
            end
        elseif i < 3
            ax.xticklabelsvisible = false
        end
        q0 = j == 1 ? "q0_" : "nzq_" 
        ls = [:solid, :dash]
        mks = [:circle, :cross]
        for (n,Em) in enumerate([2.0, 2.2])
            Estr = "Em" * replace(string(Em), "."=>"p")
            for (k,ΩR) in enumerate(ΩRvals)
                Rstr = "R" * replace(string(ΩR), "."=>"p")
                pvals = zeros(4)
                for (l,σM) in enumerate(σMvals)
                    σMstr = replace(string(σM), "."=>"p")
                    p = h5read(joinpath(@__DIR__, "../../mode_weight/disorder/$Estr/$Rstr/sm$σMstr/out.h5"), "$(q0)$(σx)_phot_cont")
                    pvals[l] = mean(p)
                end
                scatter!(ax, σMvals, pvals, color=Cycled(k), marker=mks[n])
                lines!(ax, σMvals, pvals, color=Cycled(k), linestyle=ls[n])
                text!(ax, 0.1, 0.4, text=L"\sigma_x = %$(σx)\;\mathrm{nm}", align=(:right, :top))
                text!(ax, 0.1, 0.5, text=glab, align=(:right, :top))
            end
        end
    end

    Label(fig[1:3,0], "Average photon probability", rotation=π/2)
    Label(fig[4,1:2], L"\sigma_M\;\mathrm{(eV)}")


    # Manually creating legend, bit of a mess don't even try
    elem_1 = [LineElement(color = :black, linestyle = nothing),MarkerElement(color = :black, marker = :circle, markersize = 10,
    strokecolor = :black)]
    elem_2 = [LineElement(color = :black, linestyle = :dash),MarkerElement(color = :black, marker = :cross, markersize = 10,
    strokecolor = :black)]
    Legend(fig[1:3, 3], [[elem_1, elem_2],
    [PolyElement(color=c, linewidth=2) for c in cgrad(:Dark2_7)[1:4]]],
    [[L"2.0", L"2.2"], [L"0.05", L"0.1", L"0.2", L"0.3"]],
    [L"\omega_M\;\mathrm{(eV)}", L"\Omega_R\;\mathrm{(eV)}"])#, orientation=:horizontal)#, titleposition=:left)

    fig
end

"""
    #fig5 but with different detuning values
"""
SI_fig17() = fig5(δ=0.2)
SI_fig18() = fig5(δ=0.0)

SI_fig19() = fig6(σMvals=[0.005, 0.1])
SI_fig20() = fig6(early=true)
SI_fig21() = fig6(zeroq=true)
SI_fig22() = fig6(σx=60)
SI_fig23() = fig6(σx=180)
SI_fig24() = fig6(Em=2.0)