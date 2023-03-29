function ideal_propagation(;Nmvals, Nc, ΩR, a, title=false)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel="Time (fs)", xticks=0:200:900)
    xlims!(ax1, 0,1000)
    ylims!(ax1, 0, 800)
    ax2 = Axis(fig[2,1], xlabel="Time (ps)", xticks = 0:5:50)
    xlims!(ax2, 0,30)
    Label(fig[1:end,0], L"d = \sqrt{\left \langle x^2 \right \rangle} / a", rotation=π/2, justification=:center, fontsize=25)

    for (i, Nm) in enumerate(Nmvals)

        clr = i != 5 ? Cycled(i) : Cycled(i+1)

        # Get ⟨x²⟩ data using an auxiliary function
        t, x2 = get_ideal_x2(Nm=Nm, Nc=Nc, ΩR=ΩR, a=a)

        # Compute d from x2
        # Shift results so they all start from zero.
        d = (sqrt.(x2) .- sqrt(x2[1])) ./ a

        # Plot on both panels 
        lines!(ax1, t.*1000, d, label = L"N_M = %$Nm", color=clr)
        lines!(ax2, t, d, label = L"N_M = %$Nm", color=clr)
    end

    if title
        ax1.title = L"N_c = %$(2*Nc+1)\;\; \Omega_R = %$(ΩR)\;\mathrm{eV}\;\; a = %$a \; \mathrm{nm}"
    end

    axislegend(ax1, position=:lt)
    fig
end

function SI_Nm1000()

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