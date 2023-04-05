function ideal_propagation(;Nmvals, Nc, ΩR, a, title=false)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 23, palette=(color=cgrad(:Dark2_7),))
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

function dis_propagation(;ΩR=0.1, σx=120, σMvals=[0.005, 0.01, 0.02, 0.04], xl=(0,2500), yl=(0,1100))

	fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)
	
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time (fs)", ylabel=L"d = \sqrt{\left \langle x^2 \right \rangle} / a", ylabelsize=30, xticks=0:500:2500)

    for σM in σMvals
        dis_propagation!(ax, ΩR=ΩR, σx=σx, σM=σM)
    end

	xlims!(ax, xl...)
	ylims!(ax, yl...)
    axislegend(ax, ax, L"\sigma_M / \Omega_R", nbanks=2, position=:lt)#, orientation=:horizontal)
    fig
end

function dis_propagation!(ax::Axis; ΩR=0.1, σx=120, σM=0.005)	    
    
	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 

    path = joinpath(@__DIR__, "../../propagation_study/disorder/Nm500_Nc500_a10_1/R$Rstr/sm$sm/out.h5")

    d = h5read(path, "$(Int(σx))_avg_d")
	tvals = h5read(path, "d_tvals") .* 1000

    lsty = σM ≥ ΩR ? :dash : :solid  
    lines!(ax, tvals, d, linewidth=3.0, linestyle=lsty, label=L"%$(σM)")
end
