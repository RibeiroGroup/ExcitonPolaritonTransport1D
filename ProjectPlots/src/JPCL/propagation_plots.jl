function ideal_propagation(;Nmvals, Nc, ΩR, a, title=false)

    # Plot global settings
    fontsize_theme = Theme(fontsize = 23) 
    set_theme!(fontsize_theme)

    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel="Time (fs)", xticks=0:200:900)
    xlims!(ax1, 0,1000)
    ylims!(ax1, 0, 800)
    ax2 = Axis(fig[2,1], xlabel="Time (ps)", xticks = 0:5:50)
    xlims!(ax2, 0,15)
    Label(fig[1:end,0], L"d = \sqrt{\left \langle x^2 \right \rangle} / a", rotation=π/2, justification=:center, fontsize=25)

    ideal_propagation!(ax1, ax2, Nmvals=Nmvals, Nc=Nc, ΩR=ΩR, a=a)

    if title
        ax1.title = L"N_c = %$(2*Nc+1)\;\; \Omega_R = %$(ΩR)\;\mathrm{eV}\;\; a = %$a \; \mathrm{nm}"
    end

    axislegend(ax1, ax1, L"N_M", position=:lt, nbanks=2)
    fig
end

function ideal_propagation!(ax1, ax2; Nmvals, Nc, ΩR, a)
    for (i, Nm) in enumerate(Nmvals)

        clr = i != 5 ? Cycled(i) : Cycled(i+1)

        # Get ⟨x²⟩ data using an auxiliary function
        t, x2 = get_ideal_x2(Nm=Nm, Nc=Nc, ΩR=ΩR, a=a)

        # Compute d from x2
        # Shift results so they all start from zero.
        d = (sqrt.(x2) .- sqrt(x2[1])) ./ a

        # Plot on both panels 
        lines!(ax1, t.*1000, d, label = L"%$Nm", color=clr, linewidth=3)
        lines!(ax2, t, d, label = L"%$Nm", color=clr)
    end
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

    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    d = h5read(path, "$(Int(σx))_avg_d")
	tvals = h5read(path, "d_tvals") .* 1000

    lsty = σM ≥ ΩR ? :dash : :solid  
    lines!(ax, tvals, d, linewidth=3.0, linestyle=lsty, label=L"%$(σM)")
end

function phot_cont(σx=120, ΩRvals = [0.05, 0.1, 0.2, 0.3], smvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5])
    fig = Figure()
    axs = [Axis(fig[1,1]), Axis(fig[1,2]), Axis(fig[2,1]), Axis(fig[2,2])]

    i = 1
    for ΩR in ΩRvals
        Rstr = "R" * replace(string(ΩR), "."=>"p")
        for sm in smvals
            smstr = "sm" * replace(string(sm), '.'=>'p')
            path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/$Rstr/$smstr/out.h5")
            zpm = h5read(path, "$(σx)_avg_zpm")
            ppm = h5read(path, "$(σx)_avg_ppm")
            npm = h5read(path, "$(σx)_avg_npm")
            tvals = h5read(path, "wvp_tvals")
            total = [zpm[i] + sum(ppm[:,i]) + sum(npm[:,i]) for i = eachindex(tvals)]
            println("max Pphot ΩR = $ΩR eV, σM = $sm eV  =  $(maximum(total))")
            scatter!(axs[i], tvals, total, label = L"\sigma_M = %$(sm)")
            axs[i].title = L"\Omega_R = %$(ΩR)\;\mathrm{eV}"
        end
        i += 1
    end

    fig
end

function max_phot_cont(σx=120, ΩRvals = [0.05, 0.1, 0.2, 0.3], smvals = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5])
    fig = Figure()
    ax = Axis(fig[1,1])

    for ΩR in ΩRvals
        Rstr = "R" * replace(string(ΩR), "."=>"p")
        maxvals = zeros(length(smvals))
        i = 1
        for sm in smvals
            smstr = "sm" * replace(string(sm), '.'=>'p')
            path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/$Rstr/$smstr/out.h5")
            zpm = h5read(path, "$(σx)_avg_zpm")
            ppm = h5read(path, "$(σx)_avg_ppm")
            npm = h5read(path, "$(σx)_avg_npm")
            total = [zpm[i] + sum(ppm[:,i]) + sum(npm[:,i]) for i = eachindex(zpm)]
            maxvals[i] = maximum(total)
            i += 1
        end
        ΩR == 0.1 ? println(maxvals) : nothing
        scatter!(ax, smvals ./ ΩR, maxvals) 
    end
    #lines!(ax, 0:0.1:10, [(0.8*exp(-1.5*α)+0.05) for α = 0:0.1:10], color=:black, linestyle=:dash)

    fig
end