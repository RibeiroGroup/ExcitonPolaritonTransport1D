function fig4()
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(600, 800))
    gd = fig[1,1] = GridLayout()

    ax1 = Axis(gd[1,1], xticks=[100, 200, 300, 400])
    ax2 = Axis(gd[2,1], xlabel=L"Wave packt initial spread - $\sigma_x$ (nm)", xticks=[100, 200, 300, 400])
    Label(gd[1:2,0], L"Initial spread velocity - $v_0$ ($\mu$m$\cdot$ps$^{-1}$)", rotation=π/2)

    hidexdecorations!(ax1, grid=false, ticks=false)

    linkaxes!(ax1, ax2)
    xlims!(ax2, 50,500)

    σMvals = [0.02, 0.04, 0.1]
    σxvals = [60, 120, 240, 360, 480]
    r1 = 0:0.005:0.5

    # q = 0
    mkers = [:circle, :diamond, :xcross]
    for (i,σM) in enumerate(σMvals)
        v0vals = zeros(5)
        for k in eachindex(σxvals)
            sm = replace(string(σM), "." => "p") 
            path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R0p1/sm$sm/out.h5")
            d = 0.01 * h5read(path, "$(σxvals[k])_avg_d")
            a,b = get_linear_fit(r1, d[1:length(r1)])
            v0vals[k] = b
        end
        scatter!(ax1, σxvals, v0vals, label="$(Int(1000*σM))%", marker=mkers[i])
        lines!(ax1, σxvals, v0vals, label="$(Int(1000*σM))%")
    end

    # q ≠ 0
    for (i,σM) in enumerate(σMvals)
        v0vals = zeros(5)
        for k in eachindex(σxvals)
            sm = replace(string(σM), "." => "p") 
            path = joinpath(@__DIR__, "../../../propagation_study/non_zero_q/Nm5000_Nc500_a10_Em2p0/R0p1/sm$sm/out.h5")
            d = 0.01 * h5read(path, "$(σxvals[k])_avg_d")
            a,b = get_linear_fit(r1, d[1:length(r1)])
            v0vals[k] = b
        end
        scatter!(ax2, σxvals, v0vals, label="$(Int(1000*σM))%", marker=mkers[i])
        lines!(ax2, σxvals, v0vals, label="$(Int(1000*σM))%")
    end

    text!(ax1, 0.1, 1, text="(a)", align=(:right, :top), space=:relative, font=:bold, fontsize=25)
    text!(ax1, 0.25, 1, text=L"\bar{q}_0 = 0", align=(:right, :top), space=:relative, font=:bold, fontsize=25)
    text!(ax2, 0.1, 1, text="(b)", align=(:right, :top), space=:relative, font=:bold, fontsize=25)
    text!(ax2, 0.25, 1, text=L"\bar{q}_0 \neq 0", align=(:right, :top), space=:relative, font=:bold, fontsize=25)
    Legend(gd[3,1], ax2, L"\sigma_M/\Omega_R", orientation=:horizontal, merge=true, titleposition=:left)

    fig
end