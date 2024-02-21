using Statistics

function SI_fig1(; ΩR=0.1, σx=120, σMvals=[0.04, 0.1, 0.2, 0.3], t=1)
    fontsize_theme = Theme(fontsize=20)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(900, 400))
    gd = fig[1, 1] = GridLayout()
    ax1 = Axis(gd[1, 1], xlabel=L"Distance ($\mu$m)", ylabel=L"|\langle n|\psi(t)\rangle|^2", title="(a) t = $(Int(1000*t)) fs")
    ax2 = Axis(gd[1, 2], yticks=[2, 4, 6, 8], xticks=[0, 1, 2, 3, 4, 5], xlabel=L"\sigma_M / \Omega_R", ylabel=L"RMSD ($\mu$m)", title="(b) t = $(Int(1000*t)) fs")
    ax3 = Axis(gd[1, 2], yticks=[0.2, 0.4, 0.6, 0.8], xticks=[0, 1, 2, 3, 4, 5], ylabel=L"\chi")
    ax3.yaxisposition = :right
    ylims!(ax2, 1, 9)
    ylims!(ax3, 0.1, 0.9)
    hidexdecorations!(ax3)

    for i in eachindex(σMvals)
        σM = σMvals[i]
        get_zoomed_wvp!(ax1, ΩR=ΩR, σx=σx, σM=σM, t=t, normalize=true, color=Makie.wong_colors()[i], label="$(round(10*σM, digits=2))")
    end

    all_σM = [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]
    rmsd = similar(all_σM)
    χ = similar(all_σM)
    for i in eachindex(all_σM)
        rmsd[i], χ[i] = get_rmsd_and_χ(ΩR, σx, all_σM[i], t)
    end

    xlims!(ax1, 25, 35)
    ylims!(ax1, 0.0, 0.03)

    scatter!(ax2, all_σM .* 10, rmsd, color=:firebrick3)
    lines!(ax2, all_σM .* 10, rmsd, color=:firebrick3)
    scatter!(ax3, all_σM .* 10, χ, color=:yellow3, marker=:diamond)
    lines!(ax3, all_σM .* 10, χ, color=:yellow3, marker=:diamond)

    axislegend(ax1, ax1, L"\sigma_M / \Omega_R", merge=true)
    axislegend(ax2, [[MarkerElement(marker=:circle, color=:firebrick3), MarkerElement(marker=:diamond, color=:yellow3)]],
        [["RMSD", "Migration Prob."]], [""])

    fig
end

function SI_fig2(;σx=120, ΩR=0.1, σMvals = [0.02, 0.04, 0.1, 0.2], show_std=false, fit=true, χmax=0.6, rmsdmax=3)
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 800))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1], ylabel=L"Escape Probability ($\chi$)")
    ax2 = Axis(gd[2,1], ylabel=L"RMSD ($\mu$m)", xlabel="Time (ps)",xticks=0:5)
    ax3 = Axis(gd[1,2])
    ax4 = Axis(gd[2,2], xlabel="Time (ps)")

    axs = [ax1, ax2, ax3, ax4]

    ## Link axes so they have the same plotting range
    linkxaxes!(ax1,ax2)
    linkxaxes!(ax3,ax4)

    xlims!(ax2, -0.1, 4)
    ylims!(ax2, 0, 8)
    xlims!(ax4, 0.0, 0.5)
    ylims!(ax4, 0, rmsdmax)
    ylims!(ax3, 0, χmax)


    fig2!(axs, σx=σx, ΩR=ΩR, σMvals=σMvals, show_std=show_std, fit=fit)

    Legend(gd[3,1:2], ax1, L"\sigma_M/\Omega_R", merge=true, orientation=:horizontal, labelsize=20, titleposition=:left)
    text!(axs[3], 0.15, 1, text="(c)", align=(:right, :top), space=:relative, font=:bold, fontsize=25)
    text!(axs[4], 0.15, 1, text="(d)", align=(:right, :top), space=:relative, font=:bold, fontsize=25)
    Label(gd[0,1:2], L"$\Omega_R = %$ΩR$ eV")

    rowgap!(gd, 2, 5)
    fig
end

function get_zoomed_wvp!(ax::Axis; ΩR=0.1, σx=120, σM=0.005, t=0, normalize=false, color=Makie.wong_colors()[1], label="")

    Rstr = replace(string(ΩR), "." => "p")
    sm = replace(string(σM), "." => "p")
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5
    wvp_tvals = vcat(r1, r2)
    wvp_idx = findfirst(x -> x == t, wvp_tvals)

    wvp = h5read(path, "$(Int(σx))_avg_wvp", (:, wvp_idx))

    P = sum(wvp)
    if normalize
        wvp = wvp ./ P
    end

    scatter!(ax, 0:0.5:49.5, wvp, color=color, label=label)
    lines!(ax, 0:0.5:49.5, wvp, color=color, label=label)
end

function get_rmsd_and_χ(ΩR, σx, σM, t)
    Rstr = replace(string(ΩR), "." => "p")
    sm = replace(string(σM), "." => "p")
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5
    r3 = 0.51:0.01:5
    wvp_tvals = vcat(r1, r2)
    d_tvals = vcat(r1, r3)
    wvp_idx = findfirst(x -> x == t, wvp_tvals)
    d_idx = findfirst(x -> x == t, d_tvals)

    d = h5read(path, "$(Int(σx))_avg_d")[d_idx] * 10 / 1000 # Convert to RMSD
    escp = wvp_escape_probability(ΩR, σx, σM, wvp_idx)

    return (d, escp)
end

function plot_dispersion(Emol)

    path = joinpath(@__DIR__, "../../../dispersion_curves/out.h5")

    LP_k = h5read(path, "$(Emol)_LP_k")
    LP = h5read(path, "$(Emol)_LP")
    UP_k = h5read(path, "$(Emol)_UP_k")
    UP = h5read(path, "$(Emol)_UP")


    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$q$ (nm$^{-1}$)", ylabel="Energy (eV)")

    scatter!(ax, LP_k, LP, markersize=3)
    scatter!(ax, UP_k, UP, markersize=3)

    xlims!(ax, -0.02, 0.02)
    ylims!(ax, Emol - 0.1, 2.5)

    fig

end

function plot_groupvelocity()
    fig = Figure()
    gd = fig[1, 1] = GridLayout()
    axs = [Axis(gd[i, j]) for i = 1:3, j = 1:2]

    Evals = [1.5, 1.8, 1.9, 2.0, 2.1, 2.2]

    for i in eachindex(axs)
        δ = round(2.0 - Evals[i], digits=2)
        axs[i].title = L"\delta = %$(δ)"
        plot_groupvelocity!(axs[i], Evals[i])
    end

    linkaxes!(axs...)

    xlims!(axs[end], -0.01, 0.01)
    ylims!(axs[end], -50, 50)

    hidexdecorations!(axs[1], label=true, ticklabels=true, ticks=false, grid=false)
    hidexdecorations!(axs[2], label=true, ticklabels=true, ticks=false, grid=false)
    hidexdecorations!(axs[4], label=true, ticklabels=true, ticks=false, grid=false)
    hidexdecorations!(axs[5], label=true, ticklabels=true, ticks=false, grid=false)

    hideydecorations!(axs[4], label=true, ticklabels=true, ticks=false, grid=false)
    hideydecorations!(axs[5], label=true, ticklabels=true, ticks=false, grid=false)
    hideydecorations!(axs[6], label=true, ticklabels=true, ticks=false, grid=false)

    Label(gd[:, 0], L"\partial E / \partial q", rotation=π / 2)
    Label(gd[4, :], L"$q$ (nm$^{-1}$)")

    colgap!(gd, 3)
    rowgap!(gd, 3)

    fig
end

function plot_groupvelocity(Emol)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$\partial$E/$\partial$q (eV$\cdot$nm)")
    plot_groupvelocity!(ax, Emol)
    fig
end

function plot_groupvelocity!(ax, Emol)
    plot_groupvelocityLP!(ax, Emol, color=Makie.wong_colors()[1])
    plot_groupvelocityUP!(ax, Emol, color=Makie.wong_colors()[2])
end

function plot_groupvelocityLP!(ax, Emol; color=Makie.wong_colors()[1], label="", mksize=5)

    ħ = 0.0006582119569509067 # eV⋅ps
    path = joinpath(@__DIR__, "../../../dispersion_curves/out.h5")

    LP_k = h5read(path, "$(Emol)_LP_k")
    LP = h5read(path, "$(Emol)_LP")

    ∂LP = zeros(length(LP) - 1)
    for i in eachindex(LP_k)

        if i == length(LP)
            break
        end

        δq = LP_k[i+1] - LP_k[i]
        δE = LP[i+1] - LP[i]

        ∂LP[i] = δE / δq
    end

    # Convert Units
    ∂LP = (∂LP ./ ħ) / 1000 # μm ⋅ ps⁻¹
    #LP_k = LP_k * 1e7 # cm⁻¹

    scatter!(ax, LP_k[1:end-1], ∂LP, markersize=mksize, color=color, label=label)
end

function plot_groupvelocityUP!(ax, Emol; color=Makie.wong_colors()[1], label="", mksize=5)

    ħ = 0.0006582119569509067 # eV⋅ps
    path = joinpath(@__DIR__, "../../../dispersion_curves/out.h5")

    UP_k = h5read(path, "$(Emol)_UP_k")
    UP = h5read(path, "$(Emol)_UP")

    ∂UP = zeros(length(UP) - 1)
    for i in eachindex(UP_k)

        if i == length(UP)
            break
        end

        δq = UP_k[i+1] - UP_k[i]
        δE = UP[i+1] - UP[i]

        ∂UP[i] = δE / δq
    end

    # Convert Units
    ∂UP = (∂UP ./ ħ) / 1000 # μm ⋅ ps⁻¹
    scatter!(ax, UP_k[1:end-1], ∂UP, markersize=mksize, color=color, label=label)
end

function plot_groupvelocityLP()

    fontsize_theme = Theme(fontsize=25)
    set_theme!(fontsize_theme)

    fig = Figure()
    gd = fig[1, 1] = GridLayout()
    ax = Axis(gd[1, 1], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$v_g = \hbar^{-1} \partial E / \partial q$ ($\mu$m$\cdot$ps$^{-1}$)")

    Evals = [1.5, 1.8, 1.9, 2.0, 2.1, 2.2]

    for i in eachindex(Evals)
        δ = round(2.0 - Evals[i], digits=2)
        plot_groupvelocityLP!(ax, Evals[i], color=Makie.wong_colors()[i], label="$δ", mksize=5)
    end

    xlims!(ax, 0, 0.02)
    ylims!(ax, 0, 50)

    axislegend(ax, ax, L"$\delta$ (eV)", position=:rt)

    fig
end

function plot_groupvelocityUP()

    fontsize_theme = Theme(fontsize=25)
    set_theme!(fontsize_theme)

    fig = Figure()
    gd = fig[1, 1] = GridLayout()
    ax = Axis(gd[1, 1], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$v_g = \hbar^{-1} \partial E / \partial q$ ($\mu$m$\cdot$ps$^{-1}$)")

    Evals = [1.5, 1.8, 1.9, 2.0, 2.1, 2.2]

    for i in eachindex(Evals)
        δ = round(2.0 - Evals[i], digits=2)
        plot_groupvelocityUP!(ax, Evals[i], color=Makie.wong_colors()[i], label="$δ", mksize=5)
    end

    xlims!(ax, 0, 0.03)
    ylims!(ax, 0, 150)

    axislegend(ax, ax, L"$\delta$ (eV)", position=:rt)

    fig
end

function plot_ω2LP()
    fontsize_theme = Theme(fontsize=25)
    set_theme!(fontsize_theme)

    fig = Figure()
    gd = fig[1, 1] = GridLayout()
    ax = Axis(gd[1, 1], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$\omega_{2LP} = \hbar^{-1} \partial^2 E / \partial q^2$ ($\mu$m$^2\cdot$ps$^{-1}$)")

    Evals = [1.5, 1.8, 1.9, 2.0, 2.1, 2.2]

    for i in eachindex(Evals)
        δ = round(2.0 - Evals[i], digits=2)
        plot_ω2LP!(ax, Evals[i], color=Makie.wong_colors()[i], label="$δ", mksize=5)
    end

    xlims!(ax, 0, 0.02)
    #ylims!(ax, 0, 150)

    axislegend(ax, ax, L"$\delta$ (eV)", position=:rb)

    fig
end
function plot_ω2LP!(ax, Emol; color=Makie.wong_colors()[1], label="", mksize=5)

    ħ = 0.0006582119569509067 # eV⋅ps
    path = joinpath(@__DIR__, "../../../dispersion_curves/out.h5")

    LP_k = h5read(path, "$(Emol)_LP_k")
    LP = h5read(path, "$(Emol)_LP")

    ∂LP = zeros(length(LP) - 1)
    for i in eachindex(LP_k)

        if i == length(LP)
            break
        end

        δq = LP_k[i+1] - LP_k[i]
        δE = LP[i+1] - LP[i]

        ∂LP[i] = δE / δq
    end

    # Convert Units
    ∂LP = (∂LP ./ ħ) / 1000 # μm ⋅ ps⁻¹
    #LP_k = LP_k * 1e7 # cm⁻¹

    ∂2LP = zeros(length(LP) - 2)
    for i in eachindex(LP_k)

        if i == (length(LP) - 1)
            break
        end

        δq = LP_k[i+1] - LP_k[i]
        δ2E = ∂LP[i+1] - ∂LP[i]

        ∂2LP[i] = (δ2E / δq) / 1000 # convert nm -> μm
    end

    scatter!(ax, LP_k[1:end-2], ∂2LP, markersize=mksize, color=color, label=label)
end

function plot_ω2UP()
    fontsize_theme = Theme(fontsize=25)
    set_theme!(fontsize_theme)

    fig = Figure()
    gd = fig[1, 1] = GridLayout()
    ax = Axis(gd[1, 1], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$\omega_{2UP} = \hbar^{-1} \partial^2 E / \partial q^2$ ($\mu$m$^2\cdot$ps$^{-1}$)")

    Evals = [1.5, 1.8, 1.9, 2.0, 2.1, 2.2]

    for i in eachindex(Evals)
        δ = round(2.0 - Evals[i], digits=2)
        plot_ω2UP!(ax, Evals[i], color=Makie.wong_colors()[i], label="$δ", mksize=5)
    end

    xlims!(ax, 0, 0.02)
    #ylims!(ax, 0, 150)

    axislegend(ax, ax, L"$\delta$ (eV)", position=:rb)

    fig
end
function plot_ω2UP!(ax, Emol; color=Makie.wong_colors()[1], label="", mksize=5)

    ħ = 0.0006582119569509067 # eV⋅ps
    path = joinpath(@__DIR__, "../../../dispersion_curves/out.h5")

    UP_k = h5read(path, "$(Emol)_UP_k")
    UP = h5read(path, "$(Emol)_UP")

    ∂UP = zeros(length(UP) - 1)
    for i in eachindex(UP_k)

        if i == length(UP)
            break
        end

        δq = UP_k[i+1] - UP_k[i]
        δE = UP[i+1] - UP[i]

        ∂UP[i] = δE / δq
    end

    # Convert Units
    ∂UP = (∂UP ./ ħ) / 1000 # μm ⋅ ps⁻¹
    #UP_k = UP_k * 1e7 # cm⁻¹

    ∂2UP = zeros(length(UP) - 2)
    for i in eachindex(UP_k)

        if i == (length(UP) - 1)
            break
        end

        δq = UP_k[i+1] - UP_k[i]
        δ2E = ∂UP[i+1] - ∂UP[i]

        ∂2UP[i] = (δ2E / δq) / 1000 # convert nm -> μm
    end

    scatter!(ax, UP_k[1:end-2], ∂2UP, markersize=mksize, color=color, label=label)
end

function plot_term_prob(;σx=120, ΩR=0.1, σMvals = [0.02, 0.04, 0.1, 0.2], N=1)
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(900, 500))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1], ylabel="Probability", xticks=0:4)
    ax2 = Axis(gd[1,2], xticks=0:4)
    xlims!(ax1, 0,5)
    xlims!(ax2, 0,5)

    rel_dis = [format("{:3d}%", Int(round(100*σM/ΩR, digits=0))) for σM in σMvals]

    clrs = Makie.wong_colors()
    for i = eachindex(σMvals) 
        c = clrs[i]
        prob_at_end_sites!(ax1, ΩR=ΩR, σx=σx, σM=σMvals[i], N=1, color=c, label=rel_dis[i])
        prob_at_end_sites!(ax2, ΩR=ΩR, σx=σx, σM=σMvals[i], N=5, color=c, label=rel_dis[i])
    end

    text!(ax1, 0.01, 0.98, text="(a) 100 sites", space=:relative, align=(:left, :top), fontsize=23)
    text!(ax2, 0.01, 0.98, text="(b) 500 sites", space=:relative, align=(:left, :top), fontsize=23)

    Label(gd[2,1:2], "Time (ps)")
    Legend(gd[1,3], ax1, L"\sigma_M/\Omega_R", merge=true)
    rowgap!(gd, 1, 1)
    fig
end

#function fig2!(axs; σx=120, ΩR=0.1, σMvals = [0.02, 0.04, 0.1, 0.2], show_std=false, fit=true, plotlabel=["(a)", "(b)"])
#    rel_dis = [format("{:3d}%", Int(round(100*σM/ΩR, digits=0))) for σM in σMvals]
#
#    clrs = Makie.wong_colors()
#    for i = eachindex(σMvals) 
#        c = clrs[i]
#        escp!(axs[1], ΩR=ΩR, σx=σx, σM=σMvals[i],color=c,label=rel_dis[i])
#        #escp_over_time!(axs[1], ΩR=ΩR, σx=σx, σM=σMvals[i],color=c,label=rel_dis[i])
#        rmsd_propagation!(axs[2], ΩR=ΩR, σx=σx, σM=σMvals[i], color=c, show_std=show_std, fit=fit, label=rel_dis[i])
#        #escp_over_time!(axs[3], ΩR=ΩR, σx=σx, σM=σMvals[i],color=c,label=rel_dis[i], mksize=3)
#        escp!(axs[3], ΩR=ΩR, σx=σx, σM=σMvals[i],color=c,label=rel_dis[i], mksize=3)
#        rmsd_propagation!(axs[4], ΩR=ΩR, σx=σx, σM=σMvals[i], color=c, show_std=show_std, fit=fit)
#    end
#
#    text!(axs[1], 0.15, 1, text=plotlabel[1], align=(:right, :top), space=:relative, font=:bold, fontsize=25)
#    text!(axs[2], 0.15, 1, text=plotlabel[2], align=(:right, :top), space=:relative, font=:bold, fontsize=25)
#end

function SI_fig5(;σx=240)
    fontsize_theme = Theme(fontsize = 23)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 800))
    gd = fig[1,1] = GridLayout()

    ax1 = Axis(gd[1,1], xticks=[-0.2, -0.1, 0.0, 0.1, 0.2], title=L"\bar{q}_0 = 0.0")
    ax2 = Axis(gd[1,2], xticks=[-0.2, -0.1, 0.0, 0.1, 0.2], title=L"$\bar{q}_0 = 0.008$ nm$^{-1}$")
    ax3 = Axis(gd[2,1], xticks=[-0.2, -0.1, 0.0, 0.1, 0.2])
    ax4 = Axis(gd[2,2], xticks=[-0.2, -0.1, 0.0, 0.1, 0.2])
    ax5 = Axis(gd[3,1], xticks=[-0.2, -0.1, 0.0, 0.1, 0.2])
    ax6 = Axis(gd[3,2], xticks=[-0.2, -0.1, 0.0, 0.1, 0.2])
    ax7 = Axis(gd[4,1], xticks=[-0.2, -0.1, 0.0, 0.1, 0.2])
    ax8 = Axis(gd[4,2], xticks=[-0.2, -0.1, 0.0, 0.1, 0.2])
    axs = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    linkaxes!(axs...)
    xlims!(ax8, -0.22, 0.22)

    for i = 1:6
        hidexdecorations!(axs[i], grid=false, ticks=false)
    end

    for i = [2, 4, 6, 8]
        hideydecorations!(axs[i], grid=false, ticks=false)
    end

    for (j,σx) in zip([1, 3, 5, 7], [60, 120, 240, 480])
        for (i,σM) in enumerate([0.005, 0.02, 0.04, 0.1])
            nzq_v0_vs_detuning!(axs[j+1], σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i))
            v0_vs_detuning!(axs[j], σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i))
        end
    end

    Label(gd[5,1:2], L"Detuning $\delta = \hbar\omega_0 - E_M$ (eV)")
    Label(gd[1:4,0], L"$v_0$ ($\mu$m$\cdot$ps$^{-1}$)", rotation=π/2)
    Legend(gd[6,1:2], ax1, L"\sigma_M/\Omega_R", orientation=:horizontal, merge=true, titleposition=:left)

    rowgap!(gd, 1, Relative(0.008))
    rowgap!(gd, 2, Relative(0.008))
    rowgap!(gd, 3, Relative(0.008))
    rowgap!(gd, 4, Relative(0.008))
    rowgap!(gd, 5, Relative(0.008))
    colgap!(gd, 1, Relative(0.008))

    for i in 1:8
        l = '`' + i
        σx = [60, 60, 120, 120, 240, 240, 480, 480][i]
        text!(axs[i], -0.1, 20, align=(:left, :bottom), text=L"(%$l) $\sigma_x = %$σx$ nm")
    end

    fig
end

function SI_extrafig5(;σx=240)
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 400))
    gd = fig[1,1] = GridLayout()

    ax1 = Axis(gd[1,1], ylabel=L"$v_0$ ($\mu$m$\cdot$ps$^{-1}$)")
    ax2 = Axis(gd[1,2], ylabel=L"Maximum RMSD ($\mu$m)")
    linkxaxes!(ax1, ax2)

    for (i,σM) in enumerate([0.01, 0.02, 0.04, 0.1])
    # Plot detuning effect on v0
        v0_vs_detuning!(ax1, σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i), Emvals=[1.8, 1.9, 2.0, 2.1, 2.2, 2.5, 3.0, 4.0, 5.0])
        maxd_vs_detuning!(ax2, σx=σx, σM=σM,label="$(10*σM)", color=Makie.Cycled(i), Emvals=[1.8, 1.9, 2.0, 2.1, 2.2, 2.5, 3.0, 4.0, 5.0])
    end

    Label(gd[2,1:2], L"Detuning $\delta = \hbar\omega_0 - E_M$ (eV)")
    Legend(gd[3,1:2], ax1, L"\sigma_M/\Omega_R", orientation=:horizontal, merge=true, titleposition=:left)

    text!(ax1, 0.01, 0.98, text="(a)", space=:relative, align=(:left, :top), font=:bold, fontsize=22)
    text!(ax2, 0.01, 0.98, text="(b)", space=:relative, align=(:left, :top), font=:bold, fontsize=22)

    rowgap!(gd, 1, Relative(0.008))
    rowgap!(gd, 2, Relative(0.008))

    fig
end

function SI_rabi_oscillations(;σx=120)


    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5
    tvals = vcat(r1, r2)
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 500))
    gd = fig[1,1] = GridLayout()

    ax1 = Axis(gd[1,1], xticks=[0, 50, 100, 150, 200, 250], title="(a)", titlealign=:left)
    ax2 = Axis(gd[2,1], xticks=[0, 20, 40, 60, 80], title="(b)", titlealign=:left)
    xlims!(ax1, 0, 260)
    xlims!(ax2, 0, 60)

    mk = vcat(repeat([:circle], 4), repeat([:diamond], 4))
    axs = vcat(repeat([ax1], 4), repeat([ax2], 4))
    cl = vcat(Makie.wong_colors()[1:7], Makie.wong_colors()[1])

    for (i,σM) in enumerate([0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3])
        sm = replace(string(σM), "." => "p") 
        path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R0p1/sm$sm/out.h5")
        wvp = h5read(path, "$(Int(σx))_avg_wvp")
        PM = [sum(wvp[:,t]) for t in axes(wvp, 2)]        
        scatter!(axs[i], 1000*tvals, PM, label="$(10*σM)", marker=mk[i], color=cl[i])
        lines!(axs[i], 1000*tvals, PM, label="$(10*σM)", color=cl[i])
    end

    Legend(gd[1,2], ax1, L"\sigma_M/\Omega_R", orientation=:vertical, merge=true, titleposition=:top, nbanks=2)
    Legend(gd[2,2], ax2, L"\sigma_M/\Omega_R", orientation=:vertical, merge=true, titleposition=:top, nbanks=2)
    Label(gd[1:2, 0], L"Excitonic Probability ($P_M$)", rotation=π/2, tellheight=false)
    Label(gd[3, 1], "Time (fs)", tellwidth=false)

    rowgap!(gd, 1, Relative(0.015))
    rowgap!(gd, 2, Relative(0.015))

    fig
end

function SI_dispersion()

    path = joinpath(@__DIR__, "../../../dispersion_curves/disorder/out.h5")
    fontsize_theme = Theme(fontsize = 22)
    set_theme!(fontsize_theme)
    fig = Figure(resolution=(800, 600))

    gd = fig[1,1] = GridLayout()

    axs = [
        Axis(gd[i,j]) for i = 1:3, j = 1:2
    ]

    σMvals = [0.01, 0.02, 0.04, 0.08, 0.1, 0.2]

    for i in eachindex(σMvals)

        σM = σMvals[i]
        sm = replace(string(σM), "." => "p")

        lt = '`' + i
        axs[i].title = L"(%$lt) $\sigma_M / \Omega_R = %$(10*σM)"

        xlims!(axs[i], -0.005, 0.005)
        ylims!(axs[i], 1.9, 2.1)

        qup = h5read(path, "$(sm)_up_q")
        up = h5read(path, "$(sm)_up")
        qlp = h5read(path, "$(sm)_lp_q")
        lp = h5read(path, "$(sm)_lp")

        x1, y1, σ1 = remove_duplicates_and_average(qup, up)
        x2, y2, σ2 = remove_duplicates_and_average(qlp, lp)

        if i != 3 && i != 6
            hidexdecorations!(axs[i], grid=false, ticks=false)
        end

        if i > 3
            hideydecorations!(axs[i], grid=false, ticks=false)
        end

        #scatter!(axs[i], qup, up, markersize=9)
        #scatter!(axs[i], qlp, lp, markersize=9)
        scatter!(axs[i], x1, y1)
        errorbars!(axs[i], x1, y1, 3*σ1, color=Makie.wong_colors()[1])

        scatter!(axs[i], x2, y2)
        errorbars!(axs[i], x2, y2, 3*σ2, color=Makie.wong_colors()[2])
        hlines!(axs[i], [2.0], linestyle=:dash, color=:gray)

        phot_q = h5read(path, "phot_q")
        phot_e = h5read(path, "phot_e")
        s = sortperm(phot_q)
        lines!(axs[i], phot_q[s], phot_e[s], color=:gray, linestyle=:dash)
    end

    Label(gd[1:3, 0], "Energy (eV)", rotation=π/2)
    Label(gd[4,1:2], L"q (nm$^{-1}$)")

    colgap!(gd, 1, Relative(0.015))
    colgap!(gd, 2, Relative(0.015))

    rowgap!(gd, 1, Relative(0.015))
    rowgap!(gd, 2, Relative(0.015))
    rowgap!(gd, 3, Relative(0.015))


    fig
end

function remove_duplicates_and_average(x, y)
    unique_x = unique(x)
    averaged_y = [mean(y[x .== xi]) for xi in unique_x]
    std_y = [std(y[x .== xi]) for xi in unique_x]
    return unique_x, averaged_y, std_y
end