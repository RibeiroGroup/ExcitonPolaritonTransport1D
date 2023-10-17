function SI_fig1(; ΩR=0.1, σx=120, σMvals=[0.04, 0.08, 0.1, 0.2], tvals=[0.1, 0.5, 1])
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(900, 400))
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1], xlabel=L"Distance ($\mu$m)", ylabel=L"|\langle n|\psi(t)\rangle|^2")
    ax2 = Axis(gd[1,2], yticks=[2, 4, 6, 8], xticks=[0, 1, 2, 3, 4, 5], xlabel=L"\sigma_M / \Omega_R", ylabel=L"RMSD ($\mu$m)")
    ax3 = Axis(gd[1,2], yticks=[0.2, 0.4, 0.6, 0.8], xticks=[0, 1, 2, 3, 4, 5], ylabel=L"\chi")
    ax3.yaxisposition = :right
    ylims!(ax2, 1,9)
    ylims!(ax3, 0.1,0.9)
    hidexdecorations!(ax3)

    t = 1

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
    scatter!(ax3, all_σM .* 10, χ, color=:yellow3, marker=:diamond)

    axislegend(ax1, ax1, L"\sigma_M / \Omega_R", merge=true)
    axislegend(ax2, [[MarkerElement(marker=:circle, color=:firebrick3), MarkerElement(marker=:diamond, color=:yellow3)]],
    [["RMSD", "Escape Prob."]], [""])

    fig
end

function get_zoomed_wvp!(ax::Axis; ΩR=0.1, σx=120, σM=0.005, t=0, normalize=false, color=Makie.wong_colors()[1], label="")

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5
    wvp_tvals = vcat(r1, r2)
    wvp_idx = findfirst(x->x==t, wvp_tvals)

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
    d_tvals = vcat(r1,r3)
    wvp_idx = findfirst(x->x==t, wvp_tvals)
    d_idx = findfirst(x->x==t, d_tvals)

    d = h5read(path, "$(Int(σx))_avg_d")[d_idx] * 10 / 1000 # Convert to RMSD
    escp = wvp_escape_probability(ΩR, σx, σM, wvp_idx)

    return (d,escp)
end

function plot_dispersion(Emol)

    path = joinpath(@__DIR__, "../../../dispersion_curves/out.h5")

    LP_k = h5read(path, "$(Emol)_LP_k")
    LP =   h5read(path, "$(Emol)_LP")
    UP_k = h5read(path, "$(Emol)_UP_k")
    UP =   h5read(path, "$(Emol)_UP")


    fig = Figure()
    ax = Axis(fig[1,1], xlabel=L"$q$ (nm$^{-1}$)", ylabel="Energy (eV)")

    scatter!(ax, LP_k, LP, markersize=3)
    scatter!(ax, UP_k, UP, markersize=3)

    xlims!(ax, -0.02, 0.02)
    ylims!(ax, Emol-0.1, 2.5)

    fig

end

function plot_groupvelocity()
    fig = Figure()
    gd = fig[1,1] = GridLayout()
    axs = [Axis(gd[i,j]) for i = 1:3, j = 1:2]

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

    Label(gd[:,0], L"\partial E / \partial q", rotation=π/2)
    Label(gd[4,:], L"$q$ (nm$^{-1}$)")

    colgap!(gd, 3)
    rowgap!(gd, 3)

    fig
end

function plot_groupvelocity(Emol)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$\partial$E/$\partial$q (eV$\cdot$nm)")
    plot_groupvelocity!(ax, Emol)
    fig
end

function plot_groupvelocity!(ax, Emol)

    path = joinpath(@__DIR__, "../../../dispersion_curves/out.h5")

    LP_k = h5read(path, "$(Emol)_LP_k")
    LP =   h5read(path, "$(Emol)_LP")
    UP_k = h5read(path, "$(Emol)_UP_k")
    UP =   h5read(path, "$(Emol)_UP")

    ∂UP = zeros(length(UP) - 1)
    for i in eachindex(UP_k)

        if i == length(UP)
            break
        end

        δq = UP_k[i+1] - UP_k[i] 
        δE = UP[i+1] - UP[i]

        ∂UP[i] = δE / δq
    end

    ∂LP = zeros(length(LP) - 1)
    for i in eachindex(LP_k)

        if i == length(LP)
            break
        end

        δq = LP_k[i+1] - LP_k[i] 
        δE = LP[i+1] - LP[i]

        ∂LP[i] = δE / δq
    end

    scatter!(ax, LP_k[1:end-1], ∂LP, markersize=3)
    scatter!(ax, UP_k[1:end-1], ∂UP, markersize=3)
end