using Makie
using HDF5
using LaTeXStrings
using Formatting

function fig2(;σx=120, ΩR=0.1, σMvals = [0.02, 0.04, 0.1, 0.2], show_std=false, fit=true)
    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(550, 800))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    ax1 = Axis(gd[1,1], ylabel=L"Escape Probability ($\chi$)")
    ax2 = Axis(gd[2,1], ylabel=L"RMSD ($\mu$m)", xlabel="Time (ps)",xticks=0:5)

    # The inset axis
    ax3 = Axis(gd[1, 1],
    width=Relative(0.4),
    height=Relative(0.4),
    halign=0.9,
    valign=0.2,
    backgroundcolor=:lightgray,
    xticklabelsize=15, xticks=[0.1, 0.3, 0.5],
    yticklabelsize=15, yticks=[0.1, 0.3, 0.5])

    ax4 = Axis(gd[2, 1],
    width=Relative(0.4),
    height=Relative(0.4),
    halign=0.9,
    valign=0.2,
    backgroundcolor=:lightgray,
    xticklabelsize=15, xticks=[0.1, 0.3, 0.5],
    yticklabelsize=15, yticks=[0, 1, 2])

    axs = [ax1, ax2, ax3, ax4]
    # Link axes so they have the same plotting range
    linkxaxes!(ax1,ax2)
    linkxaxes!(ax3,ax4)

    xlims!(ax2, -0.1, 4)
    ylims!(ax2, 0, 5)
    xlims!(ax4, 0.0, 0.5)
    ylims!(ax4, 0, 3)

    # Hide redundant axis info

    #rel_dis = ["$(Int(100*σM/ΩR))%" for σM in σMvals]
    #rel_dis = ["$(round(100*σM/ΩR, digits=2))%" for σM in σMvals]
    rel_dis = [format("{:3d}%", Int(100*σM/ΩR)) for σM in σMvals]

    clrs = Makie.wong_colors()
    for i = eachindex(σMvals) 
        c = clrs[i]
        escp_over_time!(axs[1], ΩR=ΩR, σx=σx, σM=σMvals[i],color=c,label=rel_dis[i])
        rmsd_propagation!(axs[2], ΩR=ΩR, σx=σx, σM=σMvals[i], color=c, show_std=show_std, fit=fit)
        escp_over_time!(axs[3], ΩR=ΩR, σx=σx, σM=σMvals[i],color=c,label=rel_dis[i], mksize=3)
        rmsd_propagation!(axs[4], ΩR=ΩR, σx=σx, σM=σMvals[i], color=c, show_std=show_std, fit=fit)
    end

    text!(ax1, 0.1, 1, text="(a)", align=(:right, :top), space=:relative, font=:bold, fontsize=25)
    text!(ax2, 0.1, 1, text="(b)", align=(:right, :top), space=:relative, font=:bold, fontsize=25)
    #linkxaxes!(ax1, ax2)

  
    Legend(gd[3,1], ax1, L"\sigma_M/\Omega_R", merge=true, orientation=:horizontal, labelsize=20, titleposition=:left)

    rowgap!(gd, 2, 5)
    fig
end


function wvp_fig2(;σxvals=[120, 480], ΩR=0.1, σMvals = [0.02, 0.04, 0.1, 0.2], show_std=false)

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 800))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()
    axs = [Axis(gd[i,j]) for i = 1:length(σMvals), j = 1:2]
    linkxaxes!(axs...)
    for i = eachindex(σMvals)
        linkyaxes!(axs[i,2], axs[i,1])
    end
    xlims!(axs[end,2], 15, 35)

    axs[1,1].title = L"$\sigma_x$ = %$(σxvals[1]) nm"
    axs[1,2].title = L"$\sigma_x$ = %$(σxvals[2]) nm"

    rel_dis = ["$(round(σM/ΩR, digits=2))" for σM in σMvals]

    clrs = Makie.wong_colors()
    for j in 1:2
        σx = σxvals[j]
        for i = eachindex(σMvals) 
            c = clrs[i]
            σM = σMvals[i]
            get_wvp!(axs[i,j], ΩR=ΩR, σx=σx, σM=σM, t=5, normalize=true, color=c)
            text!(axs[i,j], 0.9, 0.9, text=L"\sigma_M / \Omega_R = %$(rel_dis[i])",align=(:right,:top), fontsize=20, space=:relative)
        end
    end

    Label(gd[length(σMvals)+1,1:2], "Distance μm")

    fig
end