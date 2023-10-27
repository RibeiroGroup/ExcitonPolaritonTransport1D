using LinearAlgebra

struct PolaritonBranch{T}
    evals::Array{T}
    qvals::Array{T}
    mol_cont::Array{T}
end

# This function returns the photonic content of the eigenstate i
function phot_cont(sys, i)
    return sum(abs2.(sys.Uix[sys.phot_range, i]))
end

# This function returns the molecular content of the eigenstate i
function mol_cont(sys, i)
    return 1 - phot_cont(sys, i)
end

function extract_LP(sys)
    N = length(sys.evals) ÷ 2
    
    LP = sys.evals[1:N]
    Π = zeros(N)
    for i in eachindex(Π)
        Π[i] = mol_cont(sys, i)
    end

    m = sortperm(sys.phot_wavevectors)
    
    return PolaritonBranch(LP[m], sys.phot_wavevectors[m], Π[m])
end

function extract_UP(sys)
    N = length(sys.evals) ÷ 2
    
    UP = sys.evals[(N+1):end]
    Π = zeros(N)
    for i in eachindex(Π)
        Π[i] = mol_cont(sys, i+N)
    end
    
    m = sortperm(sys.phot_wavevectors)
    
    return PolaritonBranch(UP[m], sys.phot_wavevectors[m], Π[m])
end

function get_excitonic_vg(XP)
    ħ = 0.0006582119569509067 # eV⋅ps
    
    #∂XP = zeros(length(XP.evals))
    #for i in eachindex(XP.qvals)

    #    if i == length(XP.evals)
    #        ∂XP[i] = ∂XP[i-1]
    #        break
    #    end

    #    δq = XP.qvals[i+1] - XP.qvals[i] 
    #    δE = XP.evals[i+1] - XP.evals[i]

    #    ∂XP[i] = XP.mol_cont[i] .* (δE / δq) / (1000*ħ) # Convert to μm/ps
    #end
    dEdq = numerical_derivative_three_point(XP.qvals, XP.evals) ./ (1000*ħ)
    
    return XP.mol_cont .* dEdq
end

function numerical_derivative_three_point(x, y)
    n = length(x)
    dy_dx = zeros(n)  # Initialize the derivative array

    for i in 2:n-1
        dy_dx[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    end

    # Calculate the endpoints using a forward and a backward difference
    dy_dx[1] = (y[2] - y[1]) / (x[2] - x[1])
    dy_dx[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])

    return dy_dx
end

function plot_excitonic_vg!(ax; sys, c)
    LP = extract_LP(sys)
    UP = extract_UP(sys)

    ∂LP = get_excitonic_vg(LP)
    ∂UP = get_excitonic_vg(UP)

    #scatter!(ax, LP.qvals, ∂LP, label="LP")
    lines!(ax, LP.qvals, ∂LP, label="LP", linewidth=3, color=c)

    #scatter!(ax, UP.qvals, ∂UP, label="UP")
    lines!(ax, UP.qvals, ∂UP, label="UP", linewidth=3, linestyle=:dash, color=c)
end

function plot_excitonic_wvp!(ax; sys, wvps, ymax)
    LP = extract_LP(sys)
    UP = extract_UP(sys)

    # Get k-space wave packet
    for w in wvps
        σx = w[1]
        q0 = w[2]
        h = w[3]
        kwvp = get_kwvp(sys, σx, q0)
        kwvp = kwvp ./ maximum(kwvp)
        if q0 > 0 
            cmap = [(:salmon3, t) for t in kwvp]
        else
            cmap = [(:steelblue2, t) for t in kwvp]
        end
        lines!(ax, LP.qvals, repeat([ymax*h], length(LP.qvals)), color=cmap, linewidth=12)
        text!(ax, 0.01, h+0.02, align=(:left, :bottom), space=:relative, text=L"$\sigma_x = %$σx$ nm", fontsize=15)
    end
end

function get_kwvp(sys, σx, q0)
    Um = get_ktrans(sys)
    wvp = create_exciton_wavepacket(5005, σx, sys, q=q0)
    kwvp = Um * sys.Uix * wvp
    kwvp = abs2.(kwvp[sys.mol_range])
    m = sortperm(sys.phot_wavevectors)
    return kwvp[m]
end

function get_ktrans(sys)
    Nm = length(sys.mol_range)
    Nc = length(sys.phot_range)
    Um = zeros(ComplexF64, Nc, Nm)
    
    for i = 1:Nc
        q = sys.phot_wavevectors[i]
        for j = 1:Nm
            r = sys.mol_positions[j]
            Um[i,j] = exp(-im*q*r)
        end
    end
    
    U = zeros(ComplexF64, Nc+Nm, Nc+Nm)
    
    U[1:Nc, 1:Nc] .= I(Nc)
    
    U[(Nc+1):end, (Nc+1):end] .= Um ./ √(Nm)
    
    return U
end

function effective_wvp_velocity(sys, σx, q0)
    kwvp = get_kwvp(sys, σx, q0)

    LP = extract_LP(sys)
    UP = extract_UP(sys)

    ∂LP = get_excitonic_vg(LP) .^2
    ∂UP = get_excitonic_vg(UP) .^2

    return sqrt(sum(kwvp .* (∂LP .+ ∂UP)))
end

# Paper figure
function plot_excvg(sysA, sysB, sysC, sysD)

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 800))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()

    ax1 = Axis(gd[1,1], xticks=[0, 0.005, 0.010])
    hidexdecorations!(ax1, grid=false, ticks=false)

    ax2 = Axis(gd[1,2], xticks=[0, 0.005, 0.010])
    hidexdecorations!(ax2, grid=false, ticks=false)
    hideydecorations!(ax2, grid=false, ticks=false)

    ax3 = Axis(gd[2,1], xticks=[0, 0.005, 0.010])
    ax4 = Axis(gd[2,2], xticks=[0, 0.005, 0.010])
    hideydecorations!(ax4, grid=false, ticks=false)

    linkxaxes!(ax4, ax3, ax2, ax1)
    linkyaxes!(ax4, ax3)
    linkyaxes!(ax2, ax1)

    ylims!(ax2, 0, 12)
    ylims!(ax4, 0, 20)
    xlims!(ax4, -0.002, 0.015)

    wvps0 = [(60, 0.0, 5/6), (120, 0.0, 4/6), (240, 0.0, 3/6), (360, 0.0, 2/6), (480, 0.0, 1/6)]
    wvps = [(60, 0.0055, 5/6), (120, 0.0055, 4/6), (240, 0.0055, 3/6), (360, 0.0055, 2/6), (480, 0.0055, 1/6)]

    c1 = Makie.wong_colors()[3]
    c2 = Makie.wong_colors()[2]
    c3 = Makie.wong_colors()[1]
    c4 = Makie.wong_colors()[4]

    # Rabi analysis
    plot_excitonic_wvp!(ax1, sys=sysA, wvps=wvps0, ymax=12)
    plot_excitonic_vg!(ax1; sys=sysA, c=c1)
    plot_excitonic_vg!(ax1, sys=sysB, c=c2)

    plot_excitonic_wvp!(ax2, sys=sysA, wvps=wvps, ymax=12)
    plot_excitonic_vg!(ax2; sys=sysA, c=c1)
    plot_excitonic_vg!(ax2, sys=sysB, c=c2)

    # detuning
    plot_excitonic_wvp!(ax3, sys=sysA, wvps=wvps0, ymax=20)
    plot_excitonic_vg!(ax3; sys=sysC, c=c3)
    plot_excitonic_vg!(ax3, sys=sysA, c=c1)
    plot_excitonic_vg!(ax3, sys=sysD, c=c4)

    plot_excitonic_wvp!(ax4, sys=sysA, wvps=wvps, ymax=20)
    plot_excitonic_vg!(ax4; sys=sysC, c=c3)
    plot_excitonic_vg!(ax4, sys=sysA, c=c1)
    plot_excitonic_vg!(ax4, sys=sysD, c=c4)

    Label(gd[3,1:2], L"$q$ (nm$^{-1}$)")
    Label(gd[1:2,0], L"$v_g^\text{exc.}$ ($\mu$m$\cdot$ps$^{-1}$)", rotation=π/2)


    n = 0
    for ax in [ax1, ax2, ax3, ax4]
        text!(ax, 0.02, 1, align=(:left, :top), space=:relative, text="($('a'+n))", font=:bold, fontsize=20)
        n += 1
    end

    group_line = [LineElement(color=:black, linestyle=:solid), LineElement(color=:black, linestyle=:dash)]
    group_colors = [PolyElement(color=c, strokecolor=:transparent) for c in [c1, c3, c4, c2]]
    Legend(gd[4,1:2], [group_line, group_colors], [["LP", "UP"], 
    [L"$\Omega_R = 0.1$ eV $\delta = 0.0$ eV", L"$\Omega_R = 0.1$ eV $\delta = 0.2$ eV", 
    L"$\Omega_R = 0.1$ eV $\delta = -0.2$ eV", L"$\Omega_R = 0.2$ eV $\delta = 0.0$ eV" ]], 
    ["", ""], orientation=:horizontal, nbanks=2, titleposition=:left)

    rowgap!(gd, 2, 5)
    rowgap!(gd, 3, 5)
    colgap!(gd, 1, 5)

    fig
end