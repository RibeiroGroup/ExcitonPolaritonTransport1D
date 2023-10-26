using PolaritonicSystems
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
    dEdq = numerical_derivative_three_point(XP.qvals, XP.evals)
    
    return (XP.mol_cont .* dEdq) ./ (1000*ħ)
end

function plot_excitonic_vg(ΩR, Em)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$v_g^\text{exc.}$ ($\mu$m$\cdot$ps$^{-1}$)")
    plot_excitonic_vg!(ax, ΩR, Em)
    axislegend(ax, merge=true)
    fig
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

function plot_excitonic_vg!(ax, ΩR, Em)
    sys = QuantumWire(
           ΩR = ΩR*eV,
           ϵ  = 3.0,
           Nc = 500,
           Nm = 1001,
           a  = 10nm,
           σa = 0,
           ωM = Em*eV,
           σM = 0.0,
           Ly = 200nm,
           Lz = 400nm,
           nz = 1,
           ny = 1
       )
    
    LP = extract_LP(sys)
    UP = extract_UP(sys)

    ∂LP = get_excitonic_vg(LP)
    ∂UP = get_excitonic_vg(UP)

    #scatter!(ax, LP.qvals, ∂LP, label="LP")
    lines!(ax, LP.qvals, ∂LP, label="LP", linewidth=3)

    #scatter!(ax, UP.qvals, ∂UP, label="UP")
    lines!(ax, UP.qvals, ∂UP, label="UP", linewidth=3)
end

function plot_excitonic_wvp!(ax, ΩR, Em, wvps)
    sys = QuantumWire(
           ΩR = ΩR*eV,
           ϵ  = 3.0,
           Nc = 500,
           Nm = 1001,
           a  = 10nm,
           σa = 0,
           ωM = Em*eV,
           σM = 0.0,
           Ly = 200nm,
           Lz = 400nm,
           nz = 1,
           ny = 1
       )
    
    LP = extract_LP(sys)
    UP = extract_UP(sys)

    # Get k-space wave packet
    for w in wvps
        σx = w[1]
        q0 = w[2]
        h = w[3]
        kwvp = get_kwvp(sys, σx, q0)
        kwvp = kwvp ./ maximum(kwvp)
        cmap = [(:steelblue2, t) for t in kwvp]
        lines!(ax, LP.qvals, repeat([h], length(LP.qvals)), color=cmap, linewidth=12)
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

# Paper figure
function plot_excvg()

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    fig = Figure(resolution=(800, 400))

    # Create grid of axes
    gd = fig[1,1] = GridLayout()

    ax1 = Axis(gd[1,1], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$v_g^\text{exc.}$ ($\mu$m$\cdot$ps$^{-1}$)")
    ax2 = Axis(gd[1,2], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$v_g^\text{exc.}$ ($\mu$m$\cdot$ps$^{-1}$)")
    ax3 = Axis(gd[1,3], xlabel=L"$q$ (nm$^{-1}$)", ylabel=L"$v_g^\text{exc.}$ ($\mu$m$\cdot$ps$^{-1}$)")

    ylims!(ax1, 0, 12)
    xlims!(ax1, 0, 0.015)

    ylims!(ax2, 0, 11)
    xlims!(ax2, 0, 0.015)

    ylims!(ax3, 0, 20)
    xlims!(ax3, 0, 0.015)

    # Rabi analysis
    plot_excitonic_wvp!(ax1, 0.1, 2.0, [(120, 0.0, 8), (480, 0.0, 4)])
    plot_excitonic_vg!(ax1, 0.2, 2.0)
    plot_excitonic_vg!(ax1, 0.1, 2.0)

    # q diff 0 analysis
    plot_excitonic_wvp!(ax2, 0.1, 2.0, [(60, 0.0, 10), (120, 0.0, 8), (240, 0.0, 6), (360, 0.0, 4), (480, 0.0, 2),
    (60, 0.0055, 9), (120, 0.0055, 7), (240, 0.0055, 5), (360, 0.0055, 3), (480, 0.0055, 1)])
    plot_excitonic_vg!(ax2, 0.1, 2.0)

    # detuning
    plot_excitonic_vg!(ax3, 0.1, 1.8)
    plot_excitonic_vg!(ax3, 0.1, 2.0)
    plot_excitonic_vg!(ax3, 0.1, 2.2)

    fig
end