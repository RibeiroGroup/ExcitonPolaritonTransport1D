using PolaritonicSystems
using MKL
using LinearAlgebra
using Unitful
using HDF5

output("Starting dispersion computations")

# Output file path
path = joinpath(@__DIR__, "out.h5")
NR = 25

function Pphot(sys, i)
    return sum(abs2.(sys.Uix[sys.phot_range, i]))
end

function max_q(sys, i) 
    vphot = abs2.(sys.Uix[sys.phot_range, i])

    _, m = findmax(vphot)

    return sys.phot_wavevectors[m]
end

for σM in [0.005, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3]


    sm = replace(string(σM), "." => "p")
    output("σx = $σM")
    up   = Float64[]
    up_q = Float64[]
    lp   = Float64[]
    lp_q = Float64[]
    dark = Float64[]

    for N in 1:NR
        output("N = $N")
        sys = QuantumWire(
            ΩR = 0.1eV,
            ϵ  = 3.0,
            Nc = 200,
            Nm = 1000,
            a  = 10nm,
            σa = 10,
            ωM = 2.0eV,
            σM = σM,
            Ly = 200nm,
            Lz = 400nm,
            nz = 1,
            ny = 1,
        )

        for i in eachindex(sys.evals)
            P = Pphot(sys, i)

            if P < 0.1
                push!(dark, sys.evals[i])
            elseif sys.evals[i] < 2.0
                push!(lp_q, max_q(sys, i))
                push!(lp, sys.evals[i])
            else
                push!(up_q, max_q(sys, i))
                push!(up, sys.evals[i])
            end
        end
    end

    h5write(path, "$(sm)_dark", dark)
    h5write(path, "$(sm)_up", up)
    h5write(path, "$(sm)_up_q", up_q)
    h5write(path, "$(sm)_lp", lp)
    h5write(path, "$(sm)_lp_q", lp_q)
end
