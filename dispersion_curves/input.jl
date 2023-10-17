using PolaritonicSystems
using MKL
using LinearAlgebra
using Unitful
using HDF5

output("Starting dispersion computations")

# Output file path
path = joinpath(@__DIR__, "out.h5")

Evals = [1.5, 1.8, 1.9, 2.0, 2.1, 2.2]

function extract_branches(sys, Emol)
    LP_k = Float64[]
    LP = Float64[]
    UP = Float64[]
    UP_k = Float64[]

    for i in eachindex(sys.evals)
        c,idx = findmax(abs2.(sys.Uix[sys.phot_range,i]))
        if c < 1e-3
            continue
        end

        if sys.evals[i] < Emol
            push!(LP_k, sys.phot_wavevectors[idx])
            push!(LP, sys.evals[i])
        elseif sys.evals[i] > Emol
            push!(UP_k, sys.phot_wavevectors[idx])
            push!(UP, sys.evals[i])
        end
    end

    pLP = sortperm(LP_k)
    pUP = sortperm(UP_k)

    return LP_k[pLP], LP[pLP], UP_k[pUP], UP[pUP]
end

for E in Evals
    sys = QuantumWire(
        ΩR = 0.1eV,
        ϵ  = 3.0,
        Nc = 200,
        Nm = 5000,
        a  = 10nm,
        σa = 0,
        ωM = E*eV,
        σM = 0,
        Ly = 200nm,
        Lz = 400nm,
        nz = 1,
        ny = 1,
        printout = true
    )
    
    LP_k, LP, UP_k, UP = extract_branches(sys, E)

    h5write(path, "$(E)_LP_k", LP_k)
    h5write(path, "$(E)_LP", LP)
    h5write(path, "$(E)_UP_k", UP_k)
    h5write(path, "$(E)_UP", UP)
end

