using PolaritonicSystems
using MKL
using LinearAlgebra
using Unitful
using HDF5
using Statistics
using UnicodePlots
using Dates

output("Initilizing Simulation.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))

path = joinpath(@__DIR__, "out.h5")

# Number of cavity modes
NC = 200

# Wavepacket initial spread values
σxvals = [120]

# Time range
tvals = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1, 1.2, 1.5, 2]

locstates = zeros(5401, length(tvals))

### INPUTS
# System inputs
sys = QuantumWire(
    ΩR = 0.2eV,
    ϵ  = 3.0,
    Nc = NC,
    Nm = 5000,
    a  = 10.0nm,
    σa = 0.0,
    ωM = 2.2eV,
    σM = 0.0,
    Ly = 200nm,
    Lz = 400nm,
    nz = 1,
    ny = 1,
    printout = true
)
h5write(path, "evecs", sys.Uix)
h5write(path, "evals", sys.evals)
h5write(path, "mol_positions", sys.mol_positions)
h5write(path, "mol_energies", sys.mol_energies)
h5write(path, "phot_wavevectors", sys.phot_wavevectors)
h5write(path, "phot_energies", sys.phot_energies)

output("\n\n --- System Summary")
output("\n Molecular Energies (eV)")
output(string(histogram(sys.mol_energies)))
output("\n Cavity Energies")
output(string(lineplot(sys.phot_energies, xlabel="q", ylabel="Energy (ev)")))
output("\n System Eigenenergies")
output(string(scatterplot(sys.evals, xlabel="Eigenstate index", ylabel="Energy (ev)")))

output("\n\n## Time evolution Started")
### TIME EVOLUTION COMPUTATION

output("\nEvolution of molecular states started")
for σi = eachindex(σxvals)

    σx = σxvals[σi]

    # Evolution for q0 = 0
    # Initial state inputs
    output("\n\n• Wavepacket created with σ = $σx nm centered at 25000.0 nm with q = 0.0")
    wvp = create_exciton_wavepacket(25000.0, σx, sys)
    new = similar(wvp)

    for i = eachindex(tvals)
        # Get state at time `t`
        time_propagate!(new, wvp, sys, tvals[i])

        # Save wavepacket in the local basis
        locstates[:,i] = abs2.(sys.Uix * new)
    end
end # σx loop

h5write(path, "wvp", locstates)

output("Done.")

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
