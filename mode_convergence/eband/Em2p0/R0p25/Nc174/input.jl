using PolaritonicSystems
using MKL
using LinearAlgebra
using Unitful
using HDF5
using Statistics
using Dates
using UnicodePlots

output("Initilizing Simulation.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))

path = joinpath(@__DIR__, "out.h5")

# Wavepacket initial spread values
σxvals = [60, 120, 180]

# Simulation time range
tvals = 0:0.01:5

dvals = zeros(length(tvals), length(σxvals))

### INPUTS
# System inputs
sys = QuantumWire(
    ΩR = 0.25eV,
    ϵ  = 3.0,
    Nc = 174,
    Nm = 5000,
    a  = 10.0nm,
    σa = 0.0,
    ωM = 2.0eV,
    σM = 0.0,
    Ly = 200nm,
    Lz = 400nm,
    nz = 1,
    ny = 1,
    printout = true
)

output("\n\n --- System Summary")
output("\n Molecular Energies (eV)")
output(string(histogram(sys.mol_energies)))
output("\n Cavity Energies")
output(string(lineplot(sys.phot_energies, xlabel="q", ylabel="Energy (ev)")))
output("\n System Eigenenergies")
output(string(scatterplot(sys.evals, xlabel="Eigenstate index", ylabel="Energy (ev)")))

output("\n\n## Time evolution Started")
### TIME EVOLUTION COMPUTATION

# Time range

output("\nInitial time:     0.0   ps")
output("Final time:       {:5.3f} ps", tvals[end])
output("Time intervals:   {:5.3f} ps", tvals[2]-tvals[1])
output("\nTotal Number of steps:  {}", length(tvals))

output("\nEvolution of molecular states started")
for σi = 1:3

    σx = σxvals[σi]

    # Initial state inputs
    output("\n\n• Wavepacket created with σ = $σx nm centered at 25000.0 nm with q = 0.0")
    wvp = create_exciton_wavepacket(25000.0, σx, sys)
    new = similar(wvp)

    for i = eachindex(tvals)
        # Get state at time `t`
        time_propagate!(new, wvp, sys, tvals[i])

        # Compute d from Δx²: d = √⟨x²⟩/a
        dvals[i, σi] = sqrt(mean_square_disp(new, 25000.0, sys)) * 0.1
    end
end # σx loop

output("Saving results...")
for σi = 1:3
    σx = σxvals[σi]
    h5write(path, "sx$(σx)_d", dvals[:,σi])
end
output("Done.")

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
