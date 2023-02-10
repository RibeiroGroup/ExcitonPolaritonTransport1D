using PolaritonicSystems
using MKL
using LinearAlgebra
using Unitful
using HDF5
using UnicodePlots
using Dates

output("Initilizing Simulation.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))

path = joinpath(@__DIR__, "out.h5")

### INPUTS
# System inputs
sys = QuantumWire(
    ΩR = 0.1eV,
    ϵ  = 3.0,
    Nc = 200,
    Nm = 5000,
    a  = 10.0nm,
    σa = 0.0,
    ωM = 1.9eV,
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

h5write(path, "eigenenergies", sys.evals)
output("""System energies saved as "eigenenergies" in $path""")

output("\n\n## Time evolution Started")
### TIME EVOLUTION COMPUTATION

# Time range
tvals = 0:0.005:5

output("\nInitial time:     0.0   ps")
output("Final time:       {:5.3f} ps", 3)
output("Time intervals:   {:5.3f} ps", 0.01)
output("\nTotal Number of steps:  {}", length(tvals))

h5write(path, "time_range", [0.0, 0.01, 3])
output("""Time range [t0, δt, tf] saved as "time_range" in $path""")

output("\nEvolution of molecular states started")
for σx = [60, 120, 180]

    # Initial state inputs
    output("\n\n• Wavepacket created with σ = $σx nm centered at 25000.0 nm with q = 0.0")
    wvp = create_exciton_wavepacket(25000.0, σx, sys)
    new = similar(wvp)

    mw_sum = zeros(length(sys.phot_range))

    output("\n Displaying results for every 100 steps")
    for i = eachindex(tvals)
        # Get state at time `t`
        time_propagate!(new, wvp, sys, tvals[i])

        # Save wavepacket in the local basis
        mw_sum += abs2.(sys.Uix * new)[sys.phot_range]
    end
    out = mw_sum ./ length(tvals)
    out = out ./ maximum(out)

    h5write(path, "$(σx)_phot_weight", out)
    output(""""Mode weight computed.\n""")
end # σx loop

output("Time evolution computations finished.")

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
