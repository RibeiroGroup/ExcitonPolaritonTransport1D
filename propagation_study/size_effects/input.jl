using PolaritonicSystems
using MKL
using Unitful
using HDF5
using UnicodePlots
using Dates

output("Initilizing Simulation.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))

### INPUTS
# System inputs
sys = QuantumWire(
    ΩR = 0.1eV,
    ϵ  = 3.0,
    Nc = 800,
    Nm = 1000,
    a  = 10nm,
    σa = 0,
    ωM = 2.0eV,
    σM = 0,
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

# Initial state inputs
wvp = create_exciton_wavepacket(5000, 60, sys) 

output("\n\n## Time evolution Started")
### TIME EVOLUTION COMPUTATION
# Time range
tvals = 0:0.001:5

output("\nInitial time:     0.0   ps")
output("Final time:       {:5.3f} ps", 5)
output("Time intervals:   {:5.3f} ps", 0.001)
output("\nTotal Number of steps:  {}", length(tvals))

# Time evolution

# Array for photon probability
phot = zeros(length(tvals))

# Array for probability of finding the exciton on the last 100 molecules
pend = zeros(length(tvals))

for i = eachindex(tvals)
    # Get state at time `t`
    new = time_propagate(wvp, sys, tvals[i])

    # Convert to local basis
    Ploc = abs2.(sys.Uix * new)

    # Compute Δx²
    phot[i] = sum(Ploc[sys.phot_range])
    pend[i] = sum(Ploc[sys.mol_range][1:100]) + sum(Ploc[sys.mol_range][900:1000])

end
output("Time evolution computations finished.")

# Save output
path = joinpath(@__DIR__, "out.h5")
h5write(path, "pphot", phot)
h5write(path, "pend", pend)
h5write(path, "time_range", [tvals[1], tvals[2]-tvals[1], tvals[end]])

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
