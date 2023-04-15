using PolaritonicSystems
using MKL
using Unitful
using HDF5
using Dates

output("Initilizing Simulation.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))

path = joinpath(@__DIR__, "out.h5")

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

h5write(path, "positions", sys.mol_positions)

# Initial state inputs
wvp = create_exciton_wavepacket(5000.0, 60, sys) 

output("\n\n## Time evolution Started")
### TIME EVOLUTION COMPUTATION
# Time range - in fs
tvals = [0, 50, 100, 200, 300, 400]

for t = tvals
    # Get state at time `t`
    new = time_propagate(wvp, sys, t/1000)

    # Convert to local basis
    loc = sys.Uix * new

    # Get |⟨n|ψ⟩|²
    exc = abs2.(loc[sys.mol_range])

    # Save results
    h5write(path, "$(t)fs_wvp", exc)
end
output("Time evolution computations finished.")
output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
