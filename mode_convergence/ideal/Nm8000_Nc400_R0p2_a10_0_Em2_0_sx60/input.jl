using MQED
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
    ΩR = 0.2eV,
    ϵ  = 3.0,
    Nc = 400,
    Nm = 8000,
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
wvp = create_exciton_wavepacket(40000.0, 60, sys) 

output("\n\n## Time evolution Started")
### TIME EVOLUTION COMPUTATION
# Time range
tvals = 0:0.01:50

output("\nInitial time:     0.0   ps")
output("Final time:       {:5.3f} ps", 50)
output("Time intervals:   {:5.3f} ps", 0.01)
output("\nTotal Number of steps:  {}", length(tvals))

# Time evolution
Δx2 = zeros(length(tvals))

output("\n Displaying results for every 100 steps")
for i = eachindex(tvals)
    # Get state at time `t`
    new = time_propagate(wvp, sys, tvals[i])

    # Compute Δx²
    Δx2[i] = mean_square_disp(new, 40000.0, sys)

    if (i-1) % 100 == 0
        output("\nExciton profile at time = $(tvals[i]) ps")
        output(string(lineplot(1:length(sys.mol_positions), get_exciton_prob(new, sys), xlabel="Molecule", ylabel="Probability")))
    end
end

output("Time evolution computations finished.")

# Save output
path = joinpath(@__DIR__, "out.h5")
h5write(path, "mean_square_disp", Δx2)
h5write(path, "time_range", [tvals[1], tvals[2]-tvals[1], tvals[end]])

output("""Mean square displacement saved as "mean_square_disp" in $path""")
output("""Time range [t0, δt, tf] saved as "time_range" in $path""")
output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
