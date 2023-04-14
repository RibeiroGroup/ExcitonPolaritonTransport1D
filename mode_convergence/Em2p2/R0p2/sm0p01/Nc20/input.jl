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

# Total number of realizations 
NR_T = 100

# Wavepacket initial spread values
σxvals = [60, 120, 180]

# Simulation time range
tvals = 0:0.01:5

mol_amplitudes = zeros(length(tvals), length(σxvals), NR_T)

for NR in 1:NR_T
    output("Starting realization $NR")

    ### INPUTS
    # System inputs
    sys = QuantumWire(
        ΩR = 0.2eV,
        ϵ  = 3.0,
        Nc = 20,
        Nm = 5000,
        a  = 10.0nm,
        σa = 1.0,
        ωM = 2.2eV,
        σM = 0.01,
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
            mol_amplitudes[i, σi, NR] = sqrt(mean_square_disp(new, 25000.0, sys)) * 0.1
        end
    end # σx loop
    output("Realization $NR finished.")
end # NR loop

output("Computing averages...")
for σi = 1:3
    σx = σxvals[σi]
    # Compute averages for various numbers of realizations
    for NofR in [20, 40, 60, 80, 100]

        r = 1:NofR
        
        h5write(path, "NR_$(NofR)_sm$(σx)_avg_d", [mean(mol_amplitudes[i,σi,r]) for i = eachindex(tvals)])
        h5write(path, "NR_$(NofR)_sm$(σx)_std_d", [std(mol_amplitudes[i,σi,r]) for i = eachindex(tvals)])
    end
end
output("Done.")

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
