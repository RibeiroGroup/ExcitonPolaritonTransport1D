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

# Total number of realizations 
NR_T = 25

# Number of cavity modes
NC = 200

# Wavepacket initial spread values
σxvals = [60, 120, 180]

MW_q0 = zeros((2NC+1), length(σxvals), NR_T)
MW_qnz = zeros((2NC+1), length(σxvals), NR_T)
EARLY_MW_q0 = zeros((2NC+1), length(σxvals), NR_T)
EARLY_MW_qnz = zeros((2NC+1), length(σxvals), NR_T)

# Time range
tvals = 0:0.005:5

phot_cont_q0 = zeros(length(tvals), length(σxvals), NR_T)
phot_cont_qnz = zeros(length(tvals), length(σxvals), NR_T)

for NR in 1:NR_T
    output("Starting realization $NR")

    ### INPUTS
    # System inputs
    sys = QuantumWire(
        ΩR = 0.3eV,
        ϵ  = 3.0,
        Nc = NC,
        Nm = 5000,
        a  = 10.0nm,
        σa = 1.0,
        ωM = 2.2eV,
        σM = 0.005,
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
    
    
    output("\nInitial time:     0.0   ps")
    output("Final time:       {:5.3f} ps", tvals[end])
    output("Time intervals:   {:5.3f} ps", tvals[1]-tvals[2])
    output("\nTotal Number of steps:  {}", length(tvals))
    
    output("\nEvolution of molecular states started")
    for σi = 1:3
    
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
            mode_weight = abs2.(sys.Uix * new)[sys.phot_range]
            MW_q0[:,σi,NR] .+= mode_weight
            phot_cont_q0[i,σi,NR] = sum(mode_weight)

            # Save early results
            if i < 102
                EARLY_MW_q0[:,σi,NR] += mode_weight
            end
        end

        # Evolution for q0 > 0
        output("\n\n• Wavepacket created with σ = $σx nm centered at 25000.0 nm with q = 0.005 cm⁻¹")
        # Initial state inputs
        wvp = create_exciton_wavepacket(25000.0, σx, sys, q=0.005654866776461627)
        new = similar(wvp)

        for i = eachindex(tvals)
            # Get state at time `t`
            time_propagate!(new, wvp, sys, tvals[i])

            # Save wavepacket in the local basis
            mode_weight = abs2.(sys.Uix * new)[sys.phot_range]
            MW_qnz[:,σi,NR] += mode_weight
            phot_cont_qnz[i,σi,NR] = sum(mode_weight)

            # Save early results
            if i < 102
                EARLY_MW_qnz[:,σi,NR] += mode_weight
            end

        end
    end # σx loop
    output("Realization $NR finished.")
end # NR loop

# Divide it by the number of time steps to get an average over time
MW_q0 ./= length(tvals)
MW_qnz ./= length(tvals)

EARLY_MW_q0 ./= 101
EARLY_MW_qnz ./= 101

output("Computing averages...")
for σi = 1:3

    σx = σxvals[σi]
    
    h5write(path, "$(σx)_avg_mode_weight", [mean(MW_q0[i,σi,:]) for i = axes(MW_q0, 1)])
    h5write(path, "$(σx)_std_mode_weight", [mean(MW_q0[i,σi,:]) for i = axes(MW_q0, 1)])

    h5write(path, "nzq_$(σx)_avg_mode_weight", [mean(MW_qnz[i,σi,:]) for i = axes(MW_qnz, 1)])
    h5write(path, "nzq_$(σx)_std_mode_weight", [mean(MW_qnz[i,σi,:]) for i = axes(MW_qnz, 1)])

    h5write(path, "EARLY_$(σx)_avg_mode_weight", [mean(EARLY_MW_q0[i,σi,:]) for i = axes(EARLY_MW_q0, 1)])
    h5write(path, "EARLY_$(σx)_std_mode_weight", [mean(EARLY_MW_q0[i,σi,:]) for i = axes(EARLY_MW_q0, 1)])

    h5write(path, "EARLY_nzq_$(σx)_avg_mode_weight", [mean(EARLY_MW_qnz[i,σi,:]) for i = axes(EARLY_MW_qnz, 1)])
    h5write(path, "EARLY_nzq_$(σx)_std_mode_weight", [mean(EARLY_MW_qnz[i,σi,:]) for i = axes(EARLY_MW_qnz, 1)])

    h5write(path, "q0_$(σx)_phot_cont", [mean(phot_cont_q0[i,σi,:]) for i = axes(phot_cont_q0, 1)])
    h5write(path, "nzq_$(σx)_phot_cont", [mean(phot_cont_qnz[i,σi,:]) for i = axes(phot_cont_qnz, 1)])
end
output("Done.")

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
