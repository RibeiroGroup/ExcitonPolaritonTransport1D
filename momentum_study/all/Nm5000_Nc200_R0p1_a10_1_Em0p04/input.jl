using PolaritonicSystems
using MKL
using LinearAlgebra
using Unitful
using HDF5
using UnicodePlots
using Dates

output("Initilizing Simulation.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))

for NR = 1:100
    path = joinpath(@__DIR__, "out$NR.h5")
    if isfile(path)
        @info "Skipping realization $NR because output HDF5 file was found."
        continue
    end

    output("Starting Realization $NR. Output file: $path")

    ### INPUTS
    # System inputs
    sys = QuantumWire(
        ΩR = 0.1eV,
        ϵ  = 3.0,
        Nc = 200,
        Nm = 5000,
        a  = 10.0nm,
        σa = 1.0,
        ωM = 2.0eV,
        σM = 0.04,
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
    tvals = 0:0.01:3

    output("\nInitial time:     0.0   ps")
    output("Final time:       {:5.3f} ps", 3)
    output("Time intervals:   {:5.3f} ps", 0.01)
    output("\nTotal Number of steps:  {}", length(tvals))

    h5write(path, "time_range", [0.0, 0.01, 3])
    output("""Time range [t0, δt, tf] saved as "time_range" in $path""")

    output("\nEvolution of molecular states started")
    res = 50
    for σx = [120]

        # Initial state inputs
        output("\n\n• Wavepacket created with σ = $σx nm centered at 25000.0 nm")
        wvp = create_exciton_wavepacket(25000.0, σx, sys)
        new = similar(wvp)

        out = zeros(5000÷res, length(tvals))

        output("\n Displaying results for every 100 steps")
        for i = eachindex(tvals)
            # Get state at time `t`
            time_propagate!(new, wvp, sys, tvals[i])

            # Get molecular wavepacket
            mol_wvp = abs2.(sys.Uix * new)[402:end]

            for (j,r) = enumerate(0:res:(5000-res))
                out[j,i] = sum(mol_wvp[(r+1):(r+res)])
            end
        end

        # Compute data for the barplot using a resolution of 100 molecules
        h5write(path, "$(σx)_wavepacket_bars", out)
        output(""""Wavepackets absolute values saved in $path as "$(σx)_wavepacket_bars" """)
    end # σx loop

    output("Time evolution computations finished.")
    output("Realization $NR finished.")
end # Realization loop

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
