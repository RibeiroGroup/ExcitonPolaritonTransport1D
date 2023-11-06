using PolaritonicSystems
using MKL
using LinearAlgebra
using Unitful
using HDF5
using Dates
using Statistics

output("Initilizing Simulation.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))

# Output file path
path = joinpath(@__DIR__, "out.h5")

# Wave packet initial spread values
σx = 120.0

# Time ranges 
tvals = collect(0:0.01:5)

l = length(tvals)

# Save time ranges into disk
h5write(path, "tvals", tvals)

# Number of realization
TOTAL_NR = 100

# Arrays where we will save array time information 
# Dimensions: Number of time steps x Number of realizations
array_s = zeros(l, TOTAL_NR)
array_pm = zeros(l, TOTAL_NR)

for NR = 1:TOTAL_NR

    output("Starting Realization $NR")

    ### INPUTS
    # System inputs
    sys = QuantumWire(
        ΩR = 0.1eV,
        ϵ  = 3.0,
        Nc = 500,
        Nm = 5000,
        a  = 10.0nm,
        σa = 1.0,
        ωM = 2.0eV,
        σM = 0.2,
        Ly = 200nm,
        Lz = 400nm,
        nz = 1,
        ny = 1,
        printout = true
    )

    output("\n\n## Starting time evolution")

    x0 = 25000
    output("\n\n• Wavepacket created with σ = $σx nm centered at $x0 nm")
    wvp = create_exciton_wavepacket(x0, σx, sys)
    new = similar(wvp)

    ### TIME EVOLUTION COMPUTATION
    for t = eachindex(tvals)
        # Get state at time `t`
        time_propagate!(new, wvp, sys, tvals[t])

        # Compute survival probability
        array_s[t, NR] = abs2(dot(wvp, new))

        # Compute molecular probability
        array_pm[t, NR] = prob_any_mol(new, sys)
    end
end # Realization loop

output("\n\n Computing averages...")
h5write(path, "$(Int(σx))_avg_s", [mean(array_s[t, :]) for t = axes(array_s, 1)])
h5write(path, "$(Int(σx))_std_s", [std(array_s[t, :]) for t = axes(array_s, 1)])

h5write(path, "$(Int(σx))_avg_pm", [mean(array_pm[t, :]) for t = axes(array_pm, 1)])
h5write(path, "$(Int(σx))_std_pm", [std(array_pm[t, :]) for t = axes(array_pm, 1)])
output("Done.")

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
