using PolaritonicSystems
using MKL
using LinearAlgebra
using Unitful
using UnicodePlots
using HDF5
using Dates
using Statistics

output("Initilizing Simulation.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))

# Output file path
path = joinpath(@__DIR__, "out.h5")

# Wave packet initial spread values
σxvals = [60.0, 120.0, 240.0, 360.0, 480.0]

# Time ranges 
# r1 for both d and wvp
r1 = 0:0.005:0.5
# r2 for wvp
r2 = 1.0:0.5:5
# r3 for d
r3 = 0.51:0.01:5

l1 = length(r1)
l2 = length(r2)
l3 = length(r3)

wvp_tvals = zeros(l1+l2)
wvp_tvals[1:l1] .= r1
wvp_tvals[(l1+1):end] .= r2

d_tvals = zeros(l1+l3)
d_tvals[1:l1] .= r1
d_tvals[(l1+1):end] .= r3

# Save time ranges into disk
h5write(path, "wvp_tvals", wvp_tvals)
h5write(path, "d_tvals", d_tvals)

# Make sure that wvp_tvals ⊂  d_tvals
@assert all([t in d_tvals for t in wvp_tvals])

# Number of realization
TOTAL_NR = 100

# Arrays where we will save array time information 
# Dimensions: Compressed Index x Time steps x Realizations x # σx
# Exciton spread
array_d = zeros(length(d_tvals), TOTAL_NR, 5)
# Molecular wavepacket (grouped by 50 mols)
array_wvp = zeros(100, length(wvp_tvals), TOTAL_NR, 5)
# Zero momentum photon mode 
array_zpm = zeros(length(wvp_tvals), TOTAL_NR, 5)
# Positive momentum photon modes  (grouped in 10 modes)
array_ppm = zeros(50, length(wvp_tvals), TOTAL_NR, 5)
# Negative momentum photon modes  (grouped in 10 modes)
array_npm = zeros(50, length(wvp_tvals), TOTAL_NR, 5)

for NR = 1:TOTAL_NR

    output("Starting Realization $NR")

    ### INPUTS
    # System inputs
    sys = QuantumWire(
        ΩR = 0.3eV,
        ϵ  = 3.0,
        Nc = 500,
        Nm = 5000,
        a  = 10.0nm,
        σa = 1.0,
        ωM = 2.0eV,
        σM = 0.3,
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

    output("\n\n## Starting time evolution")

    ### TIME EVOLUTION COMPUTATION

    # Initial state inputs

    for (iσx, σx) in enumerate(σxvals)
        x0 = 25000
        output("\n\n• Wavepacket created with σ = $σx nm centered at $x0 nm")
        wvp = create_exciton_wavepacket(x0, σx, sys)
        new = similar(wvp)
        uncnew = similar(new)


        # Short time evolution
        for t = eachindex(wvp_tvals)
            # Get state at time `t`
            time_propagate!(new, wvp, sys, wvp_tvals[t])

            # Get new state in the uncoupled basis
            uncnew .= sys.Uix * new

            # Save zero momentum photon component

            # Save momentum components in groups of 50
            for (k,i) in enumerate(sys.phot_range)
                if k == 1
                    array_zpm[t, NR, iσx] = abs2(uncnew[i])
                elseif k%2 == 0
                    # Compute group index
                    keven = k ÷ 2 
                    idx = (keven-1) ÷ 10 + 1
                    array_ppm[idx, t, NR, iσx] += abs2(uncnew[i])
                else
                    # Compute group index
                    kodd = (k-1) ÷ 2 
                    idx = (kodd-1) ÷ 10 + 1
                    array_npm[idx, t, NR, iσx] += abs2(uncnew[i])
                end
            end

            # Save molecular part of the wvp |Cm|² in groups of 50
            for (k,i) in enumerate(sys.mol_range)
                # Compute group index
                idx = (k-1) ÷ 50 + 1
                array_wvp[idx, t, NR, iσx] += abs2(uncnew[i])
            end

            # Since theres an overlap between time steps for wvp computation and d computation
            # we also compute d while computing wvp
            # Loop through molecules
            x2 = 0.0
            Pmol = 0.0
            for (i,x) = zip(sys.mol_range, sys.mol_positions)
                # Contribution of each molecule is P * (x-x0)^2, where P is the probability of the state
                # being localized at that molecule and x is its position
                x2 += abs2(uncnew[i]) * (x - x0)^2

                # Probabilty of finding a molecule excited
                Pmol += abs2(uncnew[i])
            end
            # Renormalize (conditional probability) with the probability of finding the molecular exciton
            # Compute d (sqrt(x2) divided by average intermolecular distance)
            
            # Since the indexes of d_tvals and wvp_tvals don't match for all entries, we need to make
            # sure we are putting the computed d in the right place. Do to that, we use `findfirst`
            td = findfirst(x->x==wvp_tvals[t], d_tvals)
            array_d[td, NR, iσx] = 0.1*sqrt(x2 / Pmol) 
        end

        # Long time evolutions
        for t in eachindex(d_tvals)
            # Skip if it has been filled during the wvp computation
            if array_d[t,NR,iσx] == 0.0
                # Get state at time `t`
                time_propagate!(new, wvp, sys, d_tvals[t])
                array_d[t, NR, iσx] = 0.1*sqrt(mean_square_disp(new, x0, sys))
            end
        end
    end # σx loop
end # Realization loop

output("\n\n Computing averages...")
for (iσx, σx) in enumerate(σxvals)
    # d 
    h5write(path, "$(Int(σx))_avg_d", [mean(array_d[t, :, iσx]) for t = axes(array_d, 1)])
    h5write(path, "$(Int(σx))_std_d", [std(array_d[t, :, iσx]) for t = axes(array_d, 1)])
    
    # Photon modes (q = 0)
    h5write(path, "$(Int(σx))_avg_zpm", [mean(array_zpm[t, :, iσx]) for t = axes(array_zpm, 1)])
    h5write(path, "$(Int(σx))_std_zpm", [std(array_zpm[t, :, iσx]) for t = axes(array_zpm, 1)])
    
    # Photon modes (q > 0)
    h5write(path, "$(Int(σx))_avg_ppm", [mean(array_ppm[i, t, :, iσx]) for i = axes(array_ppm, 1), t = axes(array_ppm,2)])
    h5write(path, "$(Int(σx))_std_ppm", [std(array_ppm[i, t, :, iσx]) for i = axes(array_ppm, 1), t = axes(array_ppm,2)])
    
    # Photon modes (q > 0)
    h5write(path, "$(Int(σx))_avg_npm", [mean(array_npm[i, t, :, iσx]) for i = axes(array_npm, 1), t = axes(array_npm,2)])
    h5write(path, "$(Int(σx))_std_npm", [std(array_npm[i, t, :, iσx]) for i = axes(array_npm, 1), t = axes(array_npm,2)])
    
    # Wavepackets
    h5write(path, "$(Int(σx))_avg_wvp", [mean(array_wvp[i, t, :, iσx]) for i = axes(array_wvp, 1), t = axes(array_wvp,2)])
    h5write(path, "$(Int(σx))_std_wvp", [std(array_wvp[i, t, :, iσx]) for i = axes(array_wvp, 1), t = axes(array_wvp,2)])
end

output("Done.")

output("Exiting code sucessfully.")
output(Dates.format(now(), "U, dd yyyy - HH:MM:SS"))
