using PolaritonicSystems
using MKL       # Fast matrix operations

# Plotting libraries
using Makie
using CairoMakie

# System inputs
sys = QuantumWire(
    ΩR = 0.05eV,
    ϵ  = 3.0,
    Nc = 100,
    Nm = 1000,
    a  = 10nm,
    σa = 0,
    ωM = 2.0eV,
    σM = 0,
    Ly = 200nm,
    Lz = 400nm,
    nz = 1,
    ny = 1,
)

# Initial state inputs
wvp = create_exciton_wavepacket(5000.0, 120, sys) 
# Time inverval in picoseconds
tvals = 0:0.01:1
# Preallocate array for photonic probabilities
pvals = zeros(length(tvals))

for i in eachindex(tvals)
    t = tvals[i]

    # Get the time-evolved wave packet
    new = time_propagate(wvp, sys, t)
    # Convert to local basis
    loc = sys.Uix * new
    # Compute probability for each photon mode
    prob = abs2.(loc[sys.phot_range])

    pvals[i] = sum(prob)
end

# Plotting commands. See Makie documentation for more details. 
fig = Figure()
ax = Axis(fig[1,1], xlabel="Time (ps)", ylabel="Probability")
lines!(ax, tvals, pvals, color=:goldenrod, label="Photon", linewidth=3)
lines!(ax, tvals, 1 .- pvals, color=:royalblue, label="Exciton", linewidth=3)
axislegend(ax)
save("example.png", fig)
