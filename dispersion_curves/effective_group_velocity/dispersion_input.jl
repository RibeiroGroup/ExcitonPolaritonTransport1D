using PolaritonicSystems
using HDF5

# System variables, grouped by (ΩR, Em) i.e., (Rabi splitting, dipole excitation energy)
sysvar = [(0.1, 2.0), (0.2, 2.0), (0.1, 1.8), (0.1, 2.2)]
σxvals = [60, 120, 240, 360, 480]
q0vals = [0, 0.0055]

for (ΩR,EM) in sysvar

    sys = QuantumWire(
        ΩR = ΩR * eV,
        ϵ  = 3.0,
        Nc = 500,
        Nm = 1001,
        a  = 10nm,
        σa = 0,
        ωM = EM * eV,
        σM = 0.0,
        Ly = 200nm,
        Lz = 400nm,
        nz = 1,
        ny = 1
    )

    LP = extract_LP(sys)
    UP = extract_UP(sys)

    effv_LP = effective_group_velocity(LP) ./ 1000 # Convert from nm/ps to μm/ps
    effv_UP = effective_group_velocity(UP) ./ 1000 # Convert from nm/ps to μm/ps

    Rstr = replace(string(ΩR), "." => "p") 
    Emstr = replace(string(EM), "." => "p") 

    h5write("group_vel.h5", "effvg_LP_R$(Rstr)_Em$(Emstr)", effv_LP)
    h5write("group_vel.h5", "effvg_UP_R$(Rstr)_Em$(Emstr)", effv_UP)
    h5write("group_vel.h5", "qvals_R$(Rstr)_Em$(Emstr)", LP.qvals)

    U = PolaritonicSystems.get_q_transformation(sys)

    for σx in σxvals
        for q0 in q0vals

            qstr = replace(string(q0), "." => "p") 
            wvp = create_exciton_wavepacket(5005, σx, sys, q=q0)
            kwvp = U * sys.Uix * wvp
            kwvp = abs2.(kwvp[sys.mol_range])
            m = sortperm(sys.phot_wavevectors)
            h5write("group_vel.h5", "kwvp_q$(qstr)_sm$(σx)_R$(Rstr)_Em$(Emstr)", kwvp[m])
        end
    end
end
