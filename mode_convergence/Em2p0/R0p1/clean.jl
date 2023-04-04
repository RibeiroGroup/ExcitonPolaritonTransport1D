using HDF5   

function cleanup()
    fnames = ["sm0p005", "sm0p02", "sm0p05", "sm0p1"]

    for f in fnames
        for Nc in ["Nc0", "Nc1", "Nc10", "Nc100", "Nc20", "Nc200", "Nc35", "Nc400", "Nc5", "Nc50", "Nc500", "Nc75", "Nc800"]

            inpath = joinpath(@__DIR__, "../R0p1/$f/$Nc/out.h5")
            outpath = joinpath(@__DIR__, "$f/$Nc/out.h5")

            println("Current Path: $inpath")
            for σx = [60, 120, 180]
                for NofR in [20, 40, 60, 80, 100]

                    d = h5read(inpath, "NR_$(NofR)_sm$(σx)_avg_mode_weight")
                    std = h5read(inpath, "NR_$(NofR)_sm$(σx)_std_mode_weight")

                    h5write(outpath, "NR_$(NofR)_sm$(σx)_avg_d", d)
                    h5write(outpath, "NR_$(NofR)_sm$(σx)_std_d", std)
                end
            end
        end
    end
end
