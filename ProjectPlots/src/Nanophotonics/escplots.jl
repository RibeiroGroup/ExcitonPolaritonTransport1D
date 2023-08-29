function escp_over_time!(ax::Axis; ΩR::Float64=0.1, σx::Int=120, σM::Float64=0.005, color=Makie.wong_colors()[1], label="", fit=false)

    # Get time values
    r1 = 0:0.005:0.5
    r2 = 1.0:0.5:5
    tvals = vcat(r1, r2)
    #tvals = r1

    escp = zeros(length(tvals))

    for i in eachindex(escp)
        escp[i] = wvp_escape_probability(ΩR, σx, σM, i)
    end

    #escp = escp ./ tvals


    scatter!(ax, tvals, escp, color=color, label=label)
    lines!(ax, tvals, escp, color=color, label=label)
end


function wvp_escape_probability(ΩR, σx, σM, idx)

	Rstr = replace(string(ΩR), "." => "p") 
    sm = replace(string(σM), "." => "p") 
    path = joinpath(@__DIR__, "../../../propagation_study/disorder/Nm5000_Nc500_a10_Em2p0/R$Rstr/sm$sm/out.h5")

    wvp = h5read(path, "$(Int(σx))_avg_wvp", (:, idx))
    wvp = wvp ./ sum(wvp)
    wvp0 = h5read(path, "$(Int(σx))_avg_wvp", (:, 1))
    wvp0 = wvp0 ./ sum(wvp0)

    # Get range to be subtracted (i.e., range of bins where 0.98 of the initial probability lies)
    n = 0
    s = 0
    r = nothing
    while s < 0.98
        n += 1
        r = (50-n+1):(50+n)
        s = sum(wvp0[r])
        if n > 20
            break
        end
    end

    # Print probability of exciton at the "boundary" for the last time step
    if idx == 110
        println(wvp[1] + wvp[end])
        println(r)
    end

    return (sum(wvp) - sum(wvp[r])) 
end
