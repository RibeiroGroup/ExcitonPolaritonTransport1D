for d in ["Nc79", "Nc127", "Nc183"]
    if isdir(d)
        println("Running $d")
        include(joinpath(@__DIR__, d, "input.jl"))
        mv("log.out", d*"/log.out")
    end
end
~         
