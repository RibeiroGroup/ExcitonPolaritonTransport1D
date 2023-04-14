for d in readdir()
    if isdir(d)
        println("Running $d")
        include(joinpath(@__DIR__, d, "input.jl"))
        mv("log.out", d*"/log.out")
    end
end
~         
