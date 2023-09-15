using CairoMakie

"""
    Generate all articles figures in the given path.
"""
function generate_figures(;path=".", quality=5, SI=false, maxnum=100)

    CairoMakie.activate!()

    # Scan module for all functions named figX, where X ∈ N
    println("Looking for plot functions 👀")
    plot_functions = []
    fstr = SI ? "SI_fig" : "fig"
    i = 1
    while i ≤ maxnum
        if isdefined(Nanophotonics, Symbol(fstr*"$i"))
            push!(plot_functions, eval(Symbol(fstr*"$i")))
        elseif i > 1 # Make sure to check for at least fig2
            break
        end
        i += 1
    end

    nf = length(plot_functions)
    println("$nf functions found!\n")

    # Call plot functions and save figures
    println("📜 Generating article figures...")
    for (i, func) in enumerate(plot_functions)
        print("($i/$nf) $(string(func)) ")
        save(joinpath(path, string(func)*".png"), func(), px_per_unit=quality)
        println("✅")
    end
    println("Done 🎉")
end