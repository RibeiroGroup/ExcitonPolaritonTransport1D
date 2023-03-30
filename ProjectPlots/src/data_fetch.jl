"""
    Convert a number to a string replacing the decimal number to a pnumber

# Example

julia> pnumber(0.1)
"0p1" 
"""
function pnumber(N::Number)
    if isinteger(N)
        return string(Int(N))
    end
    return replace(string(N), '.'=>'p')
end

"""
    Get mean square displacement data for ideal propagation (no disorder), along with time range.
"""
function get_ideal_x2(;Nm::Int, Nc::Int, ΩR, a)
    fname  = "Nm$(Nm)/Nc$(Nc)/R$(pnumber(ΩR))/a$a/Em2p0_sx60/out.h5"
    path = joinpath(@__DIR__,"../../propagation_study/ideal/"*fname)
    
    x2 = h5read(path, "mean_square_disp")
    (t0, δt, tf) = h5read(path, "time_range")

    return t0:δt:tf, x2
end

function max_cavity_energy(;Nm::Int, Nc::Int, a)
    fname  = "Nm$(Nm)/Nc$(Nc)/R0p05/a$a/Em2p0_sx60/log.out"
    path = joinpath(@__DIR__,"../../propagation_study/ideal/"*fname)

    max_e = 0.0
    open(path) do io
        m = match(r"Maximum Cavity energy:\s+?(\d+?\.\d+)", read(io, String))
        if isnothing(m)
            println("Couldn't read 'Maximum Cavity energy' from the output file at Nm=$Nm, Nc=$Nc, ΩR=$ΩRm, a=$a")
        else
            max_e = parse(Float64, m.captures[1])
        end
    end

    return max_e
end