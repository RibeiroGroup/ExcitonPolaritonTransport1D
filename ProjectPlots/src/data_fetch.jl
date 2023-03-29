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