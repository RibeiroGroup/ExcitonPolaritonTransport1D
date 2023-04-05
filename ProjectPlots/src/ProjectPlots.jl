module ProjectPlots
using HDF5
using Makie
using WGLMakie

export ideal_propagation
export generate_figures

include("data_fetch.jl")
include("propagation_plots.jl")
include("mode_conv_plots.jl")
include("mode_weight_plots.jl")
include("generate.jl")

end # module
