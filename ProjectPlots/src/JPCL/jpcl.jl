"""
    Functions in this module were used to process and plot results presented in the Journal of Physical Chemistry Letters.
    DOI https://pubs.acs.org/doi/10.1021/acs.jpclett.3c01082
"""
module JPCL

using ProjectPlots
using LaTeXStrings
using Makie
using WGLMakie
using Statistics
using HDF5

    include("../data_fetch.jl")
    include("propagation_plots.jl")
    include("mode_conv_plots.jl")
    include("mode_weight_plots.jl")
    include("generate.jl")

end