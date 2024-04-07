module NonRelativisticMergers

using DifferentialEquations
using Makie

export Grid1D

export Euler1D
export Constant

export MUSCL
export NaiveRS

export solveup

export plot_euler

include("grid.jl")
include("reconstruction.jl")
include("riemann-solver.jl")

include("1d-euler.jl")

include("plot.jl")

end
