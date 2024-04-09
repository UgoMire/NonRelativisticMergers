module NonRelativisticMergers

using DifferentialEquations
using Makie

export Grid1D

export Euler1D

# Reconstruction methods.
export Constant
export MUSCL
export KT

# Riemann solvers.
export NaiveRS
export HLLC

export solveup

export plot_euler
export plot_reconstruction

include("grid.jl")
include("reconstruction.jl")
include("riemann-solver.jl")

include("1d-euler.jl")

include("plot.jl")

end
