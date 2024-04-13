module NonRelativisticMergers

using Makie
using OrdinaryDiffEq
using PreallocationTools

export FDProblem

# Grids.
export Grid1D
export Grid2D

# Fluid dynamics models.
export Euler

# Reconstruction methods.
export Constant
export MUSCL
export KT

# Riemann solvers.
export NaiveRS
export HLLC

export solve

export plot_euler
export plot_reconstruction

include("types.jl")
include("reconstructor.jl")
include("riemann-solver.jl")

include("1d-euler.jl")

include("plot.jl")

end
