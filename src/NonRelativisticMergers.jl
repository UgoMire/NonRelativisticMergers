module NonRelativisticMergers

using DifferentialEquations
using Makie
using PreallocationTools

export FDProblem

# Grids.
export Grid1D

# Fluid dynamics models.
export Euler1D

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
