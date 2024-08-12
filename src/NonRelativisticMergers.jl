module NonRelativisticMergers

using FFTW
using Makie
using OrdinaryDiffEq
using PreallocationTools
using LinearAlgebra

export FDProblem

# Grids.
export Grid1D
export Grid2D

# Fluid dynamics models.
export Euler
export EulerStaticGravity
export EulerSelfGravity

# Reconstruction methods.
export Constant
export MUSCL
export KT

# Riemann solvers.
export NaiveRS
export HLLC

export solve

export plot_euler
export plot_euler2d
export plot_reconstruction

include("types.jl")

include("reconstructor.jl")

include("riemann-solver/naive.jl")
include("riemann-solver/hllc-1d.jl")
include("riemann-solver/hllc-2d.jl")

include("poisson.jl")

include("model/1d-euler.jl")
include("model/2d-euler.jl")

include("plot.jl")

end
