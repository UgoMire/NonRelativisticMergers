abstract type Grid end
abstract type FDModel end
abstract type Reconstructor end
abstract type RiemannSolver end

struct FDProblem{GD<:Grid,MD<:FDModel,RC<:Reconstructor,RS<:RiemannSolver}
    grid::GD
    model::MD
    reconstructor::RC
    riemannsolver::RS
end

struct Grid1D <: Grid
    xmin::Float64
    xmax::Float64
    Nx::Int64
    Δx::Float64
    xl::Vector{Float64}

    function Grid1D(; xmin = 0, xmax = 1, Nx = 100)
        xl = range(xmin, xmax, Nx)
        Δx = xl[2] - xl[1]
        new(xmin, xmax, Nx, Δx, xl)
    end
end

struct Grid2D <: Grid
    xmin::Float64
    xmax::Float64
    Nx::Int64
    Δx::Float64
    xl::Vector{Float64}

    ymin::Float64
    ymax::Float64
    Ny::Int64
    Δy::Float64
    yl::Vector{Float64}

    function Grid2D(; xmin = 0, xmax = 1, ymin = 0, ymax = 1, Nx = 100, Ny = 100)
        xl = range(xmin, xmax, Nx)
        Δx = xl[2] - xl[1]

        yl = range(ymin, ymax, Ny)
        Δy = yl[2] - yl[1]

        new(xmin, xmax, Nx, Δx, xl, ymin, ymax, Ny, Δy, yl)
    end
end

"""
    Euler(; γ = 5 / 3) <: FDModel

Represent Euler's equation, with a polytropic equation of state ``P = ρ e (γ - 1)``.
"""
struct Euler <: FDModel
    γ::Float64

    Euler(; γ = 5 / 3) = new(γ)
end

"""
    EulerStaticGravity(; γ = 5 / 3) <: FDModel

Represent Euler's equation in a static external gravitational potential, with a polytropic 
equation of state ``P = ρ e (γ - 1)``.
"""
struct EulerStaticGravity <: FDModel
    γ::Float64

    EulerStaticGravity(; γ = 5 / 3) = new(γ)
end

"""
    EulerSelfGravity(; γ = 5 / 3, G = 1) <: FDModel

Represent Euler's equation with a polytropic equation of state ``P = ρ e (γ -
1)`` subjected to its own gravitational potential.

The gravitational potential is obtained by solving the Poisson equation `∇²ϕ(x)
= 4π G ρ(x)` at each time step.
"""
struct EulerSelfGravity <: FDModel
    γ::Float64
    G::Float64

    EulerSelfGravity(; γ = 5 / 3, G = 1) = new(γ, G)
end

struct Constant <: Reconstructor end
struct MUSCL <: Reconstructor
    κ::Float64

    MUSCL(κ = 1 / 3) = new(κ)
end
struct KT <: Reconstructor
    θ::Float64

    KT(θ = 2) = new(θ)
end

struct NaiveRS <: RiemannSolver end
struct HLLC <: RiemannSolver end
struct Roe <: RiemannSolver end