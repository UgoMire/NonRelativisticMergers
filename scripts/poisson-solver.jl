using NonRelativisticMergers
using GLMakie
using FFTW

import NonRelativisticMergers: solve_poisson!

(; xl, Nx, Δx) = gd = Grid1D(; xmin = -1, xmax = 1, Nx = 250)
# model = EulerSelfGravity(gd; G = 1, ϵ = 1)
prob = FDProblem(gd, EulerSelfGravity(), KT(), HLLC())

α = 1e2
ρ = map(x -> exp(-α * x^2), xl)

ϕ = zeros(Nx)

@btime solve_poisson!(prob, ϕ, ρ)

f = lines(xl, ρ)
lines!(xl, ϕ)
f

## Benchmarking the solver.
cache = (;
    kx = 2π * fftfreq(Nx, 1 / Δx),
    ρ_complex = zeros(ComplexF64, Nx),
    potentialstore_complex = zeros(ComplexF64, Nx),
    ρhat_cache = zeros(ComplexF64, Nx),
    uhat_cache = zeros(ComplexF64, Nx),
    planned_fft = plan_fft(ρ),
    planned_ifft = plan_ifft(ρ),
)

@btime solve_poisson!(prob, ϕ, ρ, cache)
# solve_poisson!(prob, ϕ, ρ, cache)

##
@btime solve_poisson!(prob, ϕ, ρ, cache)