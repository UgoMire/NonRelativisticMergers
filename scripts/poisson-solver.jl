using NonRelativisticMergers
using GLMakie
using FFTW

import NonRelativisticMergers: solve_poisson!, solve_poisson

(; xl, Nx, Δx) = gd = Grid1D(; xmin = -1, xmax = 1, Nx = 500)
prob = FDProblem(gd, EulerSelfGravity(), KT(), HLLC())

α = 1e2
ρ = map(x -> exp(-α * x^2), xl)

@time ϕ = solve_poisson(prob, ρ; fft_type = :fft);

# f = lines(xl, ρ)
# lines!(xl, ϕ)
# f

## Benchmarking the solver.
cache = NonRelativisticMergers.setup_fft_cache(gd);

@btime solve_poisson!(prob, ϕ, ρ, cache)

##
@btime solve_poisson!(prob, ϕ, ρ, cache)

## Solve in 2d.
(; xl, Nx, Δx, yl, Ny, Δy) = gd = Grid2D(;
    xmin = -1, xmax = 1, Nx = 100, ymin = -1, ymax = 1, Ny = 100)
prob = FDProblem(gd, EulerSelfGravity(; ϵ = 5, G = 1), KT(), HLLC())

α = 1e2
ρ = [exp(-α * (x^2 + y^2)) + exp(-α * ((x - 0.5)^2 + (y - 0.5)^2)) for x in xl, y in xl]

@time ϕ = NonRelativisticMergers.solve_poisson(prob, ρ; fft_type = :fft);

f = Figure(; size = (1000, 400))
ax1 = Axis3(f[1, 1])
surface!(ax1, xl, yl, ρ)
ax2 = Axis3(f[1, 2])
surface!(ax2, xl, yl, ϕ)
f

## Benchmarking the solver.
cache = NonRelativisticMergers.setup_fft_cache(gd);
ϕ = zeros(Nx, Ny)

@btime solve_poisson!(prob, ϕ, ρ, cache)