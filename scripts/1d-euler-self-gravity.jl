using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid1D(; xmin = -5, xmax = 5, Nx = 100)

# reconstructor = Constant()
# reconstructor = MUSCL()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

# model = Euler()
model = EulerSelfGravity(; γ = 5 / 3, G = 1, ϵ = 4)

prob = FDProblem(gd, model, reconstructor, riemannsolver)

# ρ0l = ones(gd.Nx)
# ρ0l = 1 .+ map(x -> 1 * exp(-100 * (x - 0.0)^2), gd.xl)
ρ0l = 1 .+ map(x -> 1.5 * exp(-20 * (x - 2)^2) + exp(-20 * (x + 2)^2), gd.xl)
# ρ0l = [0.4 < x < 0.6 ? 0.5 : 0.05 for x in gd.xl]

v0l = zeros(gd.Nx)
# v0l = map(x -> 0.0 * (-exp(-10 * (x - 2)^2) + exp(-10 * (x + 2)^2)), gd.xl)
# v0l = ones(gd.Nx)
# v0l = map(x -> 0.01 * exp(-100 * (x - 0.5)^2), gd.xl)

p0l = ones(gd.Nx)
# p0l = zeros(gd.Nx)
# p0l = 1 .+ map(x -> 0 * exp(-100 * (x - 0.0)^2), gd.xl)
# p0l = 1 .+ map(x -> 0.1 * exp(-100 * (x - 0.5)^2), gd.xl)
# p0l = 1 .+ map(x -> 1.5 * exp(-100 * (x - 2)^2) + exp(-100 * (x + 2)^2), gd.xl)
# p0l = [0.4 < x < 0.6 ? 0.5 : 0.1 for x in gd.xl]

tspan = (0, 10)

sol = solve(prob, ρ0l, v0l, p0l, tspan)

plot_euler(prob, sol)

## Benchmarking the RHS.
u0 = NonRelativisticMergers.setup_initial_state(prob, ρ0l, v0l, p0l)
du = similar(u0)

fv_cache = NonRelativisticMergers.setup_finite_volume_cache(prob)

fft_cache = NonRelativisticMergers.setup_fft_cache(gd)

p = (; prob, fv_cache, fft_cache)

@benchmark NonRelativisticMergers.euler1d_self_gravity!(du, u0, p, 0)
# @time NonRelativisticMergers.euler1d_self_gravity!(du, u0, p, 0)

## Hydrostatic initial conditions.
gd = Grid1D(; xmin = -3, xmax = 3, Nx = 500)

model = EulerSelfGravity(; γ = 4 / 3, G = 1, ϵ = 10)

prob = FDProblem(gd, model)

# ρ0l = ones(gd.Nx)
# ρ1l = [0.5 + 1.0 * exp(-10 * (x - 1)^2) for x in gd.xl]
ρ1l = [0.2 - tanh(10 * (x - 1 - 0.6)) + tanh(10 * (x - 1 + 0.6)) for x in gd.xl]
p1l = get_static_pressure(prob, ρ1l)
# v1l = [-0.1 * exp(-5 * (x - 0.5)^2) for x in gd.xl]
v1l = -0.1 * ones(gd.Nx)
# v1l = zeros(gd.Nx)
# v1l = [-0.01 - 0.1 * (-tanh(10 * (x - 1 - 1)) + tanh(10 * (x - 1 + 1))) for x in gd.xl]

ρ2l = [0.5 + 1.8 * exp(-10 * (x + 1)^2) for x in gd.xl]
p2l = get_static_pressure(prob, ρ2l)
# v2l = [0.1 * exp(-5 * (x + 0.5)^2) for x in gd.xl]

ρ0l = ρ1l
p0l = p1l .+ 0.5
v0l = v1l

# ρ0l = ρ1l .+ ρ2l
# p0l = p1l .+ p2l .+ 0.5
# v0l = v1l .+ v2l

tspan = (0, 5)

sol = solve(prob, ρ0l, v0l, p0l, tspan)

plot_euler(prob, sol)
