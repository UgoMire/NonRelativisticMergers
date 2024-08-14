using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid1D(; xmin = -5, xmax = 5, Nx = 500)

# ρ0l = ones(gd.Nx)
# ρ0l = 0.1 .+ map(x -> 1.8 * exp(-100 * (x - 0.0)^2), gd.xl)
ρ0l = 1 .+ map(x -> 1.5 * exp(-100 * (x - 2)^2) + exp(-100 * (x + 2)^2), gd.xl)
# ρ0l = [0.4 < x < 0.6 ? 0.5 : 0.05 for x in gd.xl]

v0l = zeros(gd.Nx)
# v0l = map(x -> 0.0 * (-exp(-10 * (x - 2)^2) + exp(-10 * (x + 2)^2)), gd.xl)
# v0l = ones(gd.Nx)
# v0l = map(x -> 0.01 * exp(-100 * (x - 0.5)^2), gd.xl)

# p0l = ones(gd.Nx)
# p0l = zeros(gd.Nx)
# p0l = 0.1 .+ map(x -> 3 * exp(-100 * (x - 0.0)^2), gd.xl)
# p0l = 1 .+ map(x -> 0.1 * exp(-100 * (x - 0.5)^2), gd.xl)
p0l = 1 .+ map(x -> 1.5 * exp(-100 * (x - 2)^2) + exp(-100 * (x + 2)^2), gd.xl)
# p0l = [0.4 < x < 0.6 ? 0.5 : 0.1 for x in gd.xl]

tspan = (0, 3)

# reconstructor = Constant()
# reconstructor = MUSCL()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

# model = Euler()
model = EulerSelfGravity(; γ = 4 / 3, G = 10, ϵ = 9)

prob = FDProblem(gd, model, reconstructor, riemannsolver)

sol = solve(prob, ρ0l, v0l, p0l, tspan)

plot_euler(prob, sol)
