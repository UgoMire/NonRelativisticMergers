using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid1D(; xmin = -1, xmax = 2, Nx = 500)

# ρ0l = ones(gd.Nx)
ρ0l = 1 .+ map(x -> 1 * exp(-100 * (x - 0.5)^2), gd.xl)
# ρ0l = [0.4 < x < 0.6 ? 0.5 : 0.05 for x in gd.xl]

v0l = zeros(gd.Nx)
# v0l = ones(gd.Nx)
# v0l = map(x -> 0.01 * exp(-100 * (x - 0.5)^2), gd.xl)

# p0l = ones(gd.Nx)
# p0l = zeros(gd.Nx)
# p0l = 0.1 .+ map(x -> 1.8 * exp(-100 * (x - 0.5)^2), gd.xl)
p0l = 1 .+ map(x -> 0.1 * exp(-100 * (x - 0.5)^2), gd.xl)
# p0l = [0.4 < x < 0.6 ? 0.5 : 0.1 for x in gd.xl]

tspan = (0, 1)

# reconstructor = Constant()
# reconstructor = MUSCL()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

model = EulerSelfGravity()

prob = FDProblem(gd, model, reconstructor, riemannsolver)

sol = solve(prob, ρ0l, v0l, p0l, tspan)

plot_euler(prob, sol)
