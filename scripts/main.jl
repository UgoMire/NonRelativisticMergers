using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid1D(; xmin = 0, xmax = 1, Nx = 500)

ρ0l = ones(gd.Nx)
# ρ0l = map(x -> 0.1 * exp(-100 * (x - 0.5)^2), gd.xl)
v0l = zeros(gd.Nx)
# v0l = map(x -> 0.01 * exp(-100 * (x - 0.5)^2), gd.xl)
# p0l = ones(gd.Nx)
p0l = map(x -> 0.01 * exp(-100 * (x - 0.5)^2), gd.xl)

tspan = (0, 3)

reconst = Constant()
# reconst = MUSCL()

# model = Euler1D(Constant(), NaiveRS(), ρ0l, v0l, p0l)
model = Euler1D(ρ0l, v0l, p0l; reconst)

sol = solveup(gd, model, tspan)

plot_euler(sol, gd)