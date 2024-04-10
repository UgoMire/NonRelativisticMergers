using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid1D(; xmin = 0, xmax = 1, Nx = 100)

# ρ0l = ones(gd.Nx)
# ρ0l = 0.1 .+ map(x -> 0.1 * exp(-100 * (x - 0.5)^2), gd.xl)
ρ0l = [0.4 < x < 0.6 ? 0.1 : 0.01 for x in gd.xl]

v0l = zeros(gd.Nx)
# v0l = ones(gd.Nx)
# v0l = map(x -> 0.01 * exp(-100 * (x - 0.5)^2), gd.xl)

# p0l = ones(gd.Nx)
# p0l = zeros(gd.Nx)
# p0l = 0.1 .+ map(x -> 0.2 * exp(-100 * (x - 0.5)^2), gd.xl)
p0l = [0.4 < x < 0.6 ? 0.1 : 0.01 for x in gd.xl]

tspan = (0, 0.5)

# reconst = Constant()
# reconst = MUSCL()
reconst = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

# model = Euler1D(Constant(), NaiveRS(), ρ0l, v0l, p0l)
model = Euler1D(ρ0l, v0l, p0l; reconst, riemannsolver)

sol = solveup(gd, model, tspan)

plot_euler(sol, gd)

##
u0 = zeros(3, gd.Nx)
u0[1, :] = ρ0l
u0[2, :] = v0l
u0[3, :] = p0l
u0 = sol.u[end]

plot_reconstruction(gd, model, u0)