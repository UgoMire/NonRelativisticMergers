using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid1D(; xmin = -1, xmax = 2, Nx = 600)

# ρ0l = ones(gd.Nx)
ρ0l = 0.1 .+ map(x -> 1 * exp(-100 * (x - 0.5)^2), gd.xl)
# ρ0l = [0.4 < x < 0.6 ? 0.5 : 0.01 for x in gd.xl]

v0l = zeros(gd.Nx)
# v0l = ones(gd.Nx)
# v0l = map(x -> 0.01 * exp(-100 * (x - 0.5)^2), gd.xl)

# p0l = ones(gd.Nx)
# p0l = zeros(gd.Nx)
p0l = 0.1 .+ map(x -> 5.8 * exp(-100 * (x - 0.5)^2), gd.xl)
# p0l = 1 .+ map(x -> 1 * exp(-100 * (x - 0.5)^2), gd.xl)
# p0l = [0.4 < x < 0.6 ? 0.8 : 0.02 for x in gd.xl]

tspan = (0, 1)

# reconstructor = Constant()
# reconstructor = MUSCL()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

# model = EulerSelfGravity()
model = Euler()

prob = FDProblem(gd, model, reconstructor, riemannsolver)

sol = solve(prob, ρ0l, v0l, p0l, tspan)

plot_euler(prob, sol)

##
u0 = zeros(3, gd.Nx)
u0[1, :] = ρ0l
u0[2, :] = v0l
u0[3, :] = p0l
u0 = sol.u[end]

plot_reconstruction(gd, model, u0)

##
# reconst = Constant()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

prob = FDProblem(gd, Euler(), reconstructor, riemannsolver)

u0 = zeros(3, gd.Nx)

γ = 5 / 3

@. u0[1, :] = ρ0l
@. u0[2, :] = v0l * ρ0l
@. u0[3, :] = p0l / (γ - 1) + 1 / 2 * ρ0l * v0l^2

du = copy(u0)
wstore = zeros(3, gd.Nx, 2)
fluxstore = zeros(3, gd.Nx)

p = (; prob, wstore, fluxstore)

@benchmark NonRelativisticMergers.euler1d!(du, u0, p, tspan[1])
# @time NonRelativisticMergers.euler1d!(du, u0, p, tspan[1])
##
@code_warntype NonRelativisticMergers.euler1d!(
    du,
    u0,
    (; gd, model, wstore, fluxstore),
    tspan[1]
)
