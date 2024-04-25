using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid2D(; xmin = -2, xmax = 2, ymin = -2, ymax = 2, Nx = 200, Ny = 200)
gd = Grid2D(; xmin = -2, xmax = 2, ymin = -2, ymax = 2, Nx = 200, Ny = 200)

ρ0 = [1.0 for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-100 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-100 * x^2) for x in gd.xl, y in gd.yl]
ρ0 = [1.0 for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-100 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-100 * x^2) for x in gd.xl, y in gd.yl]
# ρ0 = [-0.2 < x < 0.2 && -0.2 < y < 0.2 ? 2.0 : 1.0 for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-100 * (x^2 + y^2 - 1)^2) for x in gd.xl, y in gd.yl]

# vx0 = [0.0 for x in gd.xl, y in gd.yl]
vx0 = [1.0 for x in gd.xl, y in gd.yl]
# vx0 = [1.0 - 10y * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# vx0 = [-y / (2π * (x^2 + y^2)) for x in gd.xl, y in gd.yl]

vy0 = [0.0 for x in gd.xl, y in gd.yl]
# vy0 = [1.0 + 10x * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# vy0 = [x / (2π * (x^2 + y^2)) for x in gd.xl, y in gd.yl]

# P0 = [1.0 for x in gd.xl, y in gd.yl]
P0 = [1.0 + 10 * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# P0 = [1.0 + exp(-100 * x^2) for x in gd.xl, y in gd.yl]
# P0 = [1.0 + exp(-100 * (x^2 + y^2 - 1)^2) for x in gd.xl, y in gd.yl]

tspan = (0, 1.150)

# reconstructor = Constant()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()
# riemannsolver = NaiveRS()
riemannsolver = HLLC()

prob = FDProblem(gd, Euler(), reconstructor, riemannsolver)

sol = solve(prob, ρ0, vx0, vy0, P0, tspan);

plot_euler2d(prob, sol)

## Benchmarking.
# reconstructor = Constant()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

prob = FDProblem(gd, Euler(), reconstructor, riemannsolver)

u0 = NonRelativisticMergers.setup_initial_state(prob, ρ0, vx0, vy0, P0)

wstore = zeros(4, prob.grid.Nx, prob.grid.Ny, 4)
xfluxstore = zeros(4, prob.grid.Nx, prob.grid.Ny)
yfluxstore = zeros(4, prob.grid.Nx, prob.grid.Ny)

du = copy(u0)
p = (; prob, wstore, xfluxstore, yfluxstore)

@btime NonRelativisticMergers.euler2d!(du, u0, p, tspan[1])
##

@code_warntype NonRelativisticMergers.euler2d!(du, u0, p, tspan[1])
@code_warntype NonRelativisticMergers.hllc_riemann_solver(
    prob,
    (; ρ = 1, vx = 1, vy = 1, P = 1),
    (; ρ = 1, vx = 1, vy = 1, P = 1),
    (; x = 1, y = 0),
)
@code_warntype NonRelativisticMergers.hllc_riemann_solver(
    prob,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    (; x = 1, y = 0),
)