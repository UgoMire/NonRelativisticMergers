using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid2D(; xmin = -2, xmax = 2, ymin = -2, ymax = 2, Nx = 100, Ny = 100)

# ρ0 = [1.0 for x in gd.xl, y in gd.yl]
ρ0 = [1.0 + 10 * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-100 * x^2) for x in gd.xl, y in gd.yl]
# ρ0 = [-0.2 < x < 0.2 && -0.2 < y < 0.2 ? 2.0 : 1.0 for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-10 * (x^2 + y^2 - 1)^2) for x in gd.xl, y in gd.yl]

vx0 = [0.0 for x in gd.xl, y in gd.yl]
# vx0 = [1.0 for x in gd.xl, y in gd.yl]
# vx0 = [1.0 - 10y * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# vx0 = [-y / (2π * (x^2 + y^2)) for x in gd.xl, y in gd.yl]

vy0 = [0.0 for x in gd.xl, y in gd.yl]
# vy0 = [1.0 + 10x * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# vy0 = [x / (2π * (x^2 + y^2)) for x in gd.xl, y in gd.yl]

# P0 = [1.0 for x in gd.xl, y in gd.yl]
P0 = [1.0 + 10 * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# P0 = [1.0 + exp(-100 * x^2) for x in gd.xl, y in gd.yl]
# P0 = [1.0 + exp(-10 * (x^2 + y^2 - 1)^2) for x in gd.xl, y in gd.yl]

tspan = (0, 4.0)

# reconstructor = Constant()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

# model = Euler()
model = EulerSelfGravity(; γ = 5 / 3, G = 1, ϵ = 3)

prob = FDProblem(gd, model, reconstructor, riemannsolver)

@time sol = solve(prob, ρ0, vx0, vy0, P0, tspan);

plot_euler(prob, sol; type = :surface)

## Benchmarking.
prob = FDProblem(Grid2D(), EulerSelfGravity(), KT(), HLLC())

(; grid) = prob
(; Nx, Ny) = grid

u0 = NonRelativisticMergers.setup_initial_state(prob, ρ0, vx0, vy0, P0)
du = copy(u0)
p = (;
    prob,
    fv_cache = NonRelativisticMergers.setup_finite_volume_cache(prob),
    fft_cache = NonRelativisticMergers.setup_fft_cache(grid, @view u0[1, :, :])
)

@time NonRelativisticMergers.euler2d_self_gravity!(du, u0, p, 0.0)
# @benchmark NonRelativisticMergers.euler2d_self_gravity!(du, u0, p, 0.0)