using NonRelativisticMergers
using GLMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

gd = Grid2D(; xmin = -10, xmax = 10, ymin = -10, ymax = 10, Nx = 2^7, Ny = 2^7)

# ρ0 = [1.0 for x in gd.xl, y in gd.yl]
# ρ0 = [0.01 + 10 * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# ρ0 = [0.1 + 10 * exp(-5 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# ρ0 = [0.1 + 10 * exp(-10 * ((x - 1)^2 + (y - 1)^2)) +
#       10 * exp(-10 * ((x + 1)^2 + (y + 1)^2))
#       for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + 10 * exp(-10 * ((x + 1)^2 + (y - 1)^2)) +
#       10 * exp(-10 * ((x + 1)^2 + (y + 1)^2)) +
#       10 * exp(-10 * ((x - 1)^2 + (y + 1)^2)) +
#       40 * exp(-10 * ((x - 1)^2 + (y - 1)^2))
#       for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-100 * x^2) for x in gd.xl, y in gd.yl]
# ρ0 = [-0.2 < x < 0.2 && -0.2 < y < 0.2 ? 2.0 : 1.0 for x in gd.xl, y in gd.yl]
# ρ0 = [1.0 + exp(-10 * (x^2 + y^2 - 1)^2) for x in gd.xl, y in gd.yl]

# vx0 = [0.0 for x in gd.xl, y in gd.yl]
# vx0 = [0.5 for x in gd.xl, y in gd.yl]
# vx0 = [exp(-10 * (x^2 + y^2)) < 0.01 ? 0 : 0.5 for x in gd.xl, y in gd.yl]
# vx0 = [1.0 - 10y * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# vx0 = [-10 * y * exp(-5 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]

# vy0 = [-0.5 for x in gd.xl, y in gd.yl]
# vy0 = [exp(-10 * (x^2 + y^2)) < 0.01 ? 0 : -0.5 for x in gd.xl, y in gd.yl]
# vy0 = [1.0 + 10x * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# vy0 = [10 * x * exp(-5 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]

# P0 = [1.0 for x in gd.xl, y in gd.yl]
# P0 = [1.0 + 10 * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
# P0 = [1.0 + 10 * exp(-10 * ((x - 1)^2 + (y - 1)^2)) +
#       10 * exp(-10 * ((x + 1)^2 + (y + 1)^2))
#       for x in gd.xl, y in gd.yl]
# P0 = [1.0 +
#       10 * (exp(-10 * ((x + 1)^2 + (y - 1)^2)) +
#        exp(-10 * ((x + 1)^2 + (y + 1)^2)) +
#        exp(-10 * ((x - 1)^2 + (y + 1)^2)) +
#        exp(-10 * ((x - 1)^2 + (y - 1)^2)))
#       for x in gd.xl, y in gd.yl]
# P0 = [1.0 + exp(-100 * x^2) for x in gd.xl, y in gd.yl]
# P0 = [1.0 + exp(-10 * (x^2 + y^2 - 1)^2) for x in gd.xl, y in gd.yl]

# reconstructor = Constant()
reconstructor = KT()

# riemannsolver = NaiveRS()
riemannsolver = HLLC()

# model = Euler()
model = EulerSelfGravity(; γ = 6 / 3, G = 2, ϵ = 0.01)

prob = FDProblem(gd, model, reconstructor, riemannsolver)

(; ρ0, vx0, vy0, P0) = NonRelativisticMergers.get_initial_condition(
    prob,
    [
        (; x = -1.5, y = 1.5, height = 5, width = 0.3, v = 2, vangle = -3π / 4),
        (; x = 1.5, y = -1.5, height = 5, width = 0.3, v = 2, vangle = π / 4),
    ],
)

tspan = (0, 50.0)

heatmap(ρ0)
heatmap(vx0)

##
@time sol = solve(prob, ρ0, vx0, vy0, P0, tspan);

# plot_euler(prob, sol; type = :heatmap)
plot_euler(prob, sol; Na = 2)

# record_euler(prob, sol, "rot.mp4")

# NonRelativisticMergers.record_euler(prob, sol, "test.mp4"; trange = 0:0.1:15)

# NonRelativisticMergers.plot_euler2(
#     prob, sol; cmap = :lipari, record_file = "test.gif", trange = 0:1:10)

## Benchmarking.
prob = FDProblem(Grid2D(), EulerSelfGravity(), KT(), HLLC())

(; grid) = prob
(; Nx, Ny) = grid

u0 = NonRelativisticMergers.setup_initial_state(prob, ρ0, vx0, vy0, P0)
du = copy(u0)
p = (;
    prob,
    fv_cache = NonRelativisticMergers.setup_finite_volume_cache(prob),
    fft_cache = NonRelativisticMergers.setup_fft_cache(grid, @view u0[1, :, :]),
)

b = @allocated NonRelativisticMergers.euler2d_self_gravity!(du, u0, p, 0.0)

# b = @benchmark NonRelativisticMergers.euler2d_self_gravity!(du, u0, p, 0.0)