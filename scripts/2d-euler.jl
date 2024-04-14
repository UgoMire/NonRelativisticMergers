using NonRelativisticMergers
using GLMakie

gd = Grid2D(; xmin = -1, xmax = 1, ymin = -1, ymax = 1, Nx = 100, Ny = 100)

ρ0 = [1.0 + exp(-100 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
vx0 = [0.0 for x in gd.xl, y in gd.yl]
vy0 = [0.0 for x in gd.xl, y in gd.yl]
P0 = [1.0 + 0.01 * exp(-100 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]

tspan = (0, 0.15)

reconstructor = Constant()
# reconstructor = KT()

riemannsolver = NaiveRS()
# riemannsolver = HLLC()

prob = FDProblem(gd, Euler(), reconstructor, riemannsolver)

sol = solve(prob, ρ0, vx0, vy0, P0, tspan);

plot_euler2d(prob, sol)

## Benchmarking.
using PreallocationTools

u0 = NonRelativisticMergers.setup_initial_state(prob, ρ0, vx0, vy0, P0)

wstore = zeros(4, prob.grid.Nx, prob.grid.Ny, 4)
xfluxstore = zeros(4, prob.grid.Nx, prob.grid.Ny)
yfluxstore = zeros(4, prob.grid.Nx, prob.grid.Ny)

du = copy(u0)
p = (; prob, wstore, xfluxstore, yfluxstore)

@btime NonRelativisticMergers.euler2d!(du, u0, p, tspan[1])