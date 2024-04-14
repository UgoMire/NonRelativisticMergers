using NonRelativisticMergers
using GLMakie

gd = Grid2D(; xmin = -1, xmax = 1, ymin = -1, ymax = 1, Nx = 100, Ny = 100)

ρ0l = [1.0 + exp(-100 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
v0l = [0.0 for x in gd.xl, y in gd.yl]
p0l = [0.0 for x in gd.xl, y in gd.yl]

tspan = (0, 5)

reconstructor = Constant()
reconstructor = KT()

riemannsolver = NaiveRS()
riemannsolver = HLLC()

prob = FDProblem(gd, Euler(), reconstructor, riemannsolver)

# sol = solve(prob, ρ0l, v0l, p0l, tspan)

## Benchmarking.
using PreallocationTools

u0 = NonRelativisticMergers.get_initial_conditions(prob, ρ0l, v0l, p0l)

wstore = DiffCache(zeros(3, prob.gd.Nx, prob.gd.Ny, 2))
fluxstore = DiffCache(zeros(3, prob.gd.Nx, prob.gd.Ny))

du = copy(u0)
p = (; prob, wstore, fluxstore)

@btime NonRelativisticMergers.euler2d!(du, u0, p, tspan[1])