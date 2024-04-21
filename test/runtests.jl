using NonRelativisticMergers
using Test

verbose = true

@testset verbose = verbose "1d Euler - Hydrostatic flow is time independent" begin
    grid = Grid1D(; xmin = 0, xmax = 1, Nx = 100)

    ρ0l = ones(grid.Nx)
    v0l = zeros(grid.Nx)
    p0l = zeros(grid.Nx)

    tspan = (0, 1)

    @testset "- $reconstructor, $riemannsolver" for reconstructor in (Constant(), KT()),
        riemannsolver in (NaiveRS(), HLLC())

        prob = FDProblem(grid, Euler(), reconstructor, riemannsolver)
        sol = solve(prob, ρ0l, v0l, p0l, tspan)

        ρend = sol.u[end][1, :]
        vend = sol.u[end][2, :] ./ ρend
        pend = @. (sol.u[end][3, :] - 0.5 * ρend * vend^2) * (prob.model.γ - 1)

        @test ρ0l ≈ ρend
        @test v0l ≈ vend
        @test p0l ≈ pend
    end
end

@testset verbose = verbose "2d flux functions are consistent" begin
    import NonRelativisticMergers: Fflux, Gflux, flux

    ρ = 2.2300
    vx = 8.9900
    vy = 7.2219
    P = 5.4320

    prob = FDProblem(Grid2D(), Euler(), KT(), HLLC())

    @test all(Fflux(prob, ρ, vx, vy, P) .≈ flux(prob, ρ, vx, vy, P, (; x = 1, y = 0)))
    @test all(Gflux(prob, ρ, vx, vy, P) .≈ flux(prob, ρ, vx, vy, P, (; x = 0, y = 1)))
end

@testset verbose = verbose "2d Euler - Hydrostatic flow is time independent" begin
    grid = Grid2D(; xmin = 0, xmax = 1, ymin = 0, ymax = 1, Nx = 100, Ny = 100)

    ρ = ones(grid.Nx, grid.Ny)
    vx = zeros(grid.Nx, grid.Ny)
    vy = zeros(grid.Nx, grid.Ny)
    p = ones(grid.Nx, grid.Ny)

    tspan = (0, 1)

    @testset "- $reconstructor, $riemannsolver" for reconstructor in (Constant(), KT()),
        riemannsolver in (NaiveRS(), HLLC())

        prob = FDProblem(grid, Euler(), reconstructor, riemannsolver)
        sol = solve(prob, ρ, vx, vy, p, tspan)

        ρend = sol.u[end][1, :, :]
        vxend = sol.u[end][2, :, :] ./ ρend
        vyend = sol.u[end][3, :, :] ./ ρend
        pend =
            @. (sol.u[end][4, :, :] - 0.5 * ρend * (vxend^2 + vyend^2)) * (prob.model.γ - 1)

        @test ρ ≈ ρend
        @test vx ≈ vxend
        @test vy ≈ vyend
        @test p ≈ pend
    end
end

@testset verbose = verbose "HLLC Riemann solver - 2d and 1d implementation agree" begin
    import NonRelativisticMergers: setup_initial_state, reconstruct!, hllc_riemann_solver

    grid2d = Grid2D(; xmin = -2, xmax = 2, ymin = -2, ymax = 2, Nx = 100, Ny = 100)
    grid1d = Grid1D(; xmin = -2, xmax = 2, Nx = 100)

    # Initial condition that are y independent.
    ρ0 = [1.0 + exp(-100 * x^2) for x in grid2d.xl, y in grid2d.yl]
    vx0 = [0.0 for x in grid2d.xl, y in grid2d.yl]
    vy0 = [0.0 for x in grid2d.xl, y in grid2d.yl]
    P0 = [1.0 + exp(-100 * x^2) for x in grid2d.xl, y in grid2d.yl]

    @testset "- $reconstructor" for reconstructor in (Constant(), KT())
        problem1d = FDProblem(Grid1D(), Euler(), reconstructor, HLLC())
        problem2d = FDProblem(Grid2D(), Euler(), reconstructor, HLLC())

        u0 = setup_initial_state(problem2d, ρ0, vx0, vy0, P0)

        # Reconstruct the initial state at the cell boundaries. 
        wreconstructed = zeros(4, grid2d.Nx, grid2d.Ny, 4)
        reconstruct!(problem2d, wreconstructed, u0)

        ρ = @view wreconstructed[1, :, :, :]
        vx = @view wreconstructed[2, :, :, :]
        vy = @view wreconstructed[3, :, :, :]
        P = @view wreconstructed[4, :, :, :]

        for i in 1:grid2d.Nx, j in 1:grid2d.Ny
            ip = i == grid2d.Nx ? 1 : i + 1

            rhoL = ρ[i, j, 2]
            uL = vx[i, j, 2]
            vL = vy[i, j, 2]
            pL = P[i, j, 2]

            rhoR = ρ[ip, j, 1]
            uR = vx[ip, j, 1]
            vR = vy[ip, j, 1]
            pR = P[ip, j, 1]

            n = (; x = 1, y = 0)

            flux2d = hllc_riemann_solver(problem2d, rhoL, uL, vL, pL, rhoR, uR, vR, pR, n)
            flux1d = hllc_riemann_solver(problem1d, rhoL, uL, pL, rhoR, uR, pR)

            @test flux2d[1] ≈ flux1d[1]
            @test flux2d[2] ≈ flux1d[2]
            @test flux2d[3] ≈ 0
            @test flux2d[4] ≈ flux1d[3]
        end
    end
end