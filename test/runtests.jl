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

@testset verbose = verbose "2d Euler - Hydrostatic flow is time independent" begin
    grid = Grid2D(; xmin = 0, xmax = 1, ymin = 0, ymax = 1, Nx = 100, Ny = 100)

    ρ = ones(grid.Nx, grid.Ny)
    vx = zeros(grid.Nx, grid.Ny)
    vy = zeros(grid.Nx, grid.Ny)
    p = ones(grid.Nx, grid.Ny)

    tspan = (0, 1)

    @testset "- $reconstructor, $riemannsolver" for reconstructor in (Constant(), KT()),
        riemannsolver in (NaiveRS(),)

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