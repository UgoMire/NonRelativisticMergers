using NonRelativisticMergers
using Test

@testset verbose = true "1d Euler - Hydrostatic flow is time independent" begin
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