function get_nabla_operator(grid::Grid1D)
    (; Nx) = grid

    ∇x = Matrix(Tridiagonal(ones(Nx - 1), zeros(Nx), -ones(Nx - 1)))
    ∇x[1, Nx] = 1
    ∇x[Nx, 1] = -1

    return ∇x
end

function get_static_pressure(prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any}, ρ)
    (; grid) = prob
    (; Nx) = grid

    ϕ = solve_poisson(prob, ρ)

    ∇x = get_nabla_operator(grid)
    ∇ϕ = ∇x * ϕ

    # lprob = LinearProblem(∇x, -ρ .* ∇ϕ)
    # lsol = LinearSolve.solve(lprob, SVDFactorization())
    # P = lsol.u

    nlprob = NonlinearProblem((P, _) -> ∇x * P + ρ .* ∇ϕ, zeros(Nx))
    nlsol = NonlinearSolve.solve(nlprob, NewtonRaphson())
    P = nlsol.u

    return P
end