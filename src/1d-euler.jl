struct Euler1D{RC<:Reconstruction,RS<:RiemannSolver}
    reconst::RC
    riemannsolver::RS

    ρ0l::Vector{Float64}
    v0l::Vector{Float64}
    p0l::Vector{Float64}

    function Euler1D(
        ρ0l,
        v0l,
        p0l;
        reconst::RC = Constant(),
        riemannsolver::RS = NaiveRS(),
    ) where {RC<:Reconstruction,RS<:RiemannSolver}
        new{RC,RS}(reconst, riemannsolver, ρ0l, v0l, p0l)
    end
end

function flux(ρ, v, p)
    e = p / (ρ * (5 / 3 - 1))

    return [ρ * v, ρ * v^2 + p, (ρ * e + 1 / 2 * ρ * v^2 + p) * v]
end

function euler1d!(du, u, p, t)
    (; gd, model) = p

    w_reconstruct = reconstruct(model.reconst, gd, u)

    F = solve_riemann(model.riemannsolver, gd, w_reconstruct)

    for i = 1:gd.Nx
        im = i == 1 ? gd.Nx : i - 1
        du[:, i] = -(F[:, i] - F[:, im]) / gd.Δx
    end
end

function solveup(gd::Grid1D, model::Euler1D, tspan)
    u0 = zeros(3, gd.Nx)

    @. u0[1, :] = model.ρ0l
    @. u0[2, :] = model.v0l
    @. u0[3, :] = model.p0l

    prob = ODEProblem(euler1d!, u0, tspan, (; gd, model))

    sol = solve(
        prob,
        Tsit5(),
        # TRBDF2(),
        # QNDF(autodiff = false),
        abstol = 1e-8,
        reltol = 1e-8,
        progress = true,
        progress_steps = 500,
    )

    return sol
end