struct Euler1D{RC<:Reconstruction,RS<:RiemannSolver}
    reconst::RC
    riemannsolver::RS

    γ::Float64

    function Euler1D(;
        γ = 5 / 3,
        reconst::RC = Constant(),
        riemannsolver::RS = NaiveRS(),
    ) where {RC<:Reconstruction,RS<:RiemannSolver}
        new{RC,RS}(reconst, riemannsolver, γ)
    end
end

function flux(model, ρ, v, p)
    γ = model.γ

    return (ρ * v, ρ * v^2 + p, (p / (γ - 1) + 1 / 2 * ρ * v^2 + p) * v)
end

function flux(model, i, ρ, v, p)
    γ = model.γ

    if i == 1
        return ρ * v
    elseif i == 2
        return ρ * v^2 + p
    elseif i == 3
        return (p / (γ - 1) + 1 / 2 * ρ * v^2 + p) * v
    else
        throw(ArgumentError("Invalid index $i"))
    end
end

function euler1d!(du, u, p, t)
    (; gd, model, wstore, fluxstore) = p

    reconstruct!(wstore, model.reconst, gd, u)

    solve_riemann_problem!(fluxstore, wstore, model.riemannsolver, model, gd)

    for j = 1:3, i = 1:gd.Nx
        im = i == 1 ? gd.Nx : i - 1

        du[j, i] = -(fluxstore[j, i] - fluxstore[j, im]) / gd.Δx
    end
end

function solveup(ρ0l, v0l, p0l, gd::Grid1D, model::Euler1D, tspan)
    u0 = zeros(3, gd.Nx)

    @. u0[1, :] = ρ0l
    @. u0[2, :] = v0l * ρ0l
    @. u0[3, :] = p0l / (model.γ - 1) + 1 / 2 * ρ0l * v0l^2

    wstore = zeros(3, gd.Nx, 2)
    fluxstore = zeros(3, gd.Nx)

    prob = ODEProblem(euler1d!, u0, tspan, (; gd, model, wstore, fluxstore))

    sol = solve(
        prob,
        Tsit5(),
        # TRBDF2(autodiff = false),
        # QNDF(autodiff = false),
        abstol = 1e-14,
        reltol = 1e-14,
        progress = true,
        progress_steps = 100,
    )

    return sol
end