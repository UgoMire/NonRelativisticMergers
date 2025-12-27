function get_nabla_operator(grid::Grid2D)
    (; Nx, Δx, Ny, Δy) = grid

    ∇x = Matrix(Tridiagonal(-ones(Nx - 1), zeros(Nx), ones(Nx - 1)))
    ∇x[1, Nx] = -1
    ∇x[Nx, 1] = 1
    ∇x ./= 2Δx

    ∇y = Matrix(Tridiagonal(-ones(Ny - 1), zeros(Ny), ones(Ny - 1)))
    ∇y[1, Ny] = -1
    ∇y[Ny, 1] = 1
    ∇y ./= 2Δy

    return (; ∇x, ∇y)
end

function get_static_pressure(prob::FDProblem{Grid2D,EulerSelfGravity,<:Any,<:Any}, ρ)
    (; grid) = prob
    (; Nx, Δx, Ny, Δy) = grid

    ϕ = solve_poisson(prob, ρ)

    (; ∇x, ∇y) = get_nabla_operator(grid)

    ∇xϕ = ∇x * ϕ
    ∇yϕ = ϕ * ∇y'

    lhs = -∇x * (ρ .* ∇xϕ) .- (ρ .* ∇yϕ) * ∇y'

    kx = 2π * rfftfreq(Nx, 1 / Δx)
    ky = 2π * fftfreq(Ny, 1 / Δy)

    lhshat = rfft(lhs)

    Phat = zeros(ComplexF64, Nx ÷ 2 + 1, Ny)

    for (i, kxi) in enumerate(kx), (j, kyj) in enumerate(ky)
        if kxi == 0 && kyj == 0
            Phat[i, j] = 0
        else
            Phat[i, j] = -lhshat[i, j] / (kxi^2 + kyj^2)
        end
    end

    P = irfft(Phat, Nx)

    Pmin = minimum(P)
    @. P += -Pmin .+ 1

    return P
end

function get_initial_condition(prob::FDProblem{Grid2D,EulerSelfGravity,<:Any,<:Any}, stars)
    (; grid) = prob
    (; xl, Nx, yl, Ny) = grid

    ρ0 = 0.01 * ones(Nx, Ny)
    vx0 = zeros(Nx, Ny)
    vy0 = zeros(Nx, Ny)

    for (; x, y, height, width, v, vangle) in stars
        ρ0 .+= [height * exp(-((xp - x)^2 + (yp - y)^2) / 2width^2) for xp in xl, yp in yl]

        vx = v * cos(vangle)
        vx0 .+= [
            exp(-((xp - x)^2 + (yp - y)^2) / 2width^2) < 1e-4 ? 0.0 : vx for xp in xl,
            yp in yl
        ]

        vy = v * sin(vangle)
        vy0 .+= [
            exp(-((xp - x)^2 + (yp - y)^2) / 2width^2) < 1e-4 ? 0.0 : vy for xp in xl,
            yp in yl
        ]
    end

    P0 = get_static_pressure(prob, ρ0)

    return (; ρ0, vx0, vy0, P0)
end