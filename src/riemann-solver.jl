function solve_riemann_problem!(
    prob::FDProblem{Grid1D,<:Any,<:Any,NaiveRS},
    fluxstore,
    wreconstructed,
)
    (; Nx) = prob.grid

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1

        Fminux = flux(
            prob,
            wreconstructed[1, i, 2],
            wreconstructed[2, i, 2],
            wreconstructed[3, i, 2],
        )
        Fplus = flux(
            prob,
            wreconstructed[1, ip, 1],
            wreconstructed[2, ip, 1],
            wreconstructed[3, ip, 1],
        )

        @. fluxstore[:, i] = 1 / 2 * (Fminux + Fplus)
    end
end

function solve_riemann_problem!(
    prob::FDProblem{Grid2D,<:Any,<:Any,NaiveRS},
    xfluxstore,
    yfluxstore,
    w,
)
    (; Nx, Ny) = prob.grid

    ρ = @view w[1, :, :, :]
    vx = @view w[2, :, :, :]
    vy = @view w[3, :, :, :]
    P = @view w[4, :, :, :]

    for i in 1:Nx, j in 1:Ny
        ip = i == Nx ? 1 : i + 1
        jp = j == Ny ? 1 : j + 1

        Fminus = Fflux(prob, ρ[i, j, 2], vx[i, j, 2], vy[i, j, 2], P[i, j, 2])
        Fplus = Fflux(prob, ρ[ip, j, 1], vx[ip, j, 1], vy[ip, j, 1], P[ip, j, 1])

        @. xfluxstore[:, i, j] = 1 / 2 * (Fminus + Fplus)
        # @. xfluxstore[:, i, j] = Fplus
        # @. xfluxstore[:, i, j] = Fminus

        Gminus = Gflux(prob, ρ[i, j, 4], vx[i, j, 4], vy[i, j, 4], P[i, j, 4])
        Gplus = Gflux(prob, ρ[i, jp, 3], vx[i, jp, 3], vy[i, jp, 3], P[i, jp, 3])

        @. yfluxstore[:, i, j] = 1 / 2 * (Gminus + Gplus)
        # @. yfluxstore[:, i, j] = Gplus
        # @. yfluxstore[:, i, j] = Gminus
    end
end

function hllc_wavespeed_estimate(
    prob::FDProblem{<:Any,Euler,<:Any,HLLC},
    rhoL,
    vL,
    pL,
    rhoR,
    vR,
    pR,
)
    (; γ) = prob.model

    cL = sqrt(γ * pL / rhoL)
    cR = sqrt(γ * pR / rhoR)

    HL = pL / rhoL * γ / (γ - 1) + 1 / 2 * vL^2
    HR = pR / rhoR * γ / (γ - 1) + 1 / 2 * vR^2

    Rρ = sqrt(rhoR / rhoL)

    vtilde = (vL + vR * Rρ) / (1 + Rρ)
    Htilde = (HL + HR * Rρ) / (1 + Rρ)
    ctilde = sqrt((γ - 1) * (Htilde - 1 / 2 * vtilde^2))

    SL = min(vL - cL, vtilde - ctilde)
    SR = max(vR + cR, vtilde + ctilde)
    SM =
        (rhoR * vR * (SR - vR) - rhoL * vL * (SL - vL) + pL - pR) /
        (rhoR * (SR - vR) - rhoL * (SL - vL))

    return (; SL, SR, SM)
end

function hllc_riemann_solver(
    prob::FDProblem{Grid1D,Euler,<:Any,HLLC},
    rhoL,
    vL,
    pL,
    rhoR,
    vR,
    pR,
)
    (; γ) = prob.model

    (; SL, SR, SM) = hllc_wavespeed_estimate(prob, rhoL, vL, pL, rhoR, vR, pR)

    if 0 < SL
        return flux(prob, rhoL, vL, pL)
    elseif SL <= 0 < SM
        FL = flux(prob, rhoL, vL, pL)
        EL = pL / (γ - 1) + 1 / 2 * rhoL * vL^2

        UL = (rhoL, rhoL * vL, EL)

        rhoLstar = rhoL * (SL - vL) / (SL - SM)
        pLstar = pL + rhoL * (vL - SL) * (vL - SM)
        rhoLvLstar = (rhoL * vL * (SL - vL) + pLstar - pL) / (SL - SM)

        ELstar = (EL * (SL - vL) - pL * vL + pLstar * SM) / (SL - SM)

        ULstar = (rhoLstar, rhoLvLstar, ELstar)

        return @. FL + SL * (ULstar - UL)
    elseif SM <= 0 <= SR
        FR = flux(prob, rhoR, vR, pR)
        ER = pR / (γ - 1) + 1 / 2 * rhoR * vR^2

        UR = (rhoR, rhoR * vR, ER)

        rhoRstar = rhoR * (SR - vR) / (SR - SM)
        pRstar = pR + rhoR * (vR - SR) * (vR - SM)
        rhoRvRstar = (rhoR * vR * (SR - vR) + pRstar - pR) / (SR - SM)

        ERstar = (ER * (SR - vR) - pR * vR + pRstar * SM) / (SR - SM)

        URstar = (rhoRstar, rhoRvRstar, ERstar)

        return @. FR + SR * (URstar - UR)
    else
        return flux(prob, rhoR, vR, pR)
    end
end

function solve_riemann_problem!(
    prob::FDProblem{Grid1D,Euler,<:Any,HLLC},
    fluxstore,
    wreconstructed,
)
    (; Nx) = prob.grid

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1

        rhoL = wreconstructed[1, i, 2]
        vL = wreconstructed[2, i, 2]
        pL = wreconstructed[3, i, 2]

        rhoR = wreconstructed[1, ip, 1]
        vR = wreconstructed[2, ip, 1]
        pR = wreconstructed[3, ip, 1]

        fluxstore[:, i] .= hllc_riemann_solver(prob, rhoL, vL, pL, rhoR, vR, pR)
    end
end

function hllc_wavespeed_estimate(
    prob::FDProblem{Grid2D,Euler,<:Any,HLLC},
    rhoL,
    uL,
    vL,
    pL,
    rhoR,
    uR,
    vR,
    pR,
    n,
)
    (; γ) = prob.model

    qL = uL * n.x + vL * n.y
    qR = uR * n.x + vR * n.y

    HL = pL / rhoL * γ / (γ - 1) + 1 / 2 * (uL^2 + vL^2)
    HR = pR / rhoR * γ / (γ - 1) + 1 / 2 * (uR^2 + vR^2)

    Rρ = sqrt(rhoR / rhoL)

    u_tilde = (uL + uR * Rρ) / (1 + Rρ)
    v_tilde = (vL + vR * Rρ) / (1 + Rρ)
    q_tilde = u_tilde * n.x + v_tilde * n.y

    H_tilde = (HL + HR * Rρ) / (1 + Rρ)

    cL = sqrt(γ * pL / rhoL)
    cR = sqrt(γ * pR / rhoR)
    c_tilde = sqrt((γ - 1) * (H_tilde - 1 / 2 * (u_tilde^2 + v_tilde^2)))

    SL = min(qL - cL, q_tilde - c_tilde)
    SR = max(qR + cR, q_tilde + c_tilde)
    SM =
        (rhoR * qR * (SR - qR) - rhoL * qL * (SL - qL) + pL - pR) /
        (rhoR * (SR - qR) - rhoL * (SL - qL))

    return (; SL, SR, SM)
end

function hllc_riemann_solver(
    prob::FDProblem{Grid2D,Euler,<:Any,HLLC},
    rhoL,
    uL,
    vL,
    pL,
    rhoR,
    uR,
    vR,
    pR,
    n,
)
    (; SL, SR, SM) = hllc_wavespeed_estimate(prob, rhoL, uL, vL, pL, rhoR, uR, vR, pR, n)

    if 0 < SL
        FL = flux(prob, rhoL, uL, vL, pL, n)

        F_HLLC = FL
    elseif SL <= 0 < SM
        FL = flux(prob, rhoL, uL, vL, pL, n)

        UL = get_conserved_variables(prob, rhoL, uL, vL, pL)
        EL = UL[4]

        qL = uL * n.x + vL * n.y

        rho_star_L = rhoL * (SL - qL) / (SL - SM)
        p_star = rhoL * (qL - SL) * (qL - SM) + pL
        rho_u_star_L = ((SL - qL) * rhoL * uL + (p_star - pL) * n.x) / (SL - SM)
        rho_v_star_L = ((SL - qL) * rhoL * vL + (p_star - pL) * n.y) / (SL - SM)

        E_star_L = ((SL - qL) * EL - pL * qL + p_star * SM) / (SL - SM)

        UL_star = (rho_star_L, rho_u_star_L, rho_v_star_L, E_star_L)

        F_HLLC = @. FL + SL * (UL_star - UL)
    elseif SM <= 0 <= SR
        FR = flux(prob, rhoR, uR, vR, pR, n)

        UR = get_conserved_variables(prob, rhoR, uR, vR, pR)
        ER = UR[4]

        qR = uR * n.x + vR * n.y

        rho_star_R = rhoR * (SR - qR) / (SR - SM)
        p_star = rhoR * (qR - SR) * (qR - SM) + pR
        rho_u_star_R = ((SR - qR) * rhoR * uR + (p_star - pR) * n.x) / (SR - SM)
        rho_v_star_R = ((SR - qR) * rhoR * vR + (p_star - pR) * n.y) / (SR - SM)

        E_star_R = ((SR - qR) * ER - pR * qR + p_star * SM) / (SR - SM)

        UR_star = (rho_star_R, rho_u_star_R, rho_v_star_R, E_star_R)

        F_HLLC = @. FR + SR * (UR_star - UR)
    elseif SR < 0
        FR = flux(prob, rhoR, uR, vR, pR, n)

        F_HLLC = FR
    else
        throw(error("HLLC Riemann solver failed."))
    end

    return F_HLLC
end

function solve_riemann_problem!(
    prob::FDProblem{Grid2D,Euler,<:Any,HLLC},
    xfluxstore,
    yfluxstore,
    wreconstructed,
)
    (; Nx, Ny) = prob.grid

    ρ = @view wreconstructed[1, :, :, :]
    vx = @view wreconstructed[2, :, :, :]
    vy = @view wreconstructed[3, :, :, :]
    P = @view wreconstructed[4, :, :, :]

    for i in 1:Nx, j in 1:Ny
        ip = i == Nx ? 1 : i + 1
        jp = j == Ny ? 1 : j + 1

        # Solve the Riemann problem on the right cell boundary.
        n = (; x = 1, y = 0)

        rhoL = ρ[i, j, 2]
        uL = vx[i, j, 2]
        vL = vy[i, j, 2]
        pL = P[i, j, 2]

        rhoR = ρ[ip, j, 1]
        uR = vx[ip, j, 1]
        vR = vy[ip, j, 1]
        pR = P[ip, j, 1]

        xfluxstore[:, i, j] .=
            hllc_riemann_solver(prob, rhoL, uL, vL, pL, rhoR, uR, vR, pR, n)

        # Solve the Riemann problem on the top cell boundary.
        n = (; x = 0, y = 1)

        rhoL = ρ[i, j, 4]
        uL = vx[i, j, 4]
        vL = vy[i, j, 4]
        pL = P[i, j, 4]

        rhoR = ρ[i, jp, 3]
        uR = vx[i, jp, 3]
        vR = vy[i, jp, 3]
        pR = P[i, jp, 3]

        yfluxstore[:, i, j] .=
            hllc_riemann_solver(prob, rhoL, uL, vL, pL, rhoR, uR, vR, pR, n)
    end
end