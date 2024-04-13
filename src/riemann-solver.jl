abstract type RiemannSolver end

struct NaiveRS <: RiemannSolver end
struct HLLC <: RiemannSolver end
struct Roe <: RiemannSolver end

function solve_riemann_problem!(fluxstore, method::NaiveRS, model, gd, wstore)
    for j = 1:3, i = 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1

        Fminux = flux(model, j, wstore[1, i, 2], wstore[2, i, 2], wstore[3, i, 2])
        Fplus = flux(model, j, wstore[1, ip, 1], wstore[2, ip, 1], wstore[3, ip, 1])

        fluxstore[j, i] = 1 / 2 * (Fminux + Fplus)
    end
end

function hllc_wavespeed_estimate(model, rhoL, vL, pL, rhoR, vR, pR)
    γ = model.γ

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

function hllc_riemann_solver(model, rhoL, vL, pL, rhoR, vR, pR)
    (; SL, SR, SM) = hllc_wavespeed_estimate(model, rhoL, vL, pL, rhoR, vR, pR)
    γ = model.γ

    if 0 < SL
        FL = flux(model, rhoL, vL, pL)
        return FL
    elseif SL <= 0 < SM
        FL = flux(model, rhoL, vL, pL)
        EL = pL / (γ - 1) + 1 / 2 * rhoL * vL^2

        UL = (rhoL, rhoL * vL, EL)

        rhoLstar = rhoL * (SL - vL) / (SL - SM)
        pLstar = pL + rhoL * (vL - SL) * (vL - SM)
        rhoLvLstar = (rhoL * vL * (SL - vL) + pLstar - pL) / (SL - SM)

        ELstar = (EL * (SL - vL) - pL * vL + pLstar * SM) / (SL - SM)

        ULstar = (rhoLstar, rhoLvLstar, ELstar)

        return @. FL + SL * (ULstar - UL)
    elseif SM <= 0 <= SR
        FR = flux(model, rhoR, vR, pR)
        ER = pR / (γ - 1) + 1 / 2 * rhoR * vR^2

        UR = (rhoR, rhoR * vR, ER)

        rhoRstar = rhoR * (SR - vR) / (SR - SM)
        pRstar = pR + rhoR * (vR - SR) * (vR - SM)
        rhoRvRstar = (rhoR * vR * (SR - vR) + pRstar - pR) / (SR - SM)

        ERstar = (ER * (SR - vR) - pR * vR + pRstar * SM) / (SR - SM)

        URstar = (rhoRstar, rhoRvRstar, ERstar)

        return @. FR + SR * (URstar - UR)
        # elseif SR < 0
    else
        FR = flux(model, rhoR, vR, pR)
        return FR
    end
end

function solve_riemann_problem!(fluxstore, wreconstructed, method::HLLC, model, gd)

    for i = 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1

        rhoL = wreconstructed[1, i, 2]
        vL = wreconstructed[2, i, 2]
        pL = wreconstructed[3, i, 2]

        rhoR = wreconstructed[1, ip, 1]
        vR = wreconstructed[2, ip, 1]
        pR = wreconstructed[3, ip, 1]

        fluxstore[:, i] .= hllc_riemann_solver(model, rhoL, vL, pL, rhoR, vR, pR)
    end
end