abstract type RiemannSolver end

struct NaiveRS <: RiemannSolver end
struct HLLC <: RiemannSolver end
struct Roe <: RiemannSolver end

function solve_riemann(method::NaiveRS, gd, w_reconstruct)
    Fl = zeros(3, gd.Nx)

    for i = 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1

        Fminus =
            flux(w_reconstruct[1, i, 2], w_reconstruct[2, i, 2], w_reconstruct[3, i, 2])
        Fplus =
            flux(w_reconstruct[1, ip, 1], w_reconstruct[2, ip, 1], w_reconstruct[3, ip, 1])

        @. Fl[:, i] = 1 / 2 * (Fplus + Fminus)
        # @. Fl[:, i] = Fminus
    end

    return Fl
end

function hllc_wavespeed_estimate(rhoL, vL, pL, rhoR, vR, pR)
    γ = 5 / 3

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

function hllc_riemann_solver(rhoL, vL, pL, rhoR, vR, pR)
    (; SL, SR, SM) = hllc_wavespeed_estimate(rhoL, vL, pL, rhoR, vR, pR)
    γ = 5 / 3

    if 0 < SL
        FL = flux(rhoL, vL, pL)
        return FL
    elseif SL <= 0 < SM
        FL = flux(rhoL, vL, pL)
        EL = pL / (γ - 1) + 1 / 2 * rhoL * vL^2

        UL = [rhoL, rhoL * vL, EL]

        rhoLstar = rhoL * (SL - vL) / (SL - SM)
        pLstar = pL + rhoL * (vL - SL) * (vL - SM)
        rhoLvLstar = (rhoL * vL * (SL - vL) + pLstar - pL) / (SL - SM)

        ELstar = (EL * (SL - vL) - pL * vL + pLstar * SM) / (SL - SM)

        ULstar = [rhoLstar, rhoLvLstar, ELstar]

        FLstar = FL + SL * (ULstar - UL)
        return FLstar
    elseif SM <= 0 <= SR
        FR = flux(rhoR, vR, pR)
        ER = pR / (γ - 1) + 1 / 2 * rhoR * vR^2

        UR = [rhoR, rhoR * vR, ER]

        rhoRstar = rhoR * (SR - vR) / (SR - SM)
        pRstar = pR + rhoR * (vR - SR) * (vR - SM)
        rhoRvRstar = (rhoR * vR * (SR - vR) + pRstar - pR) / (SR - SM)

        ERstar = (ER * (SR - vR) - pR * vR + pRstar * SM) / (SR - SM)

        URstar = [rhoRstar, rhoRvRstar, ERstar]

        FRstar = FR + SR * (URstar - UR)
        return FRstar
    elseif SR < 0
        FR = flux(rhoR, vR, pR)
        return FR
    end
end

function solve_riemann(method::HLLC, gd, w_reconstruct)
    Fl = zeros(3, gd.Nx)

    for i = 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1

        rhoL = w_reconstruct[1, i, 2]
        vL = w_reconstruct[2, i, 2]
        pL = w_reconstruct[3, i, 2]

        rhoR = w_reconstruct[1, ip, 1]
        vR = w_reconstruct[2, ip, 1]
        pR = w_reconstruct[3, ip, 1]

        F = hllc_riemann_solver(rhoL, vL, pL, rhoR, vR, pR)

        @. Fl[:, i] = F
    end

    return Fl
end