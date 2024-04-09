abstract type RiemannSolver end

struct NaiveRS <: RiemannSolver end
struct HLLC <: RiemannSolver end
struct Roe <: RiemannSolver end

function solve_riemann(method::NaiveRS, gd, w_reconstruct)
    Fl = zeros(3, gd.Nx)

    for i = 1:gd.Nx
        Fminus =
            flux(w_reconstruct[1, i, 1], w_reconstruct[2, i, 1], w_reconstruct[3, i, 1])
        Fplus = flux(w_reconstruct[1, i, 2], w_reconstruct[2, i, 2], w_reconstruct[3, i, 2])

        @. Fl[:, i] = 1 / 2 * (Fplus + Fminus)
    end

    return Fl
end

function solve_riemann(method::HLLC, gd, w_reconstruct)
    Fl = zeros(3, gd.Nx)

    for i = 1:gd.Nx
        ρL = w_reconstruct[1, i, 1]
        ρR = w_reconstruct[1, i, 2]
        vL = w_reconstruct[2, i, 1]
        vR = w_reconstruct[2, i, 2]
        pL = w_reconstruct[3, i, 1]
        pR = w_reconstruct[3, i, 2]

        cL = sqrt(5 / 3 * pL / ρL)
        cR = sqrt(5 / 3 * pR / ρR)

        HL = pL / ρL * 5 / 3 / (5 / 3 - 1) + 1 / 2 * vL^2
        HR = pR / ρR * 5 / 3 / (5 / 3 - 1) + 1 / 2 * vR^2

        Rρ = sqrt(ρR / ρL)

        vtilde = (vL + vR * Rρ) / (1 + Rρ)
        Htilde = (HL + HR * Rρ) / (1 + Rρ)
        ctilde = sqrt((5 / 3 - 1) * (Htilde - 1 / 2 * vtilde^2))

        # β = sqrt((5 / 3 - 1) / (2 * 5 / 3))
        β = 1

        SL = min(vL - β * cL, vtilde - ctilde)
        SR = max(vR + β * cR, vtilde + ctilde)
        SM =
            (ρR * vR * (SR - vR) - ρL * vL * (SL - vL) + pL - pR) /
            (ρR * (SR - vR) - ρL * (SL - vL))

        if 0 < SL
            @show "1"
            FL = flux(ρL, vL, pL)
            @. Fl[:, i] = FL
        elseif SL <= 0 < SM
            FL = flux(ρL, vL, pL)
            EL = pL / (5 / 3 - 1) + 1 / 2 * ρL * vL^2

            UL = [ρL, ρL * SL, EL]

            ρLstar = ρL * (SL - vL) / (SL - SM)
            pLstar = pL + ρL * (vL - SL) * (vL - SM)
            ρLvLstar = (ρL * vL * (SL - vL) + pLstar - pL) / (SL - SM)

            ELstar = (EL * (SL - vL) - pL * vL + pLstar * SM) / (SL - SM)

            ULstar = [ρLstar, ρLvLstar, ELstar]

            @. Fl[:, i] = FL + SL * (ULstar - UL)
        elseif SM <= 0 <= SR
            FR = flux(ρR, vR, pR)
            ER = pR / (5 / 3 - 1) + 1 / 2 * ρR * vR^2

            UR = [ρR, ρR * SR, ER]

            ρRstar = ρR * (SR - vR) / (SR - SM)
            pRstar = pR + ρR * (vR - SR) * (vR - SM)
            ρRvRstar = (ρR * vR * (SR - vR) + pRstar - pR) / (SR - SM)

            ERstar = (ER * (SR - vR) - pR * vR + pRstar * SM) / (SR - SM)

            URstar = [ρRstar, ρRvRstar, ERstar]

            @. Fl[:, i] = FR + SR * (URstar - UR)
        elseif SR < 0
            FR = flux(ρR, vR, pR)
            @. Fl[:, i] = SR
        end
    end

    return Fl
end