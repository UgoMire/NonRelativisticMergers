function solve_poisson!(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ
)
    (; grid, model) = prob
    (; Nx, Δx) = grid
    (; G) = model

    kx = 2π * fftfreq(Nx, 1 / Δx)

    ρhat = fft(ρ)

    uhat = zeros(Complex{Float64}, Nx)

    ϵ = 10
    for (i, kxi) in enumerate(kx)
        if kxi == 0
            uhat[i] = 0
        else
            uhat[i] = -4π * G * ρhat[i] / (kxi^2 + ϵ^2)
        end
    end

    potentialstore .= real(ifft(uhat))
end

function solve_poisson(prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any}, ρ)
    (; grid) = prob
    (; Nx) = grid

    potentialstore = zeros(Nx)

    solve_poisson!(prob, potentialstore, ρ)

    return potentialstore
end

function setup_fft_cache(grid::Grid1D)
    (; Nx, Δx) = grid

    return (;
        kx = 2π * fftfreq(Nx, 1 / Δx),
        ρ_complex = zeros(ComplexF64, Nx),
        potentialstore_complex = zeros(ComplexF64, Nx),
        ρhat_cache = zeros(ComplexF64, Nx),
        uhat_cache = zeros(ComplexF64, Nx),
        planned_fft = plan_fft(zeros(Nx)),
        planned_ifft = plan_ifft(zeros(Nx))
    )
end

function solve_poisson!(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        (;
            kx,
            planned_fft,
            planned_ifft,
            ρhat_cache,
            uhat_cache,
            ρ_complex,
            potentialstore_complex
        )
)
    (; model) = prob
    (; G, ϵ) = model

    ρ_complex .= ρ
    potentialstore_complex .= potentialstore

    mul!(ρhat_cache, planned_fft, ρ_complex)

    for (i, kxi) in enumerate(kx)
        if kxi == 0
            uhat_cache[i] = 0
        else
            uhat_cache[i] = -4π * G * ρhat_cache[i] / (kxi^2 + ϵ^2)
        end
    end

    mul!(potentialstore_complex, planned_ifft, uhat_cache)

    @. potentialstore = real(potentialstore_complex)
end
