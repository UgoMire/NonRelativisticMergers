"""
    solve_poisson(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        ρ
)

FFT Poisson equation solver in 1D. This implementation is not optimized but
maybe clearer to understand.
"""
function solve_poisson(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        ρ;
        fft_type = :fft
)
    (; grid, model) = prob
    (; Nx, Δx) = grid
    (; G, ϵ) = model

    if fft_type == :fft
        kx = 2π * fftfreq(Nx, 1 / Δx)

        ρhat = fft(ρ)

        uhat = zeros(Complex{Float64}, Nx)
    elseif fft_type == :rfft
        kx = 2π * rfftfreq(Nx, 1 / Δx)

        ρhat = rfft(ρ)

        uhat = zeros(Complex{Float64}, Nx ÷ 2 + 1)
    end

    for (i, kxi) in enumerate(kx)
        if kxi == 0
            uhat[i] = 0
        else
            uhat[i] = -4π * G * ρhat[i] / (kxi^2 + ϵ^2)
        end
    end

    if fft_type == :fft
        return real(ifft(uhat))
    elseif fft_type == :rfft
        return real(irfft(uhat, Nx))
    end
end

"""
    solve_poisson(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        ρ
)

FFT Poisson equation solver in 2D. This implementation is not optimized but
maybe clearer to understand.
"""
function solve_poisson(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        ρ;
        fft_type = :rfft
)
    (; grid, model) = prob
    (; Nx, Ny, Δx, Δy) = grid
    (; G, ϵ) = model

    if fft_type == :fft
        kx = 2π * fftfreq(Nx, 1 / Δx)
        ky = 2π * fftfreq(Ny, 1 / Δy)

        ρhat = fft(ρ)

        uhat = zeros(Complex{Float64}, Nx, Ny)
    elseif fft_type == :rfft
        kx = 2π * rfftfreq(Nx, 1 / Δx)
        ky = 2π * rfftfreq(Ny, 1 / Δy)

        ρhat = rfft(ρ)

        uhat = zeros(Complex{Float64}, Nx ÷ 2 + 1, Ny)
    end

    for (i, kxi) in enumerate(kx), (j, kyj) in enumerate(ky)
        if kxi == 0 && kyj == 0
            uhat[i, j] = 0
        else
            uhat[i, j] = -4π * G * ρhat[i, j] / (kxi^2 + kyj^2 + ϵ^2)
        end
    end

    if fft_type == :fft
        return real(ifft(uhat))
    elseif fft_type == :rfft
        return real(irfft(uhat, Nx))
    end
end

function setup_fft_cache(grid::Grid1D)
    (; Nx, Δx) = grid

    return (;
        kx = 2π * rfftfreq(Nx, 1 / Δx),
        ρ̂_cache = zeros(ComplexF64, Nx ÷ 2 + 1),
        û_cache = zeros(ComplexF64, Nx ÷ 2 + 1),
        planned_fft = plan_rfft(zeros(Nx))
    )
end

"""
    solve_poisson!(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        cache
)

FFT Poisson equation solver in 1D. This implementation is do not allocate thanks
to the cache.
"""
function solve_poisson!(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        (; kx, planned_fft, ρ̂_cache, û_cache)
)
    (; model) = prob
    (; G, ϵ) = model

    mul!(ρ̂_cache, planned_fft, ρ)

    for (i, kxi) in enumerate(kx)
        if kxi == 0
            û_cache[i] = 0
        else
            û_cache[i] = -4π * G * ρ̂_cache[i] / (kxi^2 + ϵ^2)
        end
    end

    mul!(potentialstore, inv(planned_fft), û_cache)
end

function setup_fft_cache(grid::Grid2D)
    (; Nx, Δx, Ny, Δy) = grid

    return (;
        kx = 2π * fftfreq(Nx, 1 / Δx),
        ky = 2π * fftfreq(Ny, 1 / Δy),
        ρ_complex = zeros(ComplexF64, Nx, Ny),
        potentialstore_complex = zeros(ComplexF64, Nx, Ny),
        ρhat_cache = zeros(ComplexF64, Nx, Ny),
        uhat_cache = zeros(ComplexF64, Nx, Ny),
        planned_fft = plan_fft(zeros(Nx, Ny)),
        planned_ifft = plan_ifft(zeros(Nx, Ny))
    )
end

"""
    solve_poisson!(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        cache
)

FFT Poisson equation solver in 1D. This implementation is do not allocate thanks
to the cache.
"""
function solve_poisson!(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        (;
            kx,
            ky,
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

    for (i, kxi) in enumerate(kx), (j, kyj) in enumerate(ky)
        if kxi == 0 && kyj == 0
            uhat_cache[i, j] = 0
        else
            uhat_cache[i, j] = -4π * G * ρhat_cache[i, j] / (kxi^2 + kyj^2 + ϵ^2)
        end
    end

    mul!(potentialstore_complex, planned_ifft, uhat_cache)

    @. potentialstore = real(potentialstore_complex)
end