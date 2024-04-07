struct Grid1D
    xmin::Float64
    xmax::Float64
    Nx::Int64
    Δx::Float64
    xl::Vector{Float64}

    function Grid1D(; xmin=0, xmax=1, Nx=100)
        xl = range(xmin, xmax, Nx)
        Δx = xl[2] - xl[1]
        new(xmin, xmax, Nx, Δx, xl)
    end
end