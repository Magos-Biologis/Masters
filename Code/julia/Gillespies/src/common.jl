"""
`T <: ParamatersSSA` specifies the parameters for stochastic sampling
"""
abstract type ParameterTypeSSA end

"""
`S <: NovelSSA`
"""
struct NovelStructSSA{T<:Real} <: ParameterTypeSSA
    n ::Real
    b ::Real
    k⁺::Vector{T}
    k⁻::Vector{T}
end

function NovelStructSSA(;
        n ::Real = 10,
        b ::Real = 0,
        k⁺::Vector{T} = [1;1],
        k⁻::Vector{T} = [1;1],
    ) ::NovelStructSSA where T<:Real
    return NovelStructSSA(n, b, k⁺, k⁻)
end

"""
`S <: ODESSA`
"""
struct DifferentialStructSSA{T<:Real,V<:Real} <: ParameterTypeSSA
    k ::Vector{V}
    n ::Vector{T}
    w ::Vector{V}
    q ::Vector{V}
    m₀::Real
end

function DifferentialStructSSA(;
        k ::Vector{V} = [1; 1],
        n ::Vector{T} = [100; 100],
        w ::Vector{V} = [0.015; 0.015],
        q ::Vector{V} = [0.8; 0.8],
        m₀::Real = 0,
    ) ::DifferentialStructSSA where {T<:Real, V<:Real}

    return DifferentialStructSSA(k, n, w, q, m₀)
end


"""
Updating the division method to work with these structs,
importantly, we had to make sure that the datatype of the
non-modified and the modified values are different
"""
Base.:/(::ParameterTypeSSA, ::Number) = println("Not defined")
Base.:/(P::NovelStructSSA, N::Number) = NovelStructSSA(P.n, P.b, /(P.k⁺, N), /(P.k⁻, N))
Base.:/(P::DifferentialStructSSA, N::Number) = DifferentialStructSSA(/(P.k,N), P.n, /(P.w, N), /(P.q, N), /(P.m₀, N))





"""
Attempting to mimic the style of implementation I notice when looking
at the codebase for larger packages like 'StatsBase'.
A lot of abstract types defined before concrete structs are made,
but I am not 100% sure why not all of them define `internal pieces'
so to speak, the method thing, p.A
"""
abstract type PropensityType end

struct StepperSSA{F, T<:Real} <: PropensityType
    prop::F
    x⃗::Vector{T}
    parameters::Dict{String, Float64}
end
