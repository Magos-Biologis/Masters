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

function NovelStructSSA(n, b, k⁺, k⁻)
    T = promote_type(typeof(n), typeof(b), eltype(k⁺), eltype(k⁻))
    return NovelStructSSA(T(n), T(b), T.(k⁺), T.(k⁻))
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
struct DifferentialStructSSA{T<:Real} <: ParameterTypeSSA
    k ::Vector{T}
    n ::Vector{T}
    w ::Vector{T}
    q ::Vector{T}
    m₀::T
end

function DifferentialStructSSA(k, n, w, q, m₀)
    T = promote_type(eltype(k), eltype(n), eltype(w), eltype(q), typeof(m₀))
    return DifferentialStructSSA(T.(k), T.(n), T.(w), T.(q), T(m₀))
end

function DifferentialStructSSA(;
        k ::Vector{<:Real} = [1.0; 1.0],
        n ::Vector{<:Real} = [100; 100],
        w ::Vector{<:Real} = [0.015; 0.015],
        q ::Vector{<:Real} = [0.8; 0.8],
        m₀::Real = 0.,
    ) ::DifferentialStructSSA

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
Also overloading the convert type to work for my structs
"""




"""
Attempting to mimic the style of implementation I notice when looking
at the codebase for larger packages like 'StatsBase'.
A lot of abstract types defined before concrete structs are made,
but I am not 100% sure why not all of them define `internal pieces'
so to speak, the method thing, p.A
"""
abstract type PropensityType end

struct StepperStruct{T<:Real, S<:Integer} <: PropensityType
    time::Vector{T}
    states::AbstractArray{S}
end


# prop::F
# x⃗::Vector{T}
# parameters::Dict{String, Float64}
