"""
`F <: ParamatersSSA` specifies the parameters for stochastic sampling
"""
abstract type ParameterStructSSA end

"""
`F <: NovelSSA`
"""
struct NovelSSA{T<:Real,V<:Real} <: ParameterStructSSA
    n::T
    b::T
    k⁺::Vector{V}
    k⁻::Vector{V}
end


"""
`F <: ODESSA`
"""
struct ODESSA{T<:Real,V<:Real} <: ParameterStructSSA
    k::Vector{V}
    n::Vector{T}
    w::Vector{V}
    q::Vector{V}
    m₀::V
end


"""
Updating the division method to work with these structs,
importantly, we had to make sure that the datatype of the
non-modified and the modified values are different
"""
Base.:/(::ParameterStructSSA, ::Number) = println("Not defined")
Base.:/(P::NovelSSA, N::Number) = NovelSSA(P.n, P.b, /(P.k⁺, N), /(P.k⁻, N))
Base.:/(P::ODESSA, N::Number) = ODESSA(/(P.k,N), P.n, /(P.w, N), /(P.q, N), /(P.m₀, N))


