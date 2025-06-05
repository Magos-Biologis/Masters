


"""
`F <: ParamatersSSA` specifies the parameters for stochastic sampling
"""
abstract type ParametersSSA end

struct NovelSSA{T<:Real} <: ParametersSSA
    k⁺::Vector{T}
    k⁻::Vector{T}
end


struct ODESSA{T<:Real} <: ParametersSSA
    k::Vector{T}
    n::Vector{T}
    w::Vector{T}
    q::Vector{T}
end


## Parameter division
# div(P::SSAParameters, s::Int)::SSAParameters = SSAParameters()
div(P::ODESSA, s::Real)::ODESSA
div(P::NovelSSA, s::Real)::NOVELSSA


