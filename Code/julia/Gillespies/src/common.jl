"""
`T <: Parameters` specifies the parameters for stochastic sampling
"""
abstract type ParameterType end

"""
`S <: NovelStruct`
"""
struct NovelStruct{T <: AbstractFloat, V <: AbstractVector{T}} <: ParameterType
    n :: T
    b :: T
    k⁺:: V
    k⁻:: V
end

# function NovelStruct(n::T, b::T, k⁺::V, k⁻::V) where {T <: AbstractFloat, V <: AbstractVector{T}}
#     lengths = Set(length.((k⁺, k⁻)))
#     length(lengths) ≡ 1 ? nothing : error("All vectors must have the same length. Got lengths: $(lengths)")
#     new{T,V}(n, b, k⁺, k⁻)
# end


function NovelStruct(;
        n :: Real = 10, b :: Real = 0.,
        k⁺:: Vector{<: T} = [1;1], k⁻:: Vector{<: T} = [1;1],
    ) :: NovelStruct where T <: Real
    nn, bb   = float.(promote(n, b))
    kk⁺, kk⁻ = float.(promote(k⁺, k⁻))

    return NovelStruct{typeof(nn), typeof(kk⁺)}(nn, bb, kk⁺, kk⁻)
end


"""
`S <: ODE`
"""
struct DifferentialStruct{T <: AbstractFloat, V <: AbstractVector{T}} <: ParameterType
    n :: V
    k :: V
    w :: V
    q :: V
    m₀:: T
end

function DifferentialStruct(;
        n ::Vector{<: T} = [100; 100],     k ::Vector{<: T} = [1.0; 1.0],
        w ::Vector{<: T} = [0.015; 0.015], q ::Vector{<: T} = [0.8; 0.8],
        m₀::Real = 0 ) :: DifferentialStruct where T <: Real

    nn, kk, ww, qq = promote(n, w, k, q)
    mm₀ = float(m₀)
    return DifferentialStruct{typeof(mm₀), typeof(nn)}(nn, kk, ww, qq, mm₀)
end


"""
Updating the division method to work with these structs,
importantly, we had to make sure that the datatype of the
non-modified and the modified values are different
"""
Base.:/(::ParameterType, ::Number) = println("Not defined")
Base.:/(P::NovelStruct, N::Number) = NovelStruct(P.n, P.b, /(P.k⁺, N), /(P.k⁻, N))
Base.:/(P::DifferentialStruct, N::Number) = DifferentialStruct(P.n, /(P.k,N), /(P.w, N), /(P.q, N), /(P.m₀, N))

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

struct StepperStruct{T <: AbstractVector{Real}, S <: AbstractArray{Integer}} <: PropensityType
    time   :: T
    states :: S
end


struct PropsAndTrans{F, V <: AbstractVector}
    propensity :: F
    transition :: V
end


mutable struct SSAOutput{I <: Integer, F <: AbstractFloat} <: PropensityType
    j :: I
    τ :: F
end

SSAOutput() = SSAOutput(0, 0.)


"""
SDE structs go here
"""

abstract type FockPlanType end

abstract type ChemicalDynamic end

"""
Struct for use in defining the langevin equations used by the FPE
"""
struct LangevinStruct{F, V <: AbstractVector{F}, M <: AbstractMatrix{F}} <: FockPlanType
    A :: V
    B :: M
end

struct ReactionStruct{F, V <: AbstractVector} <: ChemicalDynamic
    t⁺:: F
    t⁻:: F
    r :: V
end

function ReactionStruct(t⁺:: Real, t⁻:: Real, r::Vector{<: Integer})

end





# prop::F
# x⃗::Vector{T}
# parameters::Dict{String, Float64}
