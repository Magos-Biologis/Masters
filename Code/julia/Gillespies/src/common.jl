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

function NovelStruct(;
        n :: Real = 10, b :: Real = 0.,
        k⁺:: Vector{<: T} = [1;1], k⁻:: Vector{<: T} = [1;1],
    ) :: NovelStruct where T <: Real
    nn,  bb  = float.(promote(n, b))
    kk⁺, kk⁻ = float.(promote(k⁺, k⁻))

    return NovelStruct{typeof(nn), typeof(kk⁺)}(nn, bb, kk⁺, kk⁻)
end


"""
`S <: DifferentialStruct`
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

    nn, kk, ww, qq = promote(n, k, w, q)
    mm₀ = float(m₀)
    return DifferentialStruct{typeof(mm₀), typeof(nn)}(nn, kk, ww, qq, mm₀)
end


"""
Updating the division method to work with these structs,
importantly, we had to make sure that the datatype of the
non-modified and the modified values are different
"""
Base.:/(::ParameterType, ::Number) = error("Not defined for this type")
Base.:/(P::NovelStruct, N::Number) = NovelStruct(P.n, P.b, /(P.k⁺, N), /(P.k⁻, N))
Base.:/(P::DifferentialStruct, N::Number) = DifferentialStruct(P.n, /(P.k, N), /(P.w, N), /(P.q, N), /(P.m₀, N))

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

struct StepperStruct{ T <: AbstractVector{<: AbstractFloat},
                      S <: AbstractArray{<: Integer} } <: PropensityType
    time   :: T
    states :: S
end

struct PropsAndTrans{F, V <: AbstractVector} <: PropensityType
    propensity :: F
    transition :: V
end

mutable struct SSAOutput{I <: Integer, F <: AbstractFloat} <: PropensityType
    j :: I
    τ :: F
end

SSAOutput() = SSAOutput(0, float(0))

"""
The benefit of now having defined an abstract supertype, is that we can easily
define functions for the base actions like 'iteration'
"""
Base.iterate(S :: T)        where T <: PropensityType = (getfield(S, 1), 1)
Base.iterate(S :: T, state) where T <: PropensityType = state == 1 ? (getfield(S, 2), 2) : nothing




"""
SDE structs go here
"""
abstract type FockPlanType end


"""
Struct for use in defining the langevin equations used by the FPE
We need one for when scalars are involved, and one for when elements of vector
spaces are involved
"""
abstract type LangevinType end

struct Langevin{F₁, F₂} <: LangevinType
    A :: F₁
    B :: F₂
end

struct ScalarLangevin{F} <: LangevinType
    A :: F
    B :: F
end

struct VectorLangevin{F, V <: AbstractVector{F}, M <: AbstractMatrix{F}} <: LangevinType
    A :: V
    B :: M
end

LangevinType(A::AbstractVector, B::AbstractMatrix) = VectorLangevin(A, B)
LangevinType(A::Number, B::Number)                 = ScalarLangevin(A, B)
LangevinType(A::F₁, B::F₂) where {F₁, F₂}          = Langevin(A, B)


"""
It keeps coming up that it is honestly easier to define a new parametric type
then to just account for the possibilities.
"""
const ReactionSpecies{N} = Union{AbstractVector{N}, N} where N <: Num

struct LangevinParams{N <: Num, D <: AbstractDict} <: LangevinType
    variables  :: ReactionSpecies{N}
    parameters :: D
end

function LangevinType(A::T₁, B::T₂, v::ReactionSpecies, d::AbstractDict) where {T₁, T₂}
    return (LangevinType(A, B), LangevinParams(v, d))
end


"""
And again we define the iteration scope
"""
Base.iterate(S :: T)        where T <: LangevinType = (getfield(S, 1), 1)
Base.iterate(S :: T, state) where T <: LangevinType = state == 1 ? (getfield(S, 2), 2) : nothing


"""
"""
abstract type AnalyticalSpace{T, N} end

const OneDimensional{T}   = AnalyticalSpace{T, 1}
const TwoDimensional{T}   = AnalyticalSpace{T, 2}
const ThreeDimensional{T} = AnalyticalSpace{T, 3}


struct OneDimGradient{T} <: OneDimensional{T}
    x :: Vector{T}
    Φ :: Vector{T}
end

struct TwoDimGradient{T} <: TwoDimensional{T}
    x :: Vector{T}
    y :: Vector{T}
    Φ :: Matrix{T}
end

struct ThreeDimGradient{T} <: ThreeDimensional{T}
    x :: Vector{T}
    y :: Vector{T}
    z :: Vector{T}
    Φ :: Array{T}
end




struct AnalyticalSol{T <: BigFloat, V <: AbstractVector{T}} <: OneDimensional{T}
    x :: V
    p :: V
end

# struct AnalyticalSol{V <: AbstractVector{BigFloat}} <: TwoDimensional
#     xs :: V
#     y :: V
# end


"""
So as to more easily control which datatypes are allowed in the ReactionStruct
"""
abstract type ChemicalKinetic end

const AbstractIntSpace{I} = Union{I, AbstractVector{I}} where I <: Integer

struct ReactionStruct{N <: Number, I <: Integer } <: ChemicalKinetic
    t⁺:: N
    t⁻:: N
    r :: AbstractIntSpace{I}
end

function ReactionStruct(; t⁺:: Number = 0, t⁻:: Number = 0, r :: AbstractIntSpace = 1)
    tt⁺, tt⁻ = promote(t⁺, t⁻)
    if typeof(r) <: AbstractVector
        return ReactionStruct{typeof(tt⁺), eltype(r)}(tt⁺, tt⁻, r)
    else
        return ReactionStruct{typeof(tt⁺), typeof(r)}(tt⁺, tt⁻, r)
    end
end



"""
Once again we define the iteration scope, a little more fancy cause its got
more than two
"""
Base.iterate(S :: T)                 where T <: ReactionStruct = (S.t⁺, 1)
function Base.iterate(S :: T, state) where T <: ReactionStruct
    state == 1 && return (S.t⁻, 2)
    state == 2 && return (S.r , 3)
    return nothing
end




#
# struct LangevinStyle <: Base.Broadcast.BroadcastStyle end
#
# BroadcastStyle(::Type{<:LangevinType}) = LangevinStyle()
# BroadcastStyle(::LangevinStyle, ::Broadcast.DefaultArrayStyle) = LangevinStyle()
# BroadcastStyle(::Broadcast.DefaultArrayStyle, ::LangevinStyle) = LangevinStyle()
#
# # Define how to broadcast operations
# function broadcasted(::LangevinStyle, f, ms::Langevin)
#     Langevin(f.(ms.A), f.(ms.B))
# end
#
# function broadcasted(::LangevinStyle, f, ms1::Langevin, ms2::Langevin)
#     Langevin(f.(ms1.A, ms2.A), f.(ms1.B, ms2.B))
# end


# prop::F
# x⃗::Vector{T}
# parameters::Dict{String, Float64}
