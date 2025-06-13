
"""
Structs go here
"""

abstract type FockPlanType end

abstract type ChemicalDynamic end

# struct SimpleChemicalStruct <: FockPlanType end


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


"""
`T <: ParamatersTypeCalc` specifies the parameters for stochastic sampling
"""
abstract type ParameterTypeCalc end

"""
`S <: NovelStructCalc`
"""
struct NovelStructCalc{T <: AbstractFloat, V <: AbstractVector{T}} <: ParameterTypeCalc
    n :: T
    b :: T
    k⁺:: V
    k⁻:: V

    function NovelStructCalc(n::T, b::T, k⁺:: V, k⁻:: V)
        lengths = [length(k⁺), length(k⁻)]
        if length(unique(lengths)) != 1
            error("All vectors must have the same length. Got lengths: $(lengths)")
        end
        new(n, b, k⁺, k⁻)
    end
end

function NovelStructCalc(;
        n :: Real = 10, b :: Real = 0,
        k⁺:: Vector{<: T} = [1;1], k⁻:: Vector{<: T} = [1;1],
    ) :: NovelStructCalc where T <: Real
    nn, bb   = float.(promote(n, b))
    kk⁺, kk⁻ = float.(promote(k⁺, k⁻))

    return NovelStructCalc{typeof(nn), typeof(kk⁺)}(nn, bb, kk⁺, kk⁻)
end
# return NovelStructCalc{Float64, Vector{Float64}}(n, b, k⁺, k⁻)

"""
`S <: DifferentialStructCalc`
"""
struct DifferentialStructCalc{T <: AbstractFloat, V <: AbstractVector{T}} <: ParameterTypeCalc
    n :: V
    k :: V
    w :: V
    q :: V
    m₀:: T
end

function DifferentialStructCalc(;
        n ::Vector{<: T} = [100; 100], k ::Vector{<: T} = [1.0; 1.0],
        w ::Vector{<: T} = [0.015; 0.015], q ::Vector{<: T} = [0.8; 0.8],
        m₀::Real = 0 ) :: DifferentialStructCalc where T <: Real

    nn, kk, ww, qq = promote(n, w, k, q)
    mm₀ = float(m₀)
    return DifferentialStructCalc{typeof(mm₀), typeof(nn)}(nn, kk, ww, qq, mm₀)
end

