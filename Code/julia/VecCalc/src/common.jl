
"""
Structs go here
"""

abstract type FockPlanType end

# struct SimpleChemicalStruct <: FockPlanType end


"""
Struct for use in defining the langevin equations used by the FPE
"""
struct LangevinStruct{F} <: FockPlanType
    A :: AbstractVector{F}
    B :: AbstractMatrix{F}
end


function LangevinEquation(Lₑ:: LangevinStruct)

    dW = @brownian a b c d
    println(dW)

    return Lₑ.A + Lₑ.B * dW

end


struct ReactionStruct{F} <: FockPlanType
    t⁺:: F
    t⁻:: F
    r :: AbstractVector{Integer}
end


# function LangevinStruct(A, B)
#     A = convert(AbstractVector{Number}, A)
#     B = convert(AbstractMatrix{Number}, B)
#     T = promote_type(eltype(A), eltype(B))
#     return LangevinStruct(T.(A), T.(B))
# end
#
# function LangevinStruct(;
#         A :: AbstractVector{T} = [0],
#         B :: AbstractMatrix{T} = [0]
#     ) :: LangevinStruct where T <: Real
#     return LangevinStruct(A, B)
# end


"""
`T <: ParamatersSSA` specifies the parameters for stochastic sampling
"""
abstract type ParameterTypeCalc end

"""
`S <: NovelStructCalc`
"""
struct NovelStructCalc{T <: AbstractFloat} <: ParameterTypeCalc
    n :: T
    b :: T
    k⁺:: AbstractVector{T}
    k⁻:: AbstractVector{T}
end




# n = convert(AbstractFloat, n)
# b = convert(AbstractFloat, b)
# k⁺= convert(AbstractVector{AbstractFloat}, k⁺)
# k⁻= convert(AbstractVector{AbstractFloat}, k⁻)

# function NovelStructCalc(n, b, k⁺, k⁻)
#     n = float(n)
#     b = float(b)
#     k⁺= float.(k⁺)
#     k⁻= float.(k⁻)
#
#     T = promote_type(typeof(n), typeof(b), eltype(k⁺), eltype(k⁻))
#     return NovelStructCalc(T(n), T(b), T.(k⁺), T.(k⁻))
# end

function NovelStructCalc(;
        n :: Number = 10,
        b :: Number = 0,
        k⁺:: Vector{T} = [1; 1],
        k⁻:: Vector{T} = [1; 1],
    ) :: NovelStructCalc where T <: Number
    return NovelStructCalc{AbstractFloat}(n, b, k⁺, k⁻)
end

"""
`S <: DifferentialStructCalc`
"""
struct DifferentialStructCalc{T <: AbstractFloat} <: ParameterTypeCalc
    n :: Vector{T}
    k :: Vector{T}
    w :: Vector{T}
    q :: Vector{T}
    m₀:: T
end

function DifferentialStructCalc(;
        n :: Vector{Real} = [100; 100],
        k :: Vector{Real} = [1.0; 1.0],
        w :: Vector{Real} = [0.015; 0.015],
        q :: Vector{Real} = [0.8; 0.8],
        m₀:: Real = 0,
    ) :: DifferentialStructCalc
    return DifferentialStructCalc{AbstractFloat}(n, k, w, q, m₀)
end


# # n = convert(AbstractVector{AbstractFloat}, n)
# # k = convert(AbstractVector{AbstractFloat}, k)
# # w = convert(AbstractVector{AbstractFloat}, w)
# # q = convert(AbstractVector{AbstractFloat}, q)
# # m₀= convert(AbstractFloat, m₀)
# function DifferentialStructCalc(n, k, w, q, m₀)
#     n = float.(n)
#     k = float.(k)
#     w = float.(w)
#     q = float.(q)
#     m₀= float(m₀)
#     T = promote_type(eltype(n), eltype(k), eltype(w), eltype(q), typeof(m₀))
#     return DifferentialStructCalc(T.(n), T.(k), T.(w), T.(q), T(m₀))
# end



