


function potentiating(L :: T, vars :: Num...) where {T <: LangevinType}

    A, B = L
    eval(quote @∇ $(vars...) end)

    ∇_B = expand_derivatives.(∇ ⋅B)
    ∇log_p = inv(B) * (A - ∇_B)
end


# abstract type AbstractShapedMatrix{T,M,N} <: AbstractArray{T,2} end
# abstract type AbstractShapedVector{T,M  } <: AbstractArray{T,1} end
#
# struct ShapedMatrix{T,M,N} #<: AbstractShapedMatrix{M,N,T}
#     array::AbstractArray{T,2}
#
#     function ShapedMatrix(array::AbstractArray{T,2}) where T
#         M, N = size(array)
#         new{T,M,N}(array)
#     end
# end
#
# const SquareMatrix{T,M} = ShapedMatrix{T,M,M}
#
#
# const   ColumnVector{T,M} = ShapedMatrix{T,M,1}
# const      RowVector{T,N} = ShapedMatrix{T,1,N}
#  This is high-key getting too ambitious for the time being, I need to
#  reign this the fuck in
#
#  Why do I not listen to myself, dear fucking lord



"""
Needing to define my own Moore-Penrose pseudo-inverse for symbolics as
it does not already exist in `Symbolics.jl`
"""
function LinearAlgebra.pinv(A :: AbstractMatrix{Num})
    m, n = size(A)
    if m ≠ n
        m < n && (return (A' *  A)^-1 * A'      )
        m > n && (return  A' * (A     * A')^-1  )
    else
        return A^-1
    end
end




# function divergence(A :: M) where M <: AbstractMatrix
#     m, n = size(A)
#     Dx = Differential(x)
#     Dy = Differential(y)
#     Dz = Differential(z)
#
#     ∇ = [Dx Dy Dz][1:m]
#
#     return ∇ * A
# end
#


# The `pinv()` of LinearAlgebra already does this
#
# function psuedo_inverse(B :: AbstractMatrix)
#     m, n = size(B)
#
#     if m ≡ n
#         (Q, R = qr(B))
#     elseif min(m, n) ≡ m
#         B⁺= (B' * B)^-1 * B'
#     elseif min(m, n) ≡ n
#         B⁺= B' * (B * B')^-1
#     end
#
# end
