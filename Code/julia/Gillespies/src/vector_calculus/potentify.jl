
function potentiating(L :: T, vars :: Num...) where {T <: LangevinType}

    A, B = L
    eval(quote @∇ $(vars...) end)

    ∇_B = expand_derivatives.(∇ ⋅B)
    ∇log_p = inv(B) * (A - ∇_B)
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
