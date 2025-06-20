

"""
Takes in the drift and diffusion elements, alongside with the parameters and
variables returning ∇log(pₛ(x))
"""
function potentiating(L :: T, Par :: N ) where {T <: LangevinType, N <: LangevinParams}
    A, B = L
    vars = Par.variables

    eval(quote @∇ $(vars...) end)
    A_∇B = 2A - expand_derivatives.( ∇·B )

    B⁺ = (vars...) -> langevin_build( B, Par )(vars...) |> pinv
    φ  = langevin_build( A_∇B, Par )

    return (vars...) -> B⁺(vars...) * φ(vars...)
end



function Symbolics.substitute(L :: LangevinType, s :: AbstractDict)
    AA = Num.(substitute(L.A, s))
    BB = Num.(substitute(L.B, s))
    return LangevinType( AA, BB )
end
function Symbolics.substitute(L :: LangevinType, P :: LangevinParams)
    return substitute(L, P.parameters)
end



function langevin_build(L :: T, Par :: N ) where {T <: LangevinType, N <: LangevinParams}
    A, B = L
    AA = langevin_build(A, Par)
    BB = langevin_build(B, Par)
    return LangevinType(AA, BB)
end
function langevin_build( m :: AbstractVecOrMat,  Par :: N ) where N <: LangevinParams
    vars, pars = Par; M = substitute(m, pars)
    MM, M! = build_function(M, vars...; expression = false)
    return MM
end
function langevin_build( m :: Num,  Par :: N ) where N <: LangevinParams
    vars, pars = Par; M = substitute(m, pars)
    MM = build_function(M, vars...; expression = false)
    return MM
end


function gradient_int_two_dim(L :: T, Par :: N;
    ) where {T <: LangevinType, N <: LangevinParams}

    local ε  = eps(Float64)
    local Δx = 0.01

    X = collect( ε:Δx:1-ε )

    f  = potentiating(L, Par)
    ∇F = f.(X, X')

    Fx = ∇F .|> first
    Fy = ∇F .|> last

    Φ = similar(Fx)
    Φ .= 0

    Φ  .= cumsum( Fx .* Δx; dims = 2 )
    Φ .+= cumsum( Fy .* Δx; dims = 1 ) .|> x -> x

    F = exp.(-Φ)

    return F ./ maximum(F)
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




# """
# Needing to define my own Moore-Penrose pseudo-inverse for symbolics as
# it does not already exist in `Symbolics.jl`
# """
# function LinearAlgebra.pinv(A :: AbstractMatrix{Num})
#     m, n = size(A)
#     if m ≠ n
#         m < n && (return (A' *  A)^-1 * A'      )
#         m > n && (return  A' * (A     * A')^-1  )
#     else
#         return A^-1
#     end
# end




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
