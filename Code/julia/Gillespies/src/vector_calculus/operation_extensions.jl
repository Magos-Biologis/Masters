
"""
A julia macro to dynamically construct the needed ∇ operator depending on
use case.
"""
macro ∇(vars...)
    names = [:($(Symbol("D" * "$var"))) for var ∈ vars]
    lines = [:($func = Differential($var)) for (func, var) ∈ zip(names, vars)]
    push!(lines, :( ∇ = reshape(eval.($names), 1, :)))
    return esc(Expr(:block, lines...))
end


"""
Adding the divergence operation to exist under the '·' binary operation
"""
function LinearAlgebra.dot(∇ :: Matrix{<: Differential}, A :: AbstractVector{<: Num})
    return sum([xᵢ|> ∇ᵢ for (xᵢ, ∇ᵢ) ∈ zip(A,∇)])
end

"""
Additionally so for matrix divergence as well
"""
function LinearAlgebra.dot(∇ :: Matrix{<: Differential}, B :: AbstractMatrix{<: Num})
    return [ ∇·A for A ∈ eachcol(B)]
end



"""
Adding the gradiant operation to exist under the '*' binary operation
"""
function Base.:*(∇ :: Matrix{<: Differential}, F :: Num)
    return nothing
end

function Base.:*(∇ :: Matrix{<: Differential}, F :: Num)
    return nothing
end

# Base.:*(∇ :: Matrix{<: Differential}, A :: AbstractMatrix{<: Num}) = sum([xᵢ|> ∇ᵢ for (∇ᵢ, xᵢ) ∈ zip(∇,A)])

# Base.:*(:: Val{:∇}, A :: Vector{<: Num}) = sum([xᵢ|> ∇ᵢ for (∇ᵢ, xᵢ) ∈ zip(∇,A)])
# Base.:*(S :: Symbol, A :: Vector{<: Num}) = *(Val(S), A)
