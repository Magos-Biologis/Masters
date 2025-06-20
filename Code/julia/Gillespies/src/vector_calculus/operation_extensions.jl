
"""
A julia macro to dynamically construct the needed ∇ operator depending on
use case.
"""
macro ∇(vars...)
    set_var = :( @variables $(vars...) )
    names = [:($(Symbol("D" * "$var"))) for var ∈ vars]
    dlines = [:($func = Differential($var)) for (func, var) ∈ zip(names, vars)]
    push!(dlines, :( ∇ = reshape(eval.($names), 1, :)))
    return esc(Expr(:block, set_var, dlines...))
end


"""
Adding the gradiant operation to exist under the '*' binary operation
"""
function Base.:*(∇ :: AbstractVecOrMat{<: Differential}, f :: Num)
    return [ f |> ∇ᵢ for ∇ᵢ ∈ ∇ ]
end


"""
Adding the divergence operation to exist under the '·' binary operation
```julia
```
"""
function LinearAlgebra.dot(∇ :: AbstractVecOrMat{<: Differential}, F :: AbstractVector{<: Num})
    return sum([ xᵢ|> ∇ᵢ for (xᵢ, ∇ᵢ) ∈ zip(F,∇) ])
end
function LinearAlgebra.dot(∇ :: Differential, F :: Num)
    return F |> ∇
end


"""
Additionally so for matrix divergence as well
"""
function LinearAlgebra.dot(∇ :: AbstractVecOrMat{<: Differential}, A :: AbstractMatrix{<: Num})
    return [ ∇·F for F ∈ eachcol(A) ]
end



"""
Defining the Laplacian
"""
function LinearAlgebra.dot(∇₁ :: AbstractVector{D}, ∇₂ :: AbstractVector{D}) where D <: Differential
    return sum([ a |> b for (a, b) ∈ zip(∇₁,∇₂) ])
end

function Base.:^(∇ :: AbstractVecOrMat{<: Differential}, n :: Integer)
    return fill(∇, n) |> splat(·)
end









# """
# Turns out I can overload macros too
# """
# macro ∇(list_of_vars)
#     vars = esc(list_of_vars)
#     parse = :( @∇  $(vars...)  )
#     return esc(parse)
# end







# """
# Adding the curl operation to exist under the '×' binary operation
# """
# function LinearAlgebra.cross(∇ :: AbstractVecOrMat{<: Differential}, F :: AbstractVector{<: Num})
#     return [ F |> ∇ᵢ for ∇ᵢ ∈ ∇ ]
# end
