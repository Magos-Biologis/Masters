
"""
A julia macro to dynamically construct the needed ∇ operator depending on
use case.
"""
macro ∇(vars...)

    # typeof(vars) <: Tuple && (println(vars); println(typeof(vars)))
    # typeof(vars) <: Symbol && (println(vars); println(typeof(vars)))
    #
    # typeof(vars) <: Tuple{Expr} && (vars = vars[1].args) # This is to make the macro more robust
    # typeof(vars) <: Tuple{Symbol} && println(vars) # This is to make the macro more robust
    # return nothing


    names = [:($(Symbol("D" * "$var"))) for var ∈ vars]
    lines = [:($func = Differential($var)) for (func, var) ∈ zip(names, vars)]
    push!(lines, :( ∇ = reshape(eval.($names), 1, :)))
    return esc(Expr(:block, lines...))
end


"""
Adding the gradiant operation to exist under the '*' binary operation
"""
function Base.:*(∇ :: AbstractVecOrMat{<: Differential}, f :: Num)
    return [ f |> ∇ᵢ for ∇ᵢ ∈ ∇ ]
end


"""
Adding the divergence operation to exist under the '·' binary operation
"""
function LinearAlgebra.dot(∇ :: AbstractVecOrMat{<: Differential}, F :: AbstractVector{<: Num})
    return sum([ xᵢ|> ∇ᵢ for (xᵢ, ∇ᵢ) ∈ zip(F,∇) ])
end

"""
Additionally so for matrix divergence as well
"""
function LinearAlgebra.dot(∇ :: AbstractVecOrMat{<: Differential}, A :: AbstractMatrix{<: Num})
    return [ ∇·F for F ∈ eachcol(A) ]
end



# """
# Adding the curl operation to exist under the '×' binary operation
# """
# function LinearAlgebra.cross(∇ :: AbstractVecOrMat{<: Differential}, F :: AbstractVector{<: Num})
#     return [ F |> ∇ᵢ for ∇ᵢ ∈ ∇ ]
# end
