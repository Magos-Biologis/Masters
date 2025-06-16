
"""
Adding the divergence operation to exist under the '· * ·' binary operation
"""
function Base.:*(∇ :: Matrix{<: Differential}, A :: AbstractVector{<: Num})
    return sum([xᵢ|> ∇ᵢ for (∇ᵢ, xᵢ) ∈ zip(∇, A)])
end

"""
Additionally so for matrix divergence as well
"""
function Base.:*(∇ :: Matrix{<: Differential}, B :: AbstractMatrix{<: Num})
    return collect([sum([xᵢ|> ∇ᵢ for (∇ᵢ, xᵢ) ∈ zip(∇,A)]) for A ∈ eachcol(B)])
end

# Base.:*(∇ :: Matrix{<: Differential}, A :: AbstractMatrix{<: Num}) = sum([xᵢ|> ∇ᵢ for (∇ᵢ, xᵢ) ∈ zip(∇,A)])

# Base.:*(:: Val{:∇}, A :: Vector{<: Num}) = sum([xᵢ|> ∇ᵢ for (∇ᵢ, xᵢ) ∈ zip(∇,A)])
# Base.:*(S :: Symbol, A :: Vector{<: Num}) = *(Val(S), A)
