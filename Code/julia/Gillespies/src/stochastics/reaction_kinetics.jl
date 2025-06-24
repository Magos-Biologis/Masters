
"""
Simple function that turns the struct into a term for the drift vector
"""
A_i(R :: ReactionStruct) = R.r * (R.t⁻ - R.t⁺)
function A_i(R :: AbstractArray{T}) where T <: ReactionStruct
    return sum([A_i(r) for r ∈ R])
end


"""
Simple function that turns the struct into a term for the diffusion matrix
"""
Bij(R :: ReactionStruct) = (R.r * R.r') * (R.t⁻ + R.t⁺)
function Bij(R :: AbstractArray{T}) where T <: ReactionStruct
    return sum([Bij(r) for r ∈ R])
end
