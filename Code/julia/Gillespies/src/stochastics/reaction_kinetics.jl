
"""
Simple function that turns the struct into a term for the drift vector
"""
function A_i(R :: ReactionStruct)
    t⁺, t⁻, r = R
    return r * (t⁻ - t⁺)
end
function A_i(R :: AbstractArray{T}) where T <: ReactionStruct
    return sum([A_i(r) for r ∈ R])
end


"""
Simple function that turns the struct into a term for the diffusion matrix
"""
function Bij(R :: ReactionStruct)
    t⁺, t⁻, r = R
    return (r * r') * (t⁻ + t⁺)
end
function Bij(R :: AbstractArray{T}) where T <: ReactionStruct
    return sum([Bij(r) for r ∈ R])
end
