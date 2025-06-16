
"""
Simple function that turns the struct into a term for the drift vector
"""
function A_i(R :: ReactionStruct)
    t⁺, t⁻, r = R
    return r * (t⁻ - t⁺)
end


"""
Simple function that turns the struct into a term for the diffusion matrix
"""
function B_ij(R :: ReactionStruct)
    t⁺, t⁻, r = R
    return (r * r') * (t⁻ - t⁺)
end
