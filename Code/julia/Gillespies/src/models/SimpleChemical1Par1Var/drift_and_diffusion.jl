"""
The Drift vector and Diffusion matrix for the simple chemical system
"""
function SimpleChemical1Par1VarAB(P :: NovelStruct; symbolic = true)

    @variables x::Real

    vars = [x]
    params :: Dict{Union{Symbol, Num}, Number} = Dict(
                  :k₁  => P.k⁺[1],
                  :k₋₁ => P.k⁻[1],
                  :n   => P.n
                 )
    @parameterification params


    r₁ = ReactionStruct(t⁺ = k₁ * x, t⁻ = k₋₁ * (1 - x), r = 1)
    A  = A_i(r₁)
    B  = Bij(r₁) / n

    # A =  r₁        .* (t⁻ - t⁺)
    # B = (r₁ * r₁') .* (t⁻ + t⁺) / n

    if symbolic
        return LangevinType(A, B), vars, params
    else
        AA = build_function(substitute(A, params), x; expression = Val{false})
        BB = build_function(substitute(B, params), x; expression = Val{false})

        return LangevinType(AA, BB)
    end
end
