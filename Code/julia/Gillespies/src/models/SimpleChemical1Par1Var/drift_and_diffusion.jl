"""
The Drift vector and Diffusion matrix for the simple chemical system
"""
function SimpleChemical1Par1VarAB(P :: NovelStruct; symbolic = true)

    x = variable(:x; T=Real)

    vars = x
    params = Dict{Union{Symbol, Num}, Number}(
                  :k₁  => P.k⁺[1], :k₋₁ => P.k⁻[1],
                  :n   => P.n
                 )

    @parameterification params


    r₁ = ReactionStruct(t⁺ = k₁ * x, t⁻ = k₋₁ * (1 - x), r = 1)
    A  = A_i(r₁)
    B  = Bij(r₁) / n


    if symbolic
        return LangevinType(A, B, [vars], params)
    else
        AA = build_function(substitute(A, params), vars; expression = Val{false})
        BB = build_function(substitute(B, params), vars; expression = Val{false})

        return LangevinType(AA, BB)
    end
end
