"""
The Drift vector and Diffusion matrix for the simple chemical system
"""
function SimpleChemical1Par1VarAB(P :: NovelStruct; symbolic = true)

    @variables x::Real
    params = Dict{Union{Symbol, Num}, Any}(
                  :k₁ => P.k⁺[1],
                  :k₋₁=> P.k⁻[1],
                  :n  => P.n
                 )
    expr = dict_to_parameters(params)
    eval(expr)

    r₁ = 1
    t⁺ = k₁  * x
    t⁻ = k₋₁ * (1 - x)

    A =  r₁        .* (t⁻ - t⁺)
    B = (r₁ * r₁') .* (t⁻ + t⁺) / n

    if symbolic
        return LangevinType(A, B)
    else
        AA = build_function(A, x; expression = Val{false})
        BB = build_function(B, x; expression = Val{false})

        return LangevinType(AA, BB)
    end
end
