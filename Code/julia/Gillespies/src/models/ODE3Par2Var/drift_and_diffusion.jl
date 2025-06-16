"""
The Drift vector and Diffusion matrix for the most basic ODE
"""
function ODE3Par2VarAB(P :: DifferentialStruct; symbolic = true)

    c = Symbolics.variables(:c, 1:2)

    params :: Dict{Union{Symbol, Num}, Number} = Dict(
                  :k₁  => P.k[1],
                  :k₂  => P.k[2],
                  :n₁  => P.n[1],
                  :n₂  => P.n[2],
                  :w₁  => P.w[1],
                  :w₂  => P.w[2],
                 )
    eval(@parameterification params)


    r₁ = ReactionStruct(t⁺ = k₁ * x, t⁻ = k₋₁ * (1 - x), r = 1)
    r₂ = ReactionStruct(t⁺ = k₁ * x, t⁻ = k₋₁ * (1 - x), r = 1)
    r₃ = ReactionStruct(t⁺ = k₁ * x, t⁻ = k₋₁ * (1 - x), r = 1)


    R = [r₁ r₂ r₃]
    A = A_i(R)
    B = Bij(R)

    # A =  r₁        .* (t⁻ - t⁺)
    # B = (r₁ * r₁') .* (t⁻ + t⁺) / n

    if symbolic
        return LangevinType(A, B)
    else
        AA = build_function(substitute(A, params), c; expression = Val{false})
        BB = build_function(substitute(B, params), c; expression = Val{false})

        return LangevinType(AA, BB)
    end
end
