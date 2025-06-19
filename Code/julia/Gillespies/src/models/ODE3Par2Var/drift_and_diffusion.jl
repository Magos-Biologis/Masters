"""
The Drift vector and Diffusion matrix for the most basic ODE
"""
function ODE3Par2VarAB(P :: DifferentialStruct; symbolic = true)

    c = Symbolics.variables(:c, 1:2)

    vars = c
    params = Dict(
                  :k₁  => P.k[1],
                  :k₂  => P.k[2],
                  :n₁  => P.n[1],
                  :n₂  => P.n[2],
                  :w₁  => P.w[1],
                  :w₂  => P.w[2],
                 )

    @parameterification params


    r₁ = ReactionStruct( t⁺ = w₁ * c[1], t⁻ = w₂ * c[2], r = [-1; 1] )
    r₂ = ReactionStruct( t⁺ = k₁ * n₁ * c[1], t⁻ = k₁ * c[1]^2, r = [ 1; 0] )
    r₃ = ReactionStruct( t⁺ = k₂ * n₂ * c[2], t⁻ = k₂ * c[2]^2, r = [ 0; 1] )


    R = [r₁ r₂ r₃]
    A = A_i(R)
    B = Bij(R)

    # A =  r₁        .* (t⁻ - t⁺)
    # B = (r₁ * r₁') .* (t⁻ + t⁺) / n

    if symbolic
        return LangevinType(A, B), vars, params
    else
        AA, AA! = build_function(substitute(A, params), c; expression = Val{false})
        BB, BB! = build_function(substitute(B, params), c; expression = Val{false})

        return LangevinType(AA, BB)
    end
end
