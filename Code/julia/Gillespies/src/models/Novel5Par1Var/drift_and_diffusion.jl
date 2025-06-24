function Novel5Par1VarAB(P :: NovelStruct; symbolic = true)
    @variables x::Real

    vars = x
    params = Dict{Union{Symbol, Num}, Number}(
                  :k₁  => P.k⁺[1],  :k₋₁ => P.k⁻[1],
                  :k₃  => P.k⁺[3],
                  :n   => P.n,      :b   => P.b
                 )
    @parameterification params


    r₁ = ReactionStruct( t⁺ = k₁  * n * x, t⁻ = k₋₁ * x * x, r = 1)
    r₂ = ReactionStruct( t⁺ = k₃  * b * x, r = -1)


    R = [r₁ r₂]
    A = A_i(R)
    B = Bij(R)

    if symbolic
        return LangevinType(A, B, vars, params)
    else
        AA, AA! = build_function(substitute(A, params), vars; expression=Val{false})
        BB, BB! = build_function(substitute(B, params), vars; expression=Val{false})

        return LangevinType(AA, BB)
    end
end

