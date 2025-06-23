function Novel5Par2VarAB(P :: NovelStruct; symbolic = true)
    @variables x::Real y::Real

    vars = [x; y]
    params = Dict{Union{Symbol, Num}, Number}(
                  :k₁  => P.k⁺[1],  :k₋₁ => P.k⁻[1],
                  :k₂  => P.k⁺[2],
                  :k₃  => P.k⁺[3],
                  :k₄  => P.k⁺[4],
                  :k₅  => P.k⁺[5],
                  :n   => P.n,      :b   => P.b
                 )
    @parameterification params


    r₁ = ReactionStruct( t⁺ = k₁ * n * x, t⁻ = k₋₁ * x * x, r = [ 1; 0])
    r₂ = ReactionStruct( t⁺ = k₂ * n * x, r = [ 0; 1])
    r₃ = ReactionStruct( t⁺ = k₃ * b * x, r = [-1; 0])
    r₄ = ReactionStruct( t⁺ = k₄ * b * x, r = [ 0;-1])
    r₅ = ReactionStruct( t⁺ = k₅ * y,     r = [ 0;-1])


    R = [r₁ r₂ r₃ r₄ r₅]
    A = A_i(R)
    B = Bij(R)

    if symbolic
        return LangevinType(A, B, vars, params)
    else
        AA, AA! = build_function(substitute(A, params), [vars]; expression=Val{false})
        BB, BB! = build_function(substitute(B, params), [vars]; expression=Val{false})

        return LangevinType(AA, BB)
    end
end
