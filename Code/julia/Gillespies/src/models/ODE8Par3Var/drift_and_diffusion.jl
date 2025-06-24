"""
The Drift vector and Diffusion matrix for the most basic ODE
"""
function ODE8Par3VarAB(P :: DifferentialStruct; symbolic = true)

    c = variables(:c, 1:2)
    m = variable(:m)

    vars = [c...; m]
    params = Dict{Union{Symbol, Num}, Number}(
                  :k₁  => P.k[1], :k₂  => P.k[2],
                  :n₁  => P.n[1], :n₂  => P.n[2],
                  :w₁  => P.w[1], :w₂  => P.w[2],
                  :q₁  => P.q[1], :q₂  => P.q[2],
                  :m₀  => P.m₀,
                 )

    @parameterification params


    r₁ = ReactionStruct( t⁺ = w₁ * c[1], t⁻ = w₂ * c[2], r = [-1; 1; 0] )
    r₂ = ReactionStruct( t⁺ = k₁ * n₁ * c[1], t⁻ = k₁ * c[1]^2, r = [ 1; 0; 0] )
    r₃ = ReactionStruct( t⁺ = k₂ * n₂ * c[2], t⁻ = k₂ * c[2]^2, r = [ 0; 1; 0] )
    r₄ = ReactionStruct( t⁺ = k₂ * c[1] * c[2], r = [-1; 0; 0] )
    r₅ = ReactionStruct( t⁺ = k₂ * c[1] * c[2], r = [ 0;-1; 0] )
    r₆ = ReactionStruct( t⁺ = m₀, r = [ 0; 0; 1] )
    r₇ = ReactionStruct( t⁺ = q₂ * m * c[2], r = [ 0; 0;-1] )
    r₈ = ReactionStruct( t⁺ = q₁ * m * c[1], r = [-1; 0; 0] )

    R = [r₁ r₂ r₃ r₄ r₅ r₆ r₇ r₈]
    A = A_i(R)
    B = Bij(R)


    if symbolic
        return LangevinType(A, B, vars, params)
    else
        AA, AA! = build_function(substitute(A, params), vars; expression = Val{false})
        BB, BB! = build_function(substitute(B, params), vars; expression = Val{false})

        return LangevinType(AA, BB)
    end
end
