"""
The Drift vector and Diffusion matrix for the simple chemical system
"""
function SimpleChemical1Par2VarAB(P :: NovelStruct; symbolic = true)

    @variables x::Real y::Real #-> So as to only make the variables at compile time

    vars = [x; y]
    params = Dict{Union{Symbol, Num}, Number}(
                  :k₁  => P.k⁺[1],
                  :k₋₁ => P.k⁻[1],
                  :n   => P.n
                 )
    @parameterification params


    r₁ = ReactionStruct( t⁺ = k₁ * x, t⁻ = k₋₁ * y, r = [ 1;-1])

    A = A_i(r₁)
    B = Bij(r₁)

    if symbolic
        return LangevinType(A, B, vars, params)
    else
        AA, AA! = build_function(substitute(A, params), [x, y]; expression=Val{false})
        BB, BB! = build_function(substitute(B, params), [x, y]; expression=Val{false})

        return LangevinType(AA, BB)
    end
end
