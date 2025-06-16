"""
The Drift vector and Diffusion matrix for the simple chemical system
"""
function SimpleChemical1Par2VarAB(P :: NovelStruct; symbolic = true)

    @variables x::Real y::Real #-> So as to only make the variables at compile time

    params :: Dict{Union{Symbol, Num}, Number} = Dict(
                  :k₁  => P.k⁺[1],
                  :k₋₁ => P.k⁻[1],
                  :n   => P.n
                 )
    eval(@parameterification params)


    r₁ = [ 1;-1]
    t⁺ = k₁  * x
    t⁻ = k₋₁ * y

    A =  r₁        .* (t⁻ - t⁺)
    B = (r₁ * r₁') .* (t⁻ + t⁺)


    if symbolic
        return LangevinType(A, B)
    else
        AA, AA! = build_function(substitute(A, params), [x, y]; expression=Val{false})
        BB, BB! = build_function(substitute(B, params), [x, y]; expression=Val{false})

        return LangevinType(AA, BB)
    end
end
