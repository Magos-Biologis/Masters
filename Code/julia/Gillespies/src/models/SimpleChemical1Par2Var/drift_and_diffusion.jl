"""
The Drift vector and Diffusion matrix for the simple chemical system
"""
function SimpleChemical1Par2VarAB(P :: NovelStruct; symbolic = false)

    @variables x::Real y::Real #-> So as to only make the variables at compile time
    # @parameters begin
        k₁ = P.k⁺[1]
        k₋₁= P.k⁻[1]
        n  = P.n
    # end

    r₁ = [ 1;-1]
    t⁺ = k₁  * x
    t⁻ = k₋₁ * y

    A =  r₁        .* (t⁻ - t⁺)
    B = (r₁ * r₁') .* (t⁻ + t⁺)


    if symbolic
        return LangevinType(A, B)
    else
        AA, AA! = build_function(A, x, y; expression=Val{false})
        BB, BB! = build_function(B, x, y; expression=Val{false})

        return LangevinType(AA, BB)
    end
end
