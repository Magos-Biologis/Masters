


"""
The Drift vector and Diffusion matrix for the simple chemical system
"""
function SimpleChemical1Par2Var(P :: NovelStructCalc{T} ) :: LangevinStruct{T} where T <: Real

    @variables x::Real y::Real

    r₁ = [-1; 1]
    t⁺ = P.k⁺[1] * x
    t⁻ = P.k⁻[2] * y

    φ = r₁ * (t⁻ - t⁺)
    A = build_function(φ, x, y)

    return eval(A)
end


