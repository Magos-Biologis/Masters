


"""
The Drift vector and Diffusion matrix for the simple chemical system
"""
function SimpleChemical1Par2VarODE(P :: NovelStruct)

    @variables x(t)::Real y(t)::Real #-> So as to only make the variables at compile time
    @parameters k₁=P.k⁺[1]  k₋₁=P.k⁻[1]

    r₁ = [-1; 1]
    t⁺ = k₁  * x
    t⁻ = k₋₁ * y

    A = r₁ .* (t⁻ - t⁺)
    B = (r₁ * r₁') .* (t⁻ - t⁺)

    return LangevinStruct(A, B)
end


# ODEProblem() = ODEProblem(f::ODEFunction,u0,tspan,p=NullParameters(),callback=CallbackSet())
