function Novel5Par1VarAB(P :: NovelStruct)

    @variables x(t)::Real
    @parameters k₁=P.k⁺[1]  k₋₁=P.k⁻[1] k₃=P.k⁺[3] n=P.n b=P.b

    r₁  = 1
    t⁺₁ = k₁  * n * x
    t⁻₁ = k₋₁ * x * x

    r₂  = -1
    t₂ = k₃ * b * x

    A =  r₁        .* (t⁻₁ - t⁺₁) +  r₂        .* t₂
    B = (r₁ * r₁') .* (t⁻₁ + t⁺₁) + (r₂ * r₂') .* t₂

    return LangevinType(A, B)
end
