"""
For a simple two part chemical species A₁ → A₂, A₁ ← A₂
"""
module SimpleChemical1Par2Var

export propensity,
       transitions


transitions::Vector{Vector{Int}} = [
                            [-1;1],
                            [1;-1],
                           ]


function propensity(p)::Function
    return xs -> begin
        x, y = xs
        a₁ = p.k⁺[1] * x
        a₂ = p.k⁻[1] * y

        return [a₁; a₂]
    end
end


end


