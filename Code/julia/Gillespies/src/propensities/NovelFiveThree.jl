"""
For a simple two part chemical species A₁ → A₂, A₁ ← A₂
"""
module NovelFiveThree

export propensity,
       transitions


transitions::Vector{Vector{Int}} = [
                            [-1;1],
                            [1;-1],
                           ]


function propensity(xs::Vector, p)
    x, y = xs
    a₁ = p.k⁺[1] * x
    a₂ = p.k⁻[1] * y

    return [a₁; a₂]
end


end


