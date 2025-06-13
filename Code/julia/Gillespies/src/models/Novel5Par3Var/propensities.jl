begin
    local function propensity(p)::Function
        return xs -> begin
            x, y = xs
            a₁ = p.k⁺[1] * x
            a₂ = p.k⁻[1] * y

            return [a₁; a₂]
        end
    end

    local transitions::Vector{Vector{Int}} = [
                                [-1;1],
                                [1;-1],
                               ]

    """
    For a simple two part chemical species A₁ → A₂, A₁ ← A₂
    """
    global Novel5Par3Var = PropsAndTrans(propensity,
                                         copy(transitions))
end
