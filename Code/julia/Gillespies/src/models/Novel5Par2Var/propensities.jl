let transitions = Vector()
    local function propensity(p)::Function
        return xs -> begin
            x, y = xs
            a₁ = p.k⁺[1] * x
            a₂ = p.k⁻[1] * y

            return [a₁; a₂]
        end
    end

    push!(transitions, [-1;1])
    push!(transitions, [1;-1])

    """
    For the novel model with 5 parameters and 2 variables
    """
    global Novel5Par2Var = PropsAndTrans(propensity,
                                         copy(transitions))
end
