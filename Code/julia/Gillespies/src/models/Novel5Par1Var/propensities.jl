let transitions = Vector()
    function propensity(p) :: Function
        return x -> begin
            a₁ = p.k⁺[1] * p.n * x
            a₂ = p.k⁻[1] * x^2

            a₃ = p.k⁺[3] * p.b * x

            return [a₁; a₂; a₃]
        end
    end

    push!(transitions,  1)
    push!(transitions, -1)
    push!(transitions, -1)

    """
    """
    global Novel5Par1Var = PropsAndTrans(propensity,
                                         copy(transitions))
end
