let transitions = Vector()
    local function propensity(p)::Function
        return xs -> begin
            x, y = xs
            a₁ = p.k⁺[1] * x
            a₂ = p.k⁻[1] * y

            return [a₁; a₂]
        end
    end

    push!(transitions, [-1; 1])
    push!(transitions, [ 1;-1])

    """
    For a simple two part chemical species A₁ → A₂, A₁ ← A₂
    """
    global Novel5Par2Var = PropsAndTrans(propensity,
                                         copy(transitions))
end


