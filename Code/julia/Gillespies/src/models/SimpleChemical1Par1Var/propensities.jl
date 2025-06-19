let transitions = Vector()
    local function propensity(p)::Function
        return x -> begin
            a₁ = p.k⁺[1] * x
            a₂ = p.k⁻[1] * (p.n - x)

            return [a₁; a₂]
        end
    end

    push!(transitions, [-1])
    push!(transitions, [ 1])

    """
    For a simple two part chemical species A₁ → A₂, A₁ ← A₂,
    such that the second species is substituted with the first
    """
    global SimpleChemical1Par1Var = PropsAndTrans(propensity,
                                                  copy(transitions))
end

