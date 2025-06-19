let transitions = Vector()
    local function propensity(p)::Function
        return cs -> begin
            c₁, c₂ = cs

            a₁  = p.w[1] * c₁
            a₋₁ = p.w[2] * c₂

            a₂  = p.k[1] * p.n[1] * c₁
            a₋₂ = p.k[1] * c₁ * c₁

            a₃  = p.k[2] * p.n[2] * c₂
            a₋₃ = p.k[2] * c₂ * c₂

            return [a₁; a₋₁; a₂; a₋₂; a₃; a₋₃]
        end
    end

    push!(transitions, [-1; 1]) # a1
    push!(transitions, [ 1;-1]) # a-1
    push!(transitions, [ 1; 0]) # a2
    push!(transitions, [-1; 0]) # a-2
    push!(transitions, [ 0; 1]) # a3
    push!(transitions, [ 0;-1]) # a-3

    """
    For a simple compartment system
    """
    global ODE3Par2Var = PropsAndTrans(propensity,
                                       copy(transitions))

end
