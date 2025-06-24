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

            a₄ = p.k[1] * c₁ * c₂
            a₅ = p.k[2] * c₁ * c₂

            a₆ = p.m0
            a₇ = p.q[2] * p.b * c₂
            a₈ = p.q[1] * p.b * c₁

            return [a₁; a₋₁; a₂; a₋₂; a₃; a₋₃; a₄; a₅; a₆; a₇; a₈]
        end
    end

    push!(transitions, [-1; 1; 0]) # a1
    push!(transitions, [ 1;-1; 0]) # a-1
    push!(transitions, [ 1; 0; 0]) # a2
    push!(transitions, [-1; 0; 0]) # a-2
    push!(transitions, [ 0; 1; 0]) # a3
    push!(transitions, [ 0;-1; 0]) # a-3
    push!(transitions, [-1; 0; 0]) # a4
    push!(transitions, [ 0;-1; 0]) # a5
    push!(transitions, [ 0; 0; 1]) # a6
    push!(transitions, [ 0; 0;-1]) # a7
    push!(transitions, [-1; 0; 0]) # a8

    """
    For a simple compartment system
    """
    global ODE8Par3Var = PropsAndTrans(propensity,
                                       copy(transitions))
end
