begin
    local function propensity(p)::Function
        return cs -> begin
            c₁, c₂ = cs

            a₁  = p.w[1] * c₁
            a₋₁ = p.w[2] * c₂

            a₂  = p.k[1] * p.n[1] * c₁
            a₋₂ = p.k[1] * c₁ * c₁

            a₃  = p.k[2] * p.n[2] * c₂
            a₋₃ = p.k[2] * c₂ * c₂

            a₄ = p.k[2] * p.n[2] * c₂

            a₅ = p.k[2] * c₂ * c₂

            a₆ = p.k[2] * p.n[2] * c₂

            a₇ = p.k[2] * c₂ * c₂

            a₈ = p.k[2] * c₂ * c₂

            return [a₁; a₋₁; a₂; a₋₂; a₃; a₋₃; a₄; a₅; a₆; a₇; a₈]
        end
    end

    local transitions::Vector{Vector{Int}} = [
                                        [-1; 1], # a1
                                        [ 1;-1], # a-1
                                        [ 1; 0], # a2
                                        [-1; 0], # a-2
                                        [ 0; 1], # a3
                                        [ 0;-1], # a-3
                                        [ 0; 1], # a3
                                        [ 0;-1], # a-3
                                        [ 0; 1], # a3
                                        [ 0;-1], # a-3
                                        [ 0; 1], # a3
                                        [ 0;-1], # a-3
                                        [ 0; 1], # a3
                                        [ 0;-1], # a-3
                                       ]


    """
    For a simple compartment system
    """
    global ODE8Par3Var = PropsAndTrans(propensity,
                                       copy(transitions))
end
