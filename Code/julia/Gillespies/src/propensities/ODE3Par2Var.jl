"""
For a simple compartment system
"""
module ODE3Par2Var

export propensity,
       transitions


transitions::Vector{Vector{Int}} = [
                            [ 1;-1], # a1
                            [-1; 1], # a-1
                            [ 1; 0], # a2
                            [-1; 0], # a-2
                            [ 0; 1], # a3
                            [ 0;-1], # a-3
                           ]


function propensity(p)::Function
    return xs -> begin
        x, y = xs

        a₁  = p.w[1] * x
        a₋₁ = p.w[2] * y

        a₂  = p.k[1] * p.n[1] * x
        a₋₂ = p.k[1] * x^2

        a₃  = p.k[2] * p.n[2] * y
        a₋₃ = p.k[2] * y^2

        return [a₁; a₋₁; a₂; a₋₂; a₃; a₋₃]
    end
end


end


