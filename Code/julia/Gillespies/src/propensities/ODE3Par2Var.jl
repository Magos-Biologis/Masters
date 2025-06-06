"""
For a simple compartment system
"""
module ODE3Par2Var

export propensity,
       transitions


transitions::Vector{Vector{Int}} = [
                            [ 1; 0],
                            [-1; 0],
                            [ 0; 1],
                            [ 0;-1],
                           ]


function propensity(p)::Function
    return xs -> begin
        x, y = xs
        a₁  = p.k[1] * p.n[1] * x
        a₋₁ = p.k[1] * x^2

        a₂  = p.k[2] * p.n[2] * y
        a₋₂ = p.k[2] * y^2

        return [a₁; a₋₁; a₂; a₋₂]
    end
end


end


