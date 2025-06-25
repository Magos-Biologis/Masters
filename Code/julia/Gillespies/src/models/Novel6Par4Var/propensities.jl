let transitions = Vector()
    local function propensity(p)::Function
        return xs -> begin
            x, y, n, b = xs
            a₁ = p.k⁺[1] * n * x
            a₋₁= p.k⁻[1] * x * x

            a₂ = p.k⁺[2] * n * x
            a₃ = p.k⁺[3] * b * x
            a₄ = p.k⁺[4] * b * y
            a₅ = p.k⁺[5] * y
            a₆ = p.k⁺[6]


            return [a₁ ;
                    a₋₁;
                    a₂ ;
                    a₃ ;
                    a₄ ;
                    a₅ ;
                    a₆ ]
        end
    end


    push!(transitions, [ 1; 0;-1; 0]) # a₁
    push!(transitions, [-1; 0; 1; 0]) # a₋₁

    push!(transitions, [ 0; 1;-1; 0]) # a₂
    push!(transitions, [-1; 0; 1; 0]) # a₃
    push!(transitions, [ 0;-1; 1;-1]) # a₄
    push!(transitions, [ 0;-1; 1; 0]) # a₅
    push!(transitions, [ 0; 0; 0; 1]) # a₆



    """
    For the novel model with 6 parameters and 4 variables
    """
    global Novel6Par4Var = PropsAndTrans(propensity,
                                         copy(transitions))
end
