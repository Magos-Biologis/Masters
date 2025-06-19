function SimpleChemical1Par1VarAs(u :: AbstractVector, p :: NovelStruct)
    x = BigFloat.(u)
    L = SimpleChemical1Par1VarAB(p; symbolic = false)


    if p.k⁺[1] == p.k⁻[1]
        C = @.  (1 - x)x
    else
        C = let a = p.k⁺[1], b = p.k⁻[1]
            @.  x - 2a / (a-b)^2 * ( (a-b)x - log(1 + (a-b)x )b )
        end
    end

    exped = @. 2 * p.n * C
    out = @. L.B(x)^-1 * exp(exped)
    return AnalyticalSol(x, out)
end
