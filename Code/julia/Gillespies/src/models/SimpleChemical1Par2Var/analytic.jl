function SimpleChemical1Par1VarAs(x :: V, y :: V, p :: NovelStruct, l :: MatrixLangevin) where V <: AbstractVector{<: AbstractFloat}
    u = BigFloat.(x)
    v = BigFloat.(y')



    if p.k⁺[1] == p.k⁻[1]
        C = @.  (1 - x)x
    else
        C = let a = p.k⁺[1], b = p.k⁻[1]
            @.  x - 2a / (a-b)^2 * ( (a-b)x - log(1 + (a-b)x )b )
        end
    end

    exped = @. 2 * p.n * C
    out = @. l.B(x)^-1 * exp(exped)
    return out

end
