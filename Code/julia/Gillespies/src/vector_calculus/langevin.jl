
function LangevinEquation(Lₑ:: LangevinStruct)

    A, B = Lₑ.A, Lₑ.B
    dW = @brownian W₁ W₂ W₃ W₄ W₅

    return A + B * dW[1:length(A)]
end
