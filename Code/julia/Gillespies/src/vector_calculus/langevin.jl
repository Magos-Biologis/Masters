
"""
A function to appropriately assign the drift vector and diffusion matrix to
a langevin equation.
"""
LangevinEquation(:: LangevinType) = error("Undefined for this type")

function LangevinEquation(Lₑ:: VectorLangevin)
    A, B = Lₑ.A, Lₑ.B
    dW = @brownian W₁ W₂ W₃ W₄ W₅
    return A + B * dW[1:length(A)]
end


function LangevinEquation(Lₑ:: ScalarLangevin)
    @brownian dW
    return Lₑ.A + Lₑ.B * dW
end







