

"""
Takes in the drift and diffusion elements, alongside with the parameters and
variables returning ∇log(pₛ(x))
"""
function potentiating(L :: T, Par :: N ) where {T <: LangevinType, N <: LangevinParams}
    A, B = L
    vars = Par.variables

    eval(quote @∇ $(vars...) end)
    A_∇B = 2A - expand_derivatives.( ∇·B )

    B⁺ = (vars...) -> langevin_build( B, Par )(vars...) |> pinv
    φ  = langevin_build( A_∇B, Par )

    return (vars...) -> B⁺(vars...) * φ(vars...)
end



function Symbolics.substitute(L :: LangevinType, s :: AbstractDict)
    AA = Num.(substitute(L.A, s))
    BB = Num.(substitute(L.B, s))
    return LangevinType( AA, BB )
end
function Symbolics.substitute(L :: LangevinType, P :: LangevinParams)
    return substitute(L, P.parameters)
end



function langevin_build(L :: T, Par :: N ) where {T <: LangevinType, N <: LangevinParams}
    A, B = L
    AA = langevin_build(A, Par)
    BB = langevin_build(B, Par)
    return LangevinType(AA, BB)
end
function langevin_build( m :: AbstractVecOrMat,  Par :: N ) where N <: LangevinParams
    vars, pars = Par; M = substitute(m, pars)
    MM, M! = build_function(M, vars...; expression = false)
    return MM
end
function langevin_build( m :: Num,  Par :: N ) where N <: LangevinParams
    vars, pars = Par; M = substitute(m, pars)
    MM = build_function(M, vars...; expression = false)
    return MM
end


const BoundOrBounds{T} = Union{Vector{T}, T} where T <: Number

function GradientSol(L :: T, Par :: N;
        west_bc :: BoundOrBounds{B₁} = eps(Float64),
        east_bc :: BoundOrBounds{B₂} = 1-eps(Float64),
        Δx      :: BoundOrBounds{B₃} = 0.005,
    ) where { T <: LangevinType, N <: LangevinParams, B₁, B₂, B₃ }


    local vars = Par.variables
    local nᵥ   = length(vars)

    length(west_bc) < nᵥ ? west_bc = fill(west_bc, nᵥ) : nothing
    length(east_bc) < nᵥ ? east_bc = fill(east_bc, nᵥ) : nothing
    length(Δx)      < nᵥ ? Δx      = fill(Δx, nᵥ)      : nothing

    ranges = [ west_bc[i]:Δx[i]:east_bc[i] for i ∈ 1:nᵥ ]

    nᵥ≡ 1 && return gradient_int_one_dim(L, Par, ranges)
    nᵥ≡ 2 && return gradient_int_two_dim(L, Par, ranges)
    # nᵥ≡ 3 && return gradient_int_thr_dim(L, Par, ranges)

    return error("Invalid Model Dimensions")
end

const Steppers = Union{StepRange, StepRangeLen}

function gradient_int_one_dim(L :: T, Par :: N, range :: Vector{S};
    ) where {T <: LangevinType, N <: LangevinParams, S <: Steppers}

    x = collect( range... )
    f = potentiating(L, Par)

    ∇F = f.(x)

    Φ  = similar(∇F)
    Φ .= cumsum( ∇F .* diff( [x[1]; x] )   )


    return (x, Φ)
end

function gradient_int_two_dim(L :: T, Par :: N, ranges :: Vector{S};
    ) where {T <: LangevinType, N <: LangevinParams, S <: Steppers}

    x = ranges[1]
    y = ranges[2]

    f  = potentiating(L, Par)
    ∇F = f.(x, y')
    Fx = ∇F .|> first
    Fy = ∇F .|> last

    Δx = sum(diff(ranges[1])) / length(ranges[1])
    Δy = sum(diff(ranges[2])) / length(ranges[2])

    Φ  = similar( Fx )
    Φ .= cumsum( Fx .* Δx; dims = 2 )
    Φ += cumsum( Fy .* Δy; dims = 1 )

    return ( x, y, Φ )
end



