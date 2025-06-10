
"""
This function modifies a given vector in place using the Dict.
By taking the value of the 'subscript' of a parameter, it places that one
in the (hopefully) correct place the vector.
Ignoring any variables that can not be assigned.
"""
function parameter_assignment!(par_vec::Vector, parameters::Dict)#, is_ode::Bool)
    vector_keys = keys(par_vec)
    for (k, v) ∈  parameters
        entry = parse(Int, match(r"(\d+)", String(k)).captures[1])
        entry ∈ vector_keys && eval(:($par_vec[$entry] = $v))
    end
end


"""
This function filters out the relevant `Dict` entries
"""
function condic(dict::Dict, param::AbstractString, len::Int)
    filtered = Dict()
    for i ∈ 1:len
        pos_key = Symbol(param * "$i")
        neg_key = Symbol(param * "_$i")
        filtered[pos_key] = get(dict, pos_key, 1.)
        filtered[neg_key] = get(dict, neg_key, 1.)
    end
    return filtered
end


"""
This function constructs the entries in the dict
"""
function filtdic(dict::Dict, param::AbstractString, len::Int; val::Real = 1)
    filtered = Dict()
    for i ∈ 1:len
        pos_key = param * "$i"
        neg_key = param * "_$i"
        filtered[pos_key] = get(dict, pos_key, val)
        filtered[neg_key] = get(dict, neg_key, val)
    end
    return filtered
end


"""
This function makes sure that the parameters fed into the simulators are
consistently sized and sorted.
It could probably be made more elegant with some for-loops
"""
function RateParameterPrimer(
        P::ParameterDict,
        ode::Bool,
        counts::ParVarStruct,
    )

    if ode
        par = P.all

        n_vec = Vector{Float64}(undef, 2); n_vec .= 1
        k_vec = Vector{Float64}(undef, 2); k_vec .= 1
        w_vec = Vector{Float64}(undef, 2); w_vec .= 1
        q_vec = Vector{Float64}(undef, 2); q_vec .= 1

        parameter_assignment!( n_vec, filtdic(par, "n", 2; val = 100) )
        parameter_assignment!( k_vec, filtdic(par, "k", 2) )
        parameter_assignment!( w_vec, filtdic(par, "w", 2) )
        parameter_assignment!( q_vec, filtdic(par, "q", 2) )

        output = n_vec, k_vec, w_vec, q_vec
    else
        pos = P.positive
        neg = P.negative

        pos_vec = Vector{Float64}(undef, counts.par); pos_vec .= 1
        neg_vec = Vector{Float64}(undef, counts.par); neg_vec .= 1

        parameter_assignment!(pos_vec, pos)
        parameter_assignment!(neg_vec, neg)

        output =  pos_vec, neg_vec
    end

    return output
end





