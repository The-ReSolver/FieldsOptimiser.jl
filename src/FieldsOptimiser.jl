module FieldsOptimiser

using NSOperators
using Projector

# TODO: use "optimtest" to make sure gradient is compatible

export optimize

# TODO: add projection step to this function
function myfinalize!(x, f, g, numiter)
    # perform projection
    

    return x, f, g
end

function init_fg(cache::Cache)
    function fg(U)
        update_v!(U, cache)
        update_p!(cache)
        localresidual!(U, cache)
        update_r!(cache)
        return ℜ(cache), dℜ!(cache) 
    end
end

# TODO: how is optimize called properly?
function optimize(U₀::V, mean::NTuple{4, Vector{T}}; algorithm=ConjugateGradient()) where {T, S<:AbstractArray{Complex{T}, 3}, V<:AbstractVector{S}}
    # initialise cache
    # TODO: initialise cache properly
    _cache = 1.0

    # define optimiser function
    function ℜdℜ(U; cache=_cache)
        r = localresidual!(U, cache)
        return ℜ(r), dℜ!(r, cache)
    end

    # run optimisation
    # FIXME: method call error?
    optimize(ℜdℜ, U₀, algorithm; finalize! = myfinalize!)
end

end
