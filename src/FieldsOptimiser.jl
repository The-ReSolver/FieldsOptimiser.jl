module FieldsOptimiser

using OptimKit
using Printf
using LinearAlgebra

using NSOperators
using Projector

# TODO: use "optimtest" to make sure gradient is compatible
# TODO: still need to figure out way to structure this to avoid using my packages
# TODO: access fields stored in cache using proper interface (methods???)
# TODO: can construct fields using type parameter (need to implement appropriate constructor tho)

export optimize

function init_myfinalize!(leray!, slipcorrector!, cache, io::IO=stdout)
    # initialise temporary field
    div_U = similar(cache.spec_cache[1])

    function myfinalize!(U, ℜ, dℜ, numiter)
        # calculate divergence
        div_U = cache.spec_cache[11] .+ cache.spec_cache[6]
        @views begin
            slip_norm = sqrt(norm(U[1, :, :])^2 + norm(U[end, :, :])^2)
            res_norm = sqrt(norm(cache.spec_cache[36][1, :, :])^2 + norm(cache.spec_cache[36][end, :, :])^2 +
                            norm(cache.spec_cache[37][1, :, :])^2 + norm(cache.spec_cache[37][end, :, :])^2 +
                            norm(cache.spec_cache[38][1, :, :])^2 + norm(cache.spec_cache[38][end, :, :])^2)
        end

        # print relevant statistics
        @printf io "||∇U|| = %0.6f     slip_norm = %0.6f     res_norm = %0.6f\n" norm(div_U) slip_norm res_norm

        # perform projections
        leray!(U)
        slipcorrector!(U)

        return U, ℜ, dℜ
    end
end

function init_ℜdℜ(cache::C) where {C}
    function ℜdℜ(U)
        update_v!(U, cache)
        update_p!(cache)
        localresidual!(U, cache)
        update_r!(cache)
        return ℜ(cache), dℜ!(cache) 
    end
end

# TODO: can the type parameter stuff be done with a different pattern?
function optimize(U₀::V, mean::NTuple{3, Vector{T}}, Re::T, Ro::T; algorithm=ConjugateGradient()) where {T, S<:AbstractArray{Complex{T}, 3}, V<:AbstractVector{S}}
    # initialise cache
    _cache = Cache(U₀[1].grid, mean..., Re, Ro)

    # initialise projections
    leray! = Leray!(U₀[1].grid)
    slipcorrector! = SlipCorrector!(U₀[1].grid)

    # initialise finalisation function

    # initialise objective function

    # run optimisation
    optimize(ℜdℜ, U₀, algorithm; finalize! = myfinalize!)
end

end
