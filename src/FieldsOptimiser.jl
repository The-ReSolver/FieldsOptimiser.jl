module FieldsOptimiser

using OptimKit
using Printf
using LinearAlgebra

using NSOperators
using Projector

include("gd.jl")

# TODO: use "optimtest" to make sure gradient is compatible
# TODO: still need to figure out way to structure this to avoid using my packages
# TODO: access fields stored in cache using proper interface (methods???)
# TODO: can construct fields using type parameter (need to implement appropriate constructor tho)

export optimize

# TODO: use functor instead of function closure to increase performance of finaliser
function init_myfinalize!(leray!, slipcorrector!, cache, io::IO=stdout)
    function myfinalize!(U, ℜ, dℜ, numiter)
        # calculate residual norm at the wall
        @views begin
            res_norm = sqrt(norm(cache.spec_cache[36][1, :, :])^2 + norm(cache.spec_cache[36][end, :, :])^2 +
                            norm(cache.spec_cache[37][1, :, :])^2 + norm(cache.spec_cache[37][end, :, :])^2 +
                            norm(cache.spec_cache[38][1, :, :])^2 + norm(cache.spec_cache[38][end, :, :])^2)
        end

        # print relevant statistics
        @printf io "res_norm = %0.6f\n" res_norm

        # perform projections
        leray!(U)
        slipcorrector!(U)

        return U, ℜ, dℜ
    end
end

struct ℜdℜClosure{C}
    cache::C
end

function (f::ℜdℜClosure)(U::AbstractVector)
    print("1: ")
    update_v!(U, f.cache)
    update_p!(f.cache)
    localresidual!(U, f.cache)
    update_r!(f.cache)
    println(U[1][2, 2, 2])
    return ℜ(f.cache), dℜ!(f.cache) 
end

# TODO: slip mean data into multiple arguments to be more consistent
function OptimKit.optimize(U₀::AbstractVector{<:AbstractArray{Complex{T}, 3}}, mean::NTuple{3, Vector{T}}, Re::T, Ro::T; alg::OptimKit.OptimizationAlgorithm=ConjugateGradient()) where {T}
    # initialise cache
    _cache = Cache(U₀[1].grid, mean..., Re, Ro)

    # initialise projections
    _leray! = Leray!(U₀)
    _slipcorrector! = SlipCorrector!(U₀)

    # initialise finalisation function
    myfinalize! = init_myfinalize!(_leray!, _slipcorrector!, _cache)

    # initialise objective function
    ℜdℜ = ℜdℜClosure(_cache)

    # run optimisation
    # OptimKit.optimize(ℜdℜ, U₀, alg; finalize! = myfinalize!)
    OptimKit.optimize(ℜdℜ, U₀, alg)
end

end
