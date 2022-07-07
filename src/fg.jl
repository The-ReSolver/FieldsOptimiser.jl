# This file contains the definition of the objective and objective gradient
# function.

struct ℜdℜClosure{Cache, S}
    cache::Cache
    leray!::Leray!{S}
    slipcorrector!::SlipCorrector!{S}
end

function ℜdℜClosure(grid, ū::Vector{T}, dūdy::Vector{T}, d2ūdy2::Vector{T}, Re::T, Ro::T) where {T}
    length(ū) == length(dūdy) == length(d2ūdy2) == size(grid)[1] || throw(ArgumentError("Arguments are incompatible sizes"))
    ℜdℜClosure(Cache(grid, ū, dūdy, d2ūdy2, Re, Ro), Leray!(grid), SlipCorrector!(grid))
end

function (f::ℜdℜClosure)(U::AbstractVector)
    update_v!(U, f.cache)
    update_p!(f.cache)
    localresidual!(U, f.cache)
    update_r!(f.cache)
    dℜ = dℜ!(f.cache)
    f.leray!(dℜ)
    f.slipcorrector!(dℜ)
    return ℜ(f.cache), dℜ
end
