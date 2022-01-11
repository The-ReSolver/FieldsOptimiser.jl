# This file constains the definitions required to interface correctly with the
# OptimKit.jl (https://github.com/Jutho/OptimKit.jl) package.

# The methods are built to work with generic sub-types of abstract arrays with
# the correct size and element type.

# TODO: do these methods need to be specialised for the type?

function finalize!(x, f, g, numiter)
    # print line of information
    println("Iter: $numiter    ℜ: $f    dℜ: $(LinearAlgebra.norm(g))")
    println()

    return x, f, g
end

retract(x, η, α) = (x + α*η, η)
inner(x, ξ1, ξ2) = LinearAlgebra.dot(ξ1, ξ2)
transport!(ξ, x, η, α, xnew) = ξ
scale!(v, α) = LinearAlgebra.rmul!(v, α)
add!(η, ξ, β) = LinearAlgebra.axpy!(β, ξ, η)
