# This file contains the definitions required to perform a gradient descent
# optimisation on a RPCF field.

function gd!(x, fg, α, maxiter, tol::Float64=1e-3)
    # initialise trace vector to hold iteration information
    trace = []

    # begin optimisation loop
    for _ in 1:maxiter
        # compute objective stuff
        fval, grad = fg(x)
        normgrad = norm(grad)

        # update trace
        push!(trace, (copy(x), fval, normgrad))

        # terminate if gradient is small
        if normgrad < tol
            return x, (fval, normgrad), trace
        end

        # update field
        x -= α*grad
    end

    return x, (fval, normgrad), trace
end
