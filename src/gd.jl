# This file contains the definitions required to perform a gradient descent
# optimisation on a RPCF field.

# NOTE: step size >1e-7 causes blow-up (without projections)
# NOTE: step size >1e-9 causes blow-up (with Leray and slip correction)
function gd!(x, fg, α, maxiter, tol::Float64=1e-3)
    # initialise trace vector to hold iteration information
    trace = []

    # perform first evaluation outside of loop
    fval, grad = fg(x)
    normgrad = norm(grad)
    push!(trace, (copy(x), fval, normgrad))

    # if initial guess is minima, return now
    if normgrad < tol
        return x, (fval, normgrad), trace
    end

    # begin optimisation loop
    for _ in 1:maxiter
        # update field
        x -= α*grad

        # compute objective stuff
        fval, grad = fg(x)
        normgrad = norm(grad)

        # update trace
        push!(trace, (copy(x), fval, normgrad))

        # terminate if gradient is small
        if normgrad < tol
            return x, (fval, normgrad), trace
        end
    end

    return x, (fval, normgrad), trace
end
