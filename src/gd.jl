# This file contains the definitions required to perform a gradient descent
# optimisation on a RPCF field.

# NOTE: step size >1e-7 causes blow-up (without projections)
# NOTE: step size >1e-9 causes blow-up (with Leray and slip correction)
function gd!(x, fg, α, maxiter, tol::Float64=1e-3; verbose::Bool=true)
    # initialise trace vector to hold iteration information
    trace = []

    # perform first evaluation outside of loop
    fval, grad = fg(x)
    normgrad = norm(grad)
    push!(trace, (copy(x), fval, normgrad))

    # print titles of states
    if verbose
        @printf "Iter      f              |∇f|           fᵢ/fᵢ₋₁\n"
        @printf "-----------------------------------------------\n"
        @printf "0         %.6e   %.6e                          \n" fval normgrad
    end

    # if initial guess is minima, return now
    if normgrad < tol
        return x, (fval, normgrad), trace
    end

    # begin optimisation loop
    for i in 1:maxiter
        # update field
        x -= α*grad

        # compute objective stuff
        fval, grad = fg(x)
        normgrad = norm(grad)

        # update trace
        push!(trace, (copy(x), fval, normgrad))

        # print state
        if verbose
            @printf "%d         %.6e   %.6e   %.6e\n" i fval normgrad fval/trace[i][2]
        end

        # terminate if gradient is small
        if normgrad < tol
            return x, (fval, normgrad), trace
        end
    end

    return x, (fval, normgrad), trace
end
