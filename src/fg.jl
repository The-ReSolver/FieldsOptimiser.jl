# This file contains the definition for the function that will return the
# residual and its gradient, based on the operators defined in NSOperators.jl.

function init_fg(cache::Cache)
    function fg(U)
        update_v!(U, cache)
        update_p!(cache)
        localresidual!(U, cache)
        update_r!(cache)
        return ℜ(cache), dℜ!(cache) 
    end
end
