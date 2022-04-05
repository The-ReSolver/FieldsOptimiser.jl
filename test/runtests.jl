using FieldsOptimiser
using Test

using ChebUtils
using Fields
using NSOperators

@testset "ℜdℜ function      " begin
    # construct velocity field (incompressible and no-slip)
    Ny = 64; Nz = 64; Nt = 64
    y = chebpts(Ny)
    Dy = chebdiff(Ny)
    Dy2 = chebddiff(Ny)
    ws = chebws(Dy)
    ω = 1.0
    β = 1.0
    grid = Grid(y, Nz, Nt, Dy, Dy2, ws, ω, β)
    u_fun(y, z, t) = sin(π*y)*exp(cos(z))*sin(t)
    v_fun(y, z, t) = (cos(π*y) + 1)*cos(z)*sin(t)
    w_fun(y, z, t) = π*sin(π*y)*sin(z)*sin(t)
    u = VectorField(PhysicalField(grid, u_fun),
                    PhysicalField(grid, v_fun),
                    PhysicalField(grid, w_fun))
    U = VectorField(grid)
    # FIXME: why does FFTW need to be imported for this line?
    # FFT! = FFTPlan!(grid; flags=FFTW.ESTIMATE)
    FFT! = FFTPlan!(grid)
    FFT!(U, u)
    Re = abs(rand())
    Ro = abs(rand())

    # construct laminar mean field
    ū_fun(y) = y
    dūdy_fun(y) = 1.0
    d2ūdy2_fun(y) = 0.0
    ū = [ū_fun(y[i]) for i in 1:Ny]
    dūdy = [dūdy_fun(y[i]) for i in 1:Ny]
    d2ūdy2 = [d2ūdy2_fun(y[i]) for i in 1:Ny]

    # initialise caches and optimisation function
    cache_opt = Cache(U[1], u[1], ū, dūdy, d2ūdy2, Re, Ro)
    cache_man = Cache(U[1], u[1], ū, dūdy, d2ūdy2, Re, Ro)
    ℜdℜ = FieldsOptimiser.init_ℜdℜ(cache_opt)

    # compute residual and gradient using optimisation function
    ℜ_opt, dℜ_opt = ℜdℜ(U)

    # compute residual and gradient using functions
    update_v!(U, cache_man)
    update_p!(cache_man)
    localresidual!(U, cache_man)
    update_r!(cache_man)
    ℜ_man, dℜ_man = ℜ(cache_man), dℜ!(cache_man)

    # compare
    @test ℜ_opt == ℜ_man
    @test dℜ_opt == dℜ_man
end

@testset "Optimize decrease " begin
    # construct velocity field (incompressible and no-slip)

    # perform single iteration of optimisation

    # check function value has decreased
end
