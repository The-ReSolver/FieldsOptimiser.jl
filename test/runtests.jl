using FieldsOptimiser

using Test
using Random
using OptimKit

using ChebUtils
using Fields
using NSOperators
using Projector
using FDGrids

include("test_gd.jl")

# FIXME: this test will now fail, and is kinda pointless
@testset "ℜdℜ function                          " begin
    # # construct velocity field (incompressible and no-slip)
    # Ny = 64; Nz = 64; Nt = 64
    # y = chebpts(Ny)
    # Dy = chebdiff(Ny)
    # Dy2 = chebddiff(Ny)
    # ws = chebws(Dy)
    # ω = 1.0
    # β = 1.0
    # grid = Grid(y, Nz, Nt, Dy, Dy2, ws, ω, β)
    # u_fun(y, z, t) = sin(π*y)*exp(cos(z))*sin(t)
    # v_fun(y, z, t) = (cos(π*y) + 1)*cos(z)*sin(t)
    # w_fun(y, z, t) = π*sin(π*y)*sin(z)*sin(t)
    # u = VectorField(PhysicalField(grid, u_fun),
    #                 PhysicalField(grid, v_fun),
    #                 PhysicalField(grid, w_fun))
    # U = VectorField(grid)
    # FFT! = FFTPlan!(grid; flags=ESTIMATE)
    # FFT!(U, u)
    # Re = abs(rand())
    # Ro = abs(rand())

    # # construct laminar mean field
    # ū_fun(y) = y
    # dūdy_fun(y) = 1.0
    # d2ūdy2_fun(y) = 0.0
    # ū = [ū_fun(y[i]) for i in 1:Ny]
    # dūdy = [dūdy_fun(y[i]) for i in 1:Ny]
    # d2ūdy2 = [d2ūdy2_fun(y[i]) for i in 1:Ny]

    # # initialise caches and optimisation function
    # cache_opt = Cache(U[1], u[1], ū, dūdy, d2ūdy2, Re, Ro)
    # cache_man = Cache(U[1], u[1], ū, dūdy, d2ūdy2, Re, Ro)
    # ℜdℜ = FieldsOptimiser.ℜdℜClosure(cache_opt)

    # # compute residual and gradient using optimisation function
    # ℜ_opt, dℜ_opt = ℜdℜ(U)

    # # compute residual and gradient using functions
    # update_v!(U, cache_man)
    # update_p!(cache_man)
    # localresidual!(U, cache_man)
    # update_r!(cache_man)
    # ℜ_man, dℜ_man = ℜ(cache_man), dℜ!(cache_man)

    # # compare
    # @test ℜ_opt == ℜ_man
    # @test dℜ_opt == dℜ_man
end

# FIXME: projection now done in objective evaluation rather than optimisation finalisation
@testset "Projection finalisation               " begin
    # # initialise velocity field (compressible and slipping)
    # Ny = 64; Nz = 64; Nt = 64
    # y = chebpts(Ny)
    # Dy = chebdiff(Ny)
    # Dy2 = chebddiff(Ny)
    # ws = chebws(Dy)
    # ω = 1.0
    # β = 1.0
    # u_fun(y, z, t) = (y^2)*exp(cos(z))*atan(sin(t))
    # v_fun(y, z, t) = sin(π*y)*exp(sin(z))*atan(cos(t))
    # w_fun(y, z, t) = cos(π*y)*cos(z)*exp(cos(t))
    # grid = Grid(y, Nz, Nt, Dy, Dy2, ws, ω, β)
    # u = VectorField(PhysicalField(grid, u_fun),
    #                 PhysicalField(grid, v_fun),
    #                 PhysicalField(grid, w_fun))
    # U = VectorField(grid)
    # FFT! = FFTPlan!(grid; flags=ESTIMATE)
    # FFT!(U, u)

    # # calculate divergence of velocity field
    # dVdy = SpectralField(grid)
    # dWdz = SpectralField(grid)
    # ddy!(U[2], dVdy)
    # ddz!(U[3], dWdz)
    # divU = dVdy + dWdz
    
    # # initialise projections and cache
    # leray! = Leray!(grid)
    # slipcorrector! = SlipCorrector!(grid)
    # cache = Cache(grid, rand(Ny), rand(Ny), rand(Ny), 1.0, 1.0)

    # # initialise finaliser function and objective function
    # tmp_io = IOBuffer()
    # myfinalize! = FieldsOptimiser.init_myfinalize!(leray!, slipcorrector!, cache, tmp_io)
    # ℜdℜ = FieldsOptimiser.ℜdℜClosure(cache)

    # # calculate residual and gradient
    # ℜ, dℜ = ℜdℜ(U)

    # # run finalizer
    # U, ℜ, dℜ = myfinalize!(U, ℜ, dℜ, 1)

    # # calculate divergence of new velocity field
    # ddy!(U[2], dVdy)
    # ddz!(U[3], dWdz)
    # divU = dVdy + dWdz
    # @test norm(divU) < 1e-9

    # # test slip
    # @test U[1][1, :, :] ≈ zeros((Nz >> 1) + 1, Nt) atol=1e-12
    # @test U[1][end, :, :] ≈ zeros((Nz >> 1) + 1, Nt) atol=1e-12
    # @test U[2][1, :, :] ≈ zeros((Nz >> 1) + 1, Nt) atol=1e-12
    # @test U[2][end, :, :] ≈ zeros((Nz >> 1) + 1, Nt) atol=1e-12
    # @test U[3][1, :, :] ≈ zeros((Nz >> 1) + 1, Nt) atol=1e-7
    # @test U[3][end, :, :] ≈ zeros((Nz >> 1) + 1, Nt) atol=1e-7
end

