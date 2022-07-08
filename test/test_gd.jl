@testset "Gradient descent decreases objective  " begin
    # # construct velocity field (incompressible and no-slip)
    # Ny = 64; Nz = 64; Nt = 64
    # y = chebpts(Ny)
    # # Dy = chebdiff(Ny); Dy2 = chebddiff(Ny)
    # Dy = DiffMatrix(y, 5, 1); Dy2 = DiffMatrix(y, 5, 2)
    # ws = chebws(Ny)
    # ω = 1.0
    # β = 1.0
    # grid = Grid(y, Nz, Nt, Dy, Dy2, ws, ω, β)
    # u_fun(y, z, t) = sin(π*y)*exp(cos(z))*sin(t)
    # v_fun(y, z, t) = (cos(π*y) + 1)*cos(z)*sin(t)
    # w_fun(y, z, t) = π*sin(π*y)*sin(z)*sin(t)
    # u₀ = VectorField(PhysicalField(grid, u_fun),
    #                 PhysicalField(grid, v_fun),
    #                 PhysicalField(grid, w_fun))
    # U₀ = VectorField(grid)
    # FFT! = FFTPlan!(grid; flags=ESTIMATE)
    # FFT!(U₀, u₀)
    # Re = 1000.0
    # Ro = 0.5

    # # construct laminar mean field
    # ū_fun(y) = y
    # dūdy_fun(y) = 1.0
    # d2ūdy2_fun(y) = 0.0
    # ū = [ū_fun(y[i]) for i in 1:Ny]
    # dūdy = [dūdy_fun(y[i]) for i in 1:Ny]
    # d2ūdy2 = [d2ūdy2_fun(y[i]) for i in 1:Ny]

    # # initialise objective function
    # ℜdℜ = ℜdℜClosure(grid, ū, dūdy, d2ūdy2, Re, Ro)

    # # perform optimisation for limited iterations
    # _, (_, _), trace = gd!(U₀, ℜdℜ, 1e-9, 1)

    # @test trace[1][2] > trace[2][2]
end

@testset "Gradient descent projection           " begin
    # construct velocity field (incompressible and no-slip)
    Ny = 64; Nz = 64; Nt = 64
    y = chebpts(Ny)
    # Dy = chebdiff(Ny); Dy2 = chebddiff(Ny)
    Dy = DiffMatrix(y, 5, 1); Dy2 = DiffMatrix(y, 5, 2)
    ws = chebws(Ny)
    ω = 1.0
    β = 1.0
    grid = Grid(y, Nz, Nt, Dy, Dy2, ws, ω, β)
    u_fun(y, z, t) = sin(π*y)*exp(cos(z))*sin(t)
    v_fun(y, z, t) = (cos(π*y) + 1)*cos(z)*sin(t)
    w_fun(y, z, t) = π*sin(π*y)*sin(z)*sin(t)
    u₀ = VectorField(PhysicalField(grid, u_fun),
                    PhysicalField(grid, v_fun),
                    PhysicalField(grid, w_fun))
    U₀ = VectorField(grid)
    FFT! = FFTPlan!(grid; flags=ESTIMATE)
    FFT!(U₀, u₀)
    Re = 1000.0
    Ro = 0.5

    # construct laminar mean field
    ū_fun(y) = y
    dūdy_fun(y) = 1.0
    d2ūdy2_fun(y) = 0.0
    ū = [ū_fun(y[i]) for i in 1:Ny]
    dūdy = [dūdy_fun(y[i]) for i in 1:Ny]
    d2ūdy2 = [d2ūdy2_fun(y[i]) for i in 1:Ny]

    # initialise objective function
    ℜdℜ = ℜdℜClosure(grid, ū, dūdy, d2ūdy2, Re, Ro)

    # perform optimisation for limited iterations
    Uf, (_, _), trace = gd!(copy(U₀), ℜdℜ, 1e-9, 1)

    # compute divergence of iterated velocity field
    dVdy = ddy!(Uf[2], SpectralField(grid))
    dWdz = ddz!(Uf[3], SpectralField(grid))
    div = dVdy .+ dWdz

    @test norm(div) < 5e-4
    @test norm(Uf[1][1, :, :], Inf) < 1e-6
    @test norm(Uf[1][end, :, :], Inf) < 1e-6
    @test norm(Uf[2][1, :, :], Inf) < 1e-6
    @test norm(Uf[2][end, :, :], Inf) < 1e-6
    @test norm(Uf[3][1, :, :], Inf) < 1e-6
    @test norm(Uf[3][end, :, :], Inf) < 1e-6
end
