@testset "Gradient descent          " begin
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
    # NOTE: step size >1e-7 causes blow-up (without projections)
    # NOTE: step size >1e-9 causes blow-up (with Leray and slip correction)
    Uf, (ℜf, dℜf), trace = gd!(copy(U₀), ℜdℜ, 1e-9, 20)

    # for i in 1:(length(trace) - 1)
    #     println(trace[i][2], " ", trace[i][3], " ", trace[i + 1][2]/trace[i][2])
    # end
end
