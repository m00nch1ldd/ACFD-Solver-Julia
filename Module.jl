module CFDVars

# file unit numbers
const fresidual = 1
const finput    = 2
const frgid     = 3
const fflow     = 4

# Global simulation state. All fields are concretely typed so the JIT can
# specialize array accesses on the (immutable) Float64/Int element types.
mutable struct SimState

    # Block / grid dimensions
    nblocks  :: Int
    NImax    :: Int
    NJmax    :: Int
    NKmax    :: Int
    Ptsmax   :: Int

    # Run control parameters
    restart  :: Int
    iter     :: Int
    nsteps   :: Int
    nprims   :: Int
    nconserv :: Int
    rk_steps :: Int
    grid2d   :: Int
    viscous  :: Int

    # Discretization scheme selectors
    dscheme  :: Int
    dschemek :: Int
    bscheme1 :: Int
    bscheme2 :: Int
    fscheme  :: Int

    # Discretization coefficients
    adisc   :: Float64
    bdisc   :: Float64
    alpha   :: Float64
    adisck  :: Float64
    bdisck  :: Float64
    alphak  :: Float64
    alpha_f :: Float64
    alpha1  :: Float64
    alpha2  :: Float64

    # Physical parameters
    time      :: Float64
    time_step :: Float64
    Lx        :: Float64
    Ly        :: Float64
    Lz        :: Float64
    Re        :: Float64
    Mach      :: Float64
    gamma     :: Float64
    prandtl   :: Float64
    T_ref     :: Float64

    # Test case flags and diagnostics
    testcase :: Int
    taylor   :: Int
    covo     :: Int
    animfreq :: Int
    tke      :: Float64
    enstpt   :: Float64

    # Grid coordinate arrays
    xgrid :: Array{Float64, 4}
    ygrid :: Array{Float64, 4}
    zgrid :: Array{Float64, 4}

    # Flow variable arrays
    # Qp = primitive    (rho, u, v, w, P, T)
    # Qc = conservative (rho, rho*u, rho*v, rho*w, rho*E)
    Qp    :: Array{Float64, 5}
    Qc    :: Array{Float64, 5}
    Qcnew :: Array{Float64, 5}
    Qcini :: Array{Float64, 5}

    # Flux arrays
    Fflux    :: Array{Float64, 5}
    Gflux    :: Array{Float64, 5}
    Hflux    :: Array{Float64, 5}
    net_flux :: Array{Float64, 5}
    fluxD    :: Array{Float64, 5}

    # Primitive variable derivatives (viscous only)
    Qpi :: Array{Float64, 5}
    Qpj :: Array{Float64, 5}
    Qpk :: Array{Float64, 5}

    # Grid derivatives (d{x,y,z}/d{i,j,k})
    xi :: Array{Float64, 4}
    yi :: Array{Float64, 4}
    zi :: Array{Float64, 4}
    xj :: Array{Float64, 4}
    yj :: Array{Float64, 4}
    zj :: Array{Float64, 4}
    xk :: Array{Float64, 4}
    yk :: Array{Float64, 4}
    zk :: Array{Float64, 4}

    # Metric terms (curvilinear -> Cartesian) and Jacobian
    ix :: Array{Float64, 4}
    iy :: Array{Float64, 4}
    iz :: Array{Float64, 4}
    jx :: Array{Float64, 4}
    jy :: Array{Float64, 4}
    jz :: Array{Float64, 4}
    kx :: Array{Float64, 4}
    ky :: Array{Float64, 4}
    kz :: Array{Float64, 4}
    Jac :: Array{Float64, 4}

    # TDMA coefficients (D = discretization, F = filter)
    AMD :: Vector{Float64}   # sub-diagonal
    ACD :: Vector{Float64}   # diagonal
    APD :: Vector{Float64}   # super-diagonal
    AMF :: Vector{Float64}
    ACF :: Vector{Float64}
    APF :: Vector{Float64}

    # Runge-Kutta and filter coefficients, residuals
    fac_qini :: Vector{Float64}
    fac_RK   :: Vector{Float64}
    fcoeff   :: Vector{Float64}
    res      :: Vector{Float64}

    # Per-block sizes
    NI :: Vector{Int}
    NJ :: Vector{Int}
    NK :: Vector{Int}

    # Turbulence diagnostics (per-point TKE and enstrophy)
    tked :: Array{Float64, 4}
    enst :: Array{Float64, 4}

    # Pre-allocated work buffers (sized to Ptsmax = max(NI, NJ, NK)).
    # Reused inside filter_*! and tdmap! to avoid per-iteration allocations.
    filter_rhs :: Vector{Float64}
    tdma_qq    :: Vector{Float64}
    tdma_ss    :: Vector{Float64}
    tdma_fei   :: Vector{Float64}
end

# One global SimState instance, used by every routine.
# `const` makes the binding fixed (the struct itself is mutable).
const G = SimState(
    # grid dimensions
    0, 0, 0, 0, 0,
    # run control
    0, 0, 0, 0, 0, 0, 0, 0,
    # discretization schemes
    0, 0, 0, 0, 0,
    # discretization coefficients
    0., 0., 0., 0., 0., 0., 0., 0., 0.,
    # physical parameters
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    # test case flags
    0, 0, 0, 0, 0., 0.,
    # grid arrays
    zeros(1,1,1,1), zeros(1,1,1,1), zeros(1,1,1,1),
    # flow arrays
    zeros(1,1,1,1,1), zeros(1,1,1,1,1), zeros(1,1,1,1,1), zeros(1,1,1,1,1),
    # flux arrays
    zeros(1,1,1,1,1), zeros(1,1,1,1,1), zeros(1,1,1,1,1), zeros(1,1,1,1,1), zeros(1,1,1,1,1),
    # viscous arrays
    zeros(1,1,1,1,1), zeros(1,1,1,1,1), zeros(1,1,1,1,1),
    # grid derivatives
    zeros(1,1,1,1), zeros(1,1,1,1), zeros(1,1,1,1),
    zeros(1,1,1,1), zeros(1,1,1,1), zeros(1,1,1,1),
    zeros(1,1,1,1), zeros(1,1,1,1), zeros(1,1,1,1),
    # metrics
    zeros(1,1,1,1), zeros(1,1,1,1), zeros(1,1,1,1),
    zeros(1,1,1,1), zeros(1,1,1,1), zeros(1,1,1,1),
    zeros(1,1,1,1), zeros(1,1,1,1), zeros(1,1,1,1),
    zeros(1,1,1,1),
    # TDMA coefficients
    zeros(1), ones(1), zeros(1),
    zeros(1), ones(1), zeros(1),
    # RK and filter coefficients
    zeros(1), zeros(1), zeros(6), zeros(1),
    # per-block sizes
    [0], [0], [0],
    # turbulence diagnostics
    zeros(1,1,1,1), zeros(1,1,1,1),
    # work buffers
    zeros(1), zeros(1), zeros(1), zeros(1),
)

end # module CFDVars
