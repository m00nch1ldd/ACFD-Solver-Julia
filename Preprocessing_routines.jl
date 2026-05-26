# Preprocessing: read input, allocate arrays, build grid, set initial state,
# compute discretization/filter/RK coefficients, and compute metric terms.

include("Module.jl")
using .CFDVars
import .CFDVars: G

include("cyclictdma.jl")
include("Solver_routines.jl")


# ---------------------------------------------------------------------------
# READ INPUT
# ---------------------------------------------------------------------------
function read_input!()

    open("input.dat", "r") do io

        # line 1: restart, nblocks
        line = split(readline(io))
        G.restart = parse(Int, line[1])
        G.nblocks = parse(Int, line[2])

        # line 2: NImax, NJmax, NKmax
        line = split(readline(io))
        G.NImax = parse(Int, line[1])
        G.NJmax = parse(Int, line[2])
        G.NKmax = parse(Int, line[3])

        # line 3: testcase, viscous
        line = split(readline(io))
        G.testcase = parse(Int, line[1])
        G.viscous  = parse(Int, line[2])

        # line 4: Re, Mach, gamma, prandtl, T_ref
        line = split(readline(io))
        G.Re      = parse(Float64, line[1])
        G.Mach    = parse(Float64, line[2])
        G.gamma   = parse(Float64, line[3])
        G.prandtl = parse(Float64, line[4])
        G.T_ref   = parse(Float64, line[5])

        # line 5: nprims, nconserv
        line = split(readline(io))
        G.nprims   = parse(Int, line[1])
        G.nconserv = parse(Int, line[2])

        # line 6: dscheme, fscheme, alpha_f
        line = split(readline(io))
        G.dscheme = parse(Int,     line[1])
        G.fscheme = parse(Int,     line[2])
        G.alpha_f = parse(Float64, line[3])

        # line 7: rk_steps, nsteps, time_step, animfreq
        line = split(readline(io))
        G.rk_steps  = parse(Int,     line[1])
        G.nsteps    = parse(Int,     line[2])
        G.time_step = parse(Float64, line[3])
        G.animfreq  = parse(Int,     line[4])

        # line 8: execution mode selector
        # exec_mode = 0 -> serial
        # exec_mode = 1 -> parallel (parallel_mode selects backend)
        # parallel_mode = 1 -> async tasks/coroutines
        # parallel_mode = 2 -> multi-threading
        line = split(readline(io))
        G.exec_mode     = parse(Int, line[1])
        G.parallel_mode = parse(Int, line[2])
    end

    # derived quantities
    if G.exec_mode == 0
        G.parallel_mode = 0
    elseif G.exec_mode == 1
        if !(G.parallel_mode in (1, 2))
            error("parallel_mode must be 1 (async tasks) or 2 (multi-threading).")
        end
    else
        error("exec_mode must be 0 (serial) or 1 (parallel).")
    end
    G.grid2d = (G.NKmax == 1) ? 1 : 0
    G.Ptsmax = max(G.NImax, G.NJmax, G.NKmax)

    println("Reading input file Done")
    println("NImax = ", G.NImax)
    println("Re    = ", G.Re)
    println("Mach  = ", G.Mach)
end

read_input!()


# ---------------------------------------------------------------------------
# ALLOCATE ARRAYS
# ---------------------------------------------------------------------------
function allocate_routine!()

    ni = G.NImax
    nj = G.NJmax
    nk = G.NKmax
    nb = G.nblocks
    np = G.nprims
    nc = G.nconserv
    pm = G.Ptsmax

    # per-block sizes
    G.NI = fill(ni, nb)
    G.NJ = fill(nj, nb)
    G.NK = fill(nk, nb)

    # grid arrays
    G.xgrid = zeros(ni, nj, nk, nb)
    G.ygrid = zeros(ni, nj, nk, nb)
    G.zgrid = zeros(ni, nj, nk, nb)

    # flow arrays
    G.Qp    = zeros(ni, nj, nk, nb, np)
    G.Qc    = zeros(ni, nj, nk, nb, nc)
    G.Qcini = zeros(ni, nj, nk, nb, nc)
    G.Qcnew = zeros(ni, nj, nk, nb, nc)

    # flux arrays
    G.Fflux    = zeros(ni, nj, nk, nb, nc)
    G.Gflux    = zeros(ni, nj, nk, nb, nc)
    G.Hflux    = zeros(ni, nj, nk, nb, nc)
    G.net_flux = zeros(ni, nj, nk, nb, nc)
    G.fluxD    = zeros(ni, nj, nk, nb, nc)

    # viscous arrays
    if G.viscous == 1
        G.Qpi = zeros(ni, nj, nk, nb, np)
        G.Qpj = zeros(ni, nj, nk, nb, np)
        G.Qpk = zeros(ni, nj, nk, nb, np)
    end

    # grid derivative arrays
    G.xi = zeros(ni, nj, nk, nb)
    G.yi = zeros(ni, nj, nk, nb)
    G.zi = zeros(ni, nj, nk, nb)
    G.xj = zeros(ni, nj, nk, nb)
    G.yj = zeros(ni, nj, nk, nb)
    G.zj = zeros(ni, nj, nk, nb)
    G.xk = zeros(ni, nj, nk, nb)
    G.yk = zeros(ni, nj, nk, nb)
    G.zk = zeros(ni, nj, nk, nb)

    # metric arrays + Jacobian
    G.ix = zeros(ni, nj, nk, nb)
    G.iy = zeros(ni, nj, nk, nb)
    G.iz = zeros(ni, nj, nk, nb)
    G.jx = zeros(ni, nj, nk, nb)
    G.jy = zeros(ni, nj, nk, nb)
    G.jz = zeros(ni, nj, nk, nb)
    G.kx = zeros(ni, nj, nk, nb)
    G.ky = zeros(ni, nj, nk, nb)
    G.kz = zeros(ni, nj, nk, nb)
    G.Jac = zeros(ni, nj, nk, nb)

    # TDMA coefficients (discretization and filter)
    G.AMD = zeros(pm)
    G.ACD = ones(pm)
    G.APD = zeros(pm)
    G.AMF = zeros(pm)
    G.ACF = ones(pm)
    G.APF = zeros(pm)

    # RK and filter coefficients, residuals
    G.fac_qini = zeros(G.rk_steps)
    G.fac_RK   = zeros(G.rk_steps)
    G.fcoeff   = zeros(6)
    G.res      = zeros(nc)

    # turbulence diagnostics
    G.tked = zeros(ni, nj, nk, nb)
    G.enst = zeros(ni, nj, nk, nb)

    # work buffers (reused inside filter_*! and tdmap!)
    G.filter_rhs = zeros(pm)
    G.tdma_qq    = zeros(pm)
    G.tdma_ss    = zeros(pm)
    G.tdma_fei   = zeros(pm)

    println("Allocating arrays Done")
    println("xgrid size = ", size(G.xgrid))
    println("Qp size    = ", size(G.Qp))
end

allocate_routine!()


# ---------------------------------------------------------------------------
# GENERATE GRID
# ---------------------------------------------------------------------------
function generate_grid!()

    # testcase = 1 : Taylor-Green vortex (3D, periodic, uniform)
    if G.testcase == 1

        G.Lx = 2.0 * 22.0 / 7.0
        G.Ly = 2.0 * 22.0 / 7.0
        G.Lz = 2.0 * 22.0 / 7.0

        @inbounds for nbl = 1:G.nblocks
            for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]
                G.xgrid[i,j,k,nbl] = G.Lx*(i-1.0)/(G.NI[nbl]-1.0)
                G.ygrid[i,j,k,nbl] = G.Ly*(j-1.0)/(G.NJ[nbl]-1.0)
                G.zgrid[i,j,k,nbl] = G.Lz*(k-1.0)/(G.NK[nbl]-1.0)
            end
        end

    # testcase = 2 : CoVo uniform Cartesian, centered at origin (2D)
    elseif G.testcase == 2

        G.Lx = 16.0
        G.Ly = 16.0
        G.Lz = 0.0

        @inbounds for nbl = 1:G.nblocks
            for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]
                G.xgrid[i,j,k,nbl] = G.Lx*(i-1.0)/(G.NI[nbl]-1.0) - G.Lx/2.0
                G.ygrid[i,j,k,nbl] = G.Ly*(j-1.0)/(G.NJ[nbl]-1.0) - G.Ly/2.0
                G.zgrid[i,j,k,nbl] = 0.0
            end
        end

    # testcase = 3 : CoVo randomly perturbed grid inside [-5,5]^2 (2D).
    # Julia rand() differs from Fortran rand(); grid is not bit-identical.
    elseif G.testcase == 3

        G.Lx = 12.0
        G.Ly = 12.0
        G.Lz = 0.0

        @inbounds for nbl = 1:G.nblocks
            for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]
                xg = G.Lx*(i-1.0)/(G.NI[nbl]-1.0) - G.Lx/2.0
                yg = G.Ly*(j-1.0)/(G.NJ[nbl]-1.0) - G.Ly/2.0
                if -5.0 <= xg <= 5.0 && -5.0 <= yg <= 5.0
                    xg += 0.2 * (G.Lx/(G.NI[nbl]-1.0)) * (rand() - 0.5)
                    yg += 0.2 * (G.Ly/(G.NJ[nbl]-1.0)) * (rand() - 0.5)
                end
                G.xgrid[i,j,k,nbl] = xg
                G.ygrid[i,j,k,nbl] = yg
                G.zgrid[i,j,k,nbl] = 0.0
            end
        end

    # testcase = 4 : CoVo sinusoidally distorted grid (2D)
    elseif G.testcase == 4

        G.Lx = 12.0
        G.Ly = 12.0
        G.Lz = 0.0
        sixπ = 6.0 * π

        @inbounds for nbl = 1:G.nblocks
            for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]
                dxi = G.Lx/(G.NI[nbl]-1.0)
                dyj = G.Ly/(G.NJ[nbl]-1.0)
                G.xgrid[i,j,k,nbl] = G.Lx*(i-1.0)/(G.NI[nbl]-1.0) - G.Lx/2.0 +
                                     dxi * 1.0 * sin(sixπ * (j-1.0)/(G.NJ[nbl]-1.0))
                G.ygrid[i,j,k,nbl] = G.Ly*(j-1.0)/(G.NJ[nbl]-1.0) - G.Ly/2.0 +
                                     dyj * 2.0 * sin(sixπ * (i-1.0)/(G.NI[nbl]-1.0))
                G.zgrid[i,j,k,nbl] = 0.0
            end
        end
    end

    println("Grid Generation Done")
    println("xgrid min = ", minimum(G.xgrid))
    println("xgrid max = ", maximum(G.xgrid))
end

generate_grid!()


# ---------------------------------------------------------------------------
# INITIAL CONDITIONS  (TGV and CoVo)
# ---------------------------------------------------------------------------
function initialize_non_dimensionalize!()

    γ = G.gamma
    M = G.Mach

    @inbounds for nbl = 1:G.nblocks
        for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]

            xl = G.xgrid[i,j,k,nbl]
            yl = G.ygrid[i,j,k,nbl]
            zl = G.zgrid[i,j,k,nbl]

            if G.testcase == 1
                # Taylor-Green vortex
                G.Qp[i,j,k,nbl,2] =  sin(xl)*cos(yl)*cos(zl)
                G.Qp[i,j,k,nbl,3] = -cos(xl)*sin(yl)*cos(zl)
                G.Qp[i,j,k,nbl,4] =  0.0
                G.Qp[i,j,k,nbl,5] = (1.0/(γ*M^2)) +
                                    (1.0/16.0)*(cos(2.0*xl)+cos(2.0*yl)) *
                                    (cos(2.0*zl)+2.0)
                G.Qp[i,j,k,nbl,6] = 1.0
                G.Qp[i,j,k,nbl,1] = γ*M^2 * G.Qp[i,j,k,nbl,5] / G.Qp[i,j,k,nbl,6]

            elseif G.testcase == 2 || G.testcase == 3 || G.testcase == 4
                # Isentropic convected vortex (shared by cases 2, 3, 4)
                vorstr = 0.02
                rad2   = xl*xl + yl*yl
                er2h   = exp(-0.5 * rad2)
                er2    = exp(-rad2)

                G.Qp[i,j,k,nbl,2] = 1.0 - vorstr*yl*er2h
                G.Qp[i,j,k,nbl,3] =       vorstr*xl*er2h
                G.Qp[i,j,k,nbl,4] = 0.0
                G.Qp[i,j,k,nbl,1] = 1.0
                G.Qp[i,j,k,nbl,5] = (1.0/(γ*M^2)) -
                                    G.Qp[i,j,k,nbl,1] * vorstr*vorstr * er2 * 0.5
                G.Qp[i,j,k,nbl,6] = γ*M^2 * G.Qp[i,j,k,nbl,5] / G.Qp[i,j,k,nbl,1]
            end

            # Conservative variables (common to all cases)
            ρ = G.Qp[i,j,k,nbl,1]
            u = G.Qp[i,j,k,nbl,2]
            v = G.Qp[i,j,k,nbl,3]
            w = G.Qp[i,j,k,nbl,4]
            T = G.Qp[i,j,k,nbl,6]

            Etotal = T / (γ*(γ-1.0)*M^2) + 0.5*(u^2 + v^2 + w^2)

            G.Qc[i,j,k,nbl,1] = ρ
            G.Qc[i,j,k,nbl,2] = ρ*u
            G.Qc[i,j,k,nbl,3] = ρ*v
            G.Qc[i,j,k,nbl,4] = ρ*w
            G.Qc[i,j,k,nbl,5] = ρ*Etotal
        end
    end

    println("Initialization Done")
    println("u max = ", maximum(G.Qp[:,:,:,:,2]))
    println("rho   = ", G.Qp[1,1,1,1,1])
end

initialize_non_dimensionalize!()


# ---------------------------------------------------------------------------
# DISCRETIZATION, FILTER, AND RK COEFFICIENTS
# ---------------------------------------------------------------------------
function discretization_filter_rk_vals!()

    # discretization coefficients
    if G.dscheme == 1
        G.alpha = 0.0;        G.adisc = 1.0;        G.bdisc = 0.0
    elseif G.dscheme == 2
        G.alpha = 0.0;        G.adisc = 4.0/3.0;    G.bdisc = -1.0/3.0
    elseif G.dscheme == 3
        G.alpha = 1.0/4.0;    G.adisc = 3.0/2.0;    G.bdisc = 0.0
    elseif G.dscheme == 4
        G.alpha = 1.0/3.0;    G.adisc = 14.0/9.0;   G.bdisc = 1.0/9.0
    end

    fill!(G.AMD, G.alpha)
    fill!(G.ACD, 1.0)
    fill!(G.APD, G.alpha)

    # RK4 coefficients
    G.fac_qini .= [0.0,     0.5,     0.5,     1.0]
    G.fac_RK   .= [1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0]

    # filter coefficients
    af = G.alpha_f
    fs = G.fscheme

    if fs == 2
        G.fcoeff[1] = 0.5 + af
        G.fcoeff[2] = 0.5 + af
        G.fcoeff[3] = 0.0
        G.fcoeff[4] = 0.0
        G.fcoeff[5] = 0.0
        G.fcoeff[6] = 0.0

    elseif fs == 4
        G.fcoeff[1] = (5.0/8.0)  + (3.0/4.0)*af
        G.fcoeff[2] = (1.0/2.0)  + af
        G.fcoeff[3] = (-1.0/8.0) + (1.0/4.0)*af
        G.fcoeff[4] = 0.0
        G.fcoeff[5] = 0.0
        G.fcoeff[6] = 0.0

    elseif fs == 6
        G.fcoeff[1] =  (11.0/16.0) + (5.0/8.0)*af
        G.fcoeff[2] =  (15.0/32.0) + (17.0/16.0)*af
        G.fcoeff[3] = -(3.0/16.0)  + (3.0/8.0)*af
        G.fcoeff[4] =  (1.0/32.0)  - (1.0/16.0)*af
        G.fcoeff[5] =  0.0
        G.fcoeff[6] =  0.0

    elseif fs == 8
        G.fcoeff[1] =  (93.0/128.0)  + (70.0/128.0)*af
        G.fcoeff[2] =  (7.0/16.0)    + (18.0/16.0)*af
        G.fcoeff[3] = -(7.0/32.0)    + (14.0/32.0)*af
        G.fcoeff[4] =  (1.0/16.0)    - (1.0/8.0)*af
        G.fcoeff[5] = -(1.0/128.0)   + (1.0/64.0)*af
        G.fcoeff[6] =  0.0

    elseif fs == 10
        G.fcoeff[1] =  (193.0/256.0) + (126.0/256.0)*af
        G.fcoeff[2] =  (105.0/256.0) + (302.0/256.0)*af
        G.fcoeff[3] =  (15.0/64.0)   * (-1.0 + 2.0*af)
        G.fcoeff[4] =  (45.0/512.0)  * ( 1.0 - 2.0*af)
        G.fcoeff[5] =  (5.0/256.0)   * (-1.0 + 2.0*af)
        G.fcoeff[6] =  (1.0/512.0)   * ( 1.0 - 2.0*af)
    end

    fill!(G.AMF, af)
    fill!(G.ACF, 1.0)
    fill!(G.APF, af)

    println("Discretization Filter RK coefficients Done")
    println("adisc  = ", G.adisc)
    println("fac_RK = ", G.fac_RK)
end

discretization_filter_rk_vals!()


# ---------------------------------------------------------------------------
# METRICS  (grid derivatives -> Jacobian and inverse metrics)
# ---------------------------------------------------------------------------
function metrics!()

    # Step 1: compute grid derivatives (xi, yi, zi, xj, yj, zj, xk, yk, zk)
    if G.dscheme == 1 || G.dscheme == 2
        discretization_i_exp_grid!(G.xgrid, G.xi)
        discretization_i_exp_grid!(G.ygrid, G.yi)
        discretization_i_exp_grid!(G.zgrid, G.zi)

        discretization_j_exp_grid!(G.xgrid, G.xj)
        discretization_j_exp_grid!(G.ygrid, G.yj)
        discretization_j_exp_grid!(G.zgrid, G.zj)

        if G.grid2d != 1
            discretization_k_exp_grid!(G.xgrid, G.xk)
            discretization_k_exp_grid!(G.ygrid, G.yk)
            discretization_k_exp_grid!(G.zgrid, G.zk)
        else
            fill!(G.xk, 0.0)
            fill!(G.yk, 0.0)
            fill!(G.zk, 1.0)
        end

    elseif G.dscheme == 3 || G.dscheme == 4
        discretization_i_comp_grid!(G.xgrid, G.xi)
        discretization_i_comp_grid!(G.ygrid, G.yi)
        discretization_i_comp_grid!(G.zgrid, G.zi)

        discretization_j_comp_grid!(G.xgrid, G.xj)
        discretization_j_comp_grid!(G.ygrid, G.yj)
        discretization_j_comp_grid!(G.zgrid, G.zj)

        if G.grid2d != 1
            discretization_k_comp_grid!(G.xgrid, G.xk)
            discretization_k_comp_grid!(G.ygrid, G.yk)
            discretization_k_comp_grid!(G.zgrid, G.zk)
        else
            fill!(G.xk, 0.0)
            fill!(G.yk, 0.0)
            fill!(G.zk, 1.0)
        end
    end

    println("xi min = ", minimum(G.xi), "  max = ", maximum(G.xi))
    println("xj min = ", minimum(G.xj), "  max = ", maximum(G.xj))
    println("xk min = ", minimum(G.xk), "  max = ", maximum(G.xk))
    println("yi min = ", minimum(G.yi), "  max = ", maximum(G.yi))
    println("yj min = ", minimum(G.yj), "  max = ", maximum(G.yj))
    println("zk min = ", minimum(G.zk), "  max = ", maximum(G.zk))

    # Periodic enforcement: last point = first point
    G.xi[G.NI[1],:,:,:] .= G.xi[1,:,:,:]
    G.yi[G.NI[1],:,:,:] .= G.yi[1,:,:,:]
    G.zi[G.NI[1],:,:,:] .= G.zi[1,:,:,:]
    G.xj[:,G.NJ[1],:,:] .= G.xj[:,1,:,:]
    G.yj[:,G.NJ[1],:,:] .= G.yj[:,1,:,:]
    G.zj[:,G.NJ[1],:,:] .= G.zj[:,1,:,:]
    G.xk[:,:,G.NK[1],:] .= G.xk[:,:,1,:]
    G.yk[:,:,G.NK[1],:] .= G.yk[:,:,1,:]
    G.zk[:,:,G.NK[1],:] .= G.zk[:,:,1,:]

    # Step 2: Jacobian + inverse-metric terms
    @inbounds for nbl = 1:G.nblocks
        for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]

            xil = G.xi[i,j,k,nbl]; yil = G.yi[i,j,k,nbl]; zil = G.zi[i,j,k,nbl]
            xjl = G.xj[i,j,k,nbl]; yjl = G.yj[i,j,k,nbl]; zjl = G.zj[i,j,k,nbl]
            xkl = G.xk[i,j,k,nbl]; ykl = G.yk[i,j,k,nbl]; zkl = G.zk[i,j,k,nbl]

            vol = xil*(yjl*zkl - ykl*zjl) -
                  xjl*(yil*zkl - ykl*zil) +
                  xkl*(yil*zjl - zil*yjl)

            G.Jac[i,j,k,nbl] = 1.0/vol

            G.ix[i,j,k,nbl] = (yjl*zkl - ykl*zjl)/vol
            G.iy[i,j,k,nbl] = (xkl*zjl - xjl*zkl)/vol
            G.iz[i,j,k,nbl] = (xjl*ykl - xkl*yjl)/vol

            G.jx[i,j,k,nbl] = (ykl*zil - yil*zkl)/vol
            G.jy[i,j,k,nbl] = (xil*zkl - xkl*zil)/vol
            G.jz[i,j,k,nbl] = (xkl*yil - xil*ykl)/vol

            G.kx[i,j,k,nbl] = (yil*zjl - yjl*zil)/vol
            G.ky[i,j,k,nbl] = (xjl*zil - xil*zjl)/vol
            G.kz[i,j,k,nbl] = (xil*yjl - xjl*yil)/vol
        end
    end

    println("Jac min = ", minimum(G.Jac))
    println("Jac max = ", maximum(G.Jac))
    println("ix  min = ", minimum(G.ix), "  max = ", maximum(G.ix))
end

metrics!()
