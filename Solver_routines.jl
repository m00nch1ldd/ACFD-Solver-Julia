# Explicit schemes for flow variables (second-order and fourth-order central)

function discretization_i_exp!(PHI, PHID, nvars)

    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0

    @inbounds for var = 1:nvars      
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for j   = 1:G.NJ[nbl]
        
        NI = G.NI[nbl]
        
        # interior points
        for i = 3:NI-2
            PHID[i,j,k,nbl,var] = 
                bb4*(PHI[i+2,j,k,nbl,var] - PHI[i-2,j,k,nbl,var]) +
                ab2*(PHI[i+1,j,k,nbl,var] - PHI[i-1,j,k,nbl,var])
        end
        
        # periodic boundary points
        PHID[1,j,k,nbl,var] =
            bb4*(PHI[3,j,k,nbl,var]    - PHI[NI-2,j,k,nbl,var]) +
            ab2*(PHI[2,j,k,nbl,var]    - PHI[NI-1,j,k,nbl,var])
        
        PHID[2,j,k,nbl,var] =
            bb4*(PHI[4,j,k,nbl,var]    - PHI[NI-1,j,k,nbl,var]) +
            ab2*(PHI[3,j,k,nbl,var]    - PHI[1,j,k,nbl,var])
        
        PHID[NI-1,j,k,nbl,var] =
            bb4*(PHI[2,j,k,nbl,var]    - PHI[NI-3,j,k,nbl,var]) +
            ab2*(PHI[NI,j,k,nbl,var]   - PHI[NI-2,j,k,nbl,var])
        
        PHID[NI,j,k,nbl,var] =
            bb4*(PHI[3,j,k,nbl,var]    - PHI[NI-2,j,k,nbl,var]) +
            ab2*(PHI[2,j,k,nbl,var]    - PHI[NI-1,j,k,nbl,var])
        
    end
    end
    end
    end
    
end

function discretization_j_exp!(PHI, PHID, nvars)

    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0

    @inbounds for var = 1:nvars      
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for i   = 1:G.NI[nbl]

        NJ = G.NJ[nbl]

        # interior points
        for j = 3:NJ-2
            PHID[i,j,k,nbl,var] =
                bb4*(PHI[i,j+2,k,nbl,var] - PHI[i,j-2,k,nbl,var]) +
                ab2*(PHI[i,j+1,k,nbl,var] - PHI[i,j-1,k,nbl,var])
        end
        
        # periodic boundary points
        PHID[i,1,k,nbl,var] =
            bb4*(PHI[i,3,k,nbl,var]    - PHI[i,NJ-2,k,nbl,var]) +
            ab2*(PHI[i,2,k,nbl,var]    - PHI[i,NJ-1,k,nbl,var])
        
        PHID[i,2,k,nbl,var] =
            bb4*(PHI[i,4,k,nbl,var]    - PHI[i,NJ-1,k,nbl,var]) +
            ab2*(PHI[i,3,k,nbl,var]    - PHI[i,1,k,nbl,var])
        
        PHID[i,NJ-1,k,nbl,var] =
            bb4*(PHI[i,2,k,nbl,var]    - PHI[i,NJ-3,k,nbl,var]) +
            ab2*(PHI[i,NJ,k,nbl,var]   - PHI[i,NJ-2,k,nbl,var])
        
        PHID[i,NJ,k,nbl,var] =
            bb4*(PHI[i,3,k,nbl,var]    - PHI[i,NJ-2,k,nbl,var]) +
            ab2*(PHI[i,2,k,nbl,var]    - PHI[i,NJ-1,k,nbl,var])
        
    end
    end
    end
    end
    
end

function discretization_k_exp!(PHI, PHID, nvars)

    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0

    @inbounds for var = 1:nvars      
    for nbl = 1:G.nblocks
    for j   = 1:G.NJ[nbl]
    for i   = 1:G.NI[nbl]
        
        NK = G.NK[nbl]
        
        # interior points
        for k = 3:NK-2
            PHID[i,j,k,nbl,var] =
                bb4*(PHI[i,j,k+2,nbl,var] - PHI[i,j,k-2,nbl,var]) +
                ab2*(PHI[i,j,k+1,nbl,var] - PHI[i,j,k-1,nbl,var])
        end
        
        # periodic boundary points
        PHID[i,j,1,nbl,var] =
            bb4*(PHI[i,j,3,nbl,var]    - PHI[i,j,NK-2,nbl,var]) +
            ab2*(PHI[i,j,2,nbl,var]    - PHI[i,j,NK-1,nbl,var])
        
        PHID[i,j,2,nbl,var] =
            bb4*(PHI[i,j,4,nbl,var]    - PHI[i,j,NK-1,nbl,var]) +
            ab2*(PHI[i,j,3,nbl,var]    - PHI[i,j,1,nbl,var])
        
        PHID[i,j,NK-1,nbl,var] =
            bb4*(PHI[i,j,2,nbl,var]    - PHI[i,j,NK-3,nbl,var]) +
            ab2*(PHI[i,j,NK,nbl,var]   - PHI[i,j,NK-2,nbl,var])
        
        PHID[i,j,NK,nbl,var] =
            bb4*(PHI[i,j,3,nbl,var]    - PHI[i,j,NK-2,nbl,var]) +
            ab2*(PHI[i,j,2,nbl,var]    - PHI[i,j,NK-1,nbl,var])
        
    end
    end
    end
    end
    
end

# 2D case: zero out the k-derivative.
function discretization_k2d_exp!(PHI, PHID, nvars)
    PHID .= 0.0
end


# Compact schemes for flow variables (Pade-type, solved via cyclic TDMA)

function discretization_i_comp!(PHI, PHID, nvars)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for var = 1:nvars
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for j   = 1:G.NJ[nbl]
        
        NI  = G.NI[nbl]
        nm1 = NI - 1
        RHS = zeros(NI)
        
        # interior points
        for i = 3:NI-2
            RHS[i] =
                bb4*(PHI[i+2,j,k,nbl,var] - PHI[i-2,j,k,nbl,var]) +
                ab2*(PHI[i+1,j,k,nbl,var] - PHI[i-1,j,k,nbl,var])
        end
        
        # periodic boundary points
        RHS[1] =
            bb4*(PHI[3,j,k,nbl,var]    - PHI[NI-2,j,k,nbl,var]) +
            ab2*(PHI[2,j,k,nbl,var]    - PHI[NI-1,j,k,nbl,var])
        
        RHS[2] =
            bb4*(PHI[4,j,k,nbl,var]    - PHI[NI-1,j,k,nbl,var]) +
            ab2*(PHI[3,j,k,nbl,var]    - PHI[1,j,k,nbl,var])
        
        RHS[NI-1] =
            bb4*(PHI[2,j,k,nbl,var]    - PHI[NI-3,j,k,nbl,var]) +
            ab2*(PHI[NI,j,k,nbl,var]   - PHI[NI-2,j,k,nbl,var])
        
        RHS[NI] =
            bb4*(PHI[3,j,k,nbl,var]    - PHI[NI-2,j,k,nbl,var]) +
            ab2*(PHI[2,j,k,nbl,var]    - PHI[NI-1,j,k,nbl,var])
        
        # solve tridiagonal system
        tdmap!(1, nm1, @view(G.APD[1:nm1]), @view(G.ACD[1:nm1]), @view(G.AMD[1:nm1]), @view(RHS[1:nm1]))
        
        PHID[1:nm1,j,k,nbl,var] .= RHS[1:nm1]
        PHID[NI,j,k,nbl,var]     = PHID[1,j,k,nbl,var]
        
    end
    end
    end
    end
    
end



function discretization_j_comp!(PHI, PHID, nvars)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for var = 1:nvars
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for i   = 1:G.NI[nbl]
        
        NJ  = G.NJ[nbl]
        nm1 = NJ - 1
        RHS = zeros(NJ)
        
        # interior points
        for j = 3:NJ-2
            RHS[j] =
                bb4*(PHI[i,j+2,k,nbl,var] - PHI[i,j-2,k,nbl,var]) +
                ab2*(PHI[i,j+1,k,nbl,var] - PHI[i,j-1,k,nbl,var])
        end
        
        # periodic boundary points
        RHS[1] =
            bb4*(PHI[i,3,k,nbl,var]    - PHI[i,NJ-2,k,nbl,var]) +
            ab2*(PHI[i,2,k,nbl,var]    - PHI[i,NJ-1,k,nbl,var])
        
        RHS[2] =
            bb4*(PHI[i,4,k,nbl,var]    - PHI[i,NJ-1,k,nbl,var]) +
            ab2*(PHI[i,3,k,nbl,var]    - PHI[i,1,k,nbl,var])
        
        RHS[NJ-1] =
            bb4*(PHI[i,2,k,nbl,var]    - PHI[i,NJ-3,k,nbl,var]) +
            ab2*(PHI[i,NJ,k,nbl,var]   - PHI[i,NJ-2,k,nbl,var])
        
        RHS[NJ] =
            bb4*(PHI[i,3,k,nbl,var]    - PHI[i,NJ-2,k,nbl,var]) +
            ab2*(PHI[i,2,k,nbl,var]    - PHI[i,NJ-1,k,nbl,var])
        
        # solve tridiagonal system
        tdmap!(1, nm1, @view(G.APD[1:nm1]), @view(G.ACD[1:nm1]), @view(G.AMD[1:nm1]), @view(RHS[1:nm1]))
        
        PHID[i,1:nm1,k,nbl,var] .= RHS[1:nm1]
        PHID[i,NJ,k,nbl,var]     = PHID[i,1,k,nbl,var]
        
    end
    end
    end
    end
    
end

function discretization_k_comp!(PHI, PHID, nvars)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for var = 1:nvars
    for nbl = 1:G.nblocks
    for j   = 1:G.NJ[nbl]
    for i   = 1:G.NI[nbl]
        
        NK  = G.NK[nbl]
        nm1 = NK - 1
        RHS = zeros(NK)
        
        # interior points
        for k = 3:NK-2
            RHS[k] =
                bb4*(PHI[i,j,k+2,nbl,var] - PHI[i,j,k-2,nbl,var]) +
                ab2*(PHI[i,j,k+1,nbl,var] - PHI[i,j,k-1,nbl,var])
        end
        
        # periodic boundary points
        RHS[1] =
            bb4*(PHI[i,j,3,nbl,var]    - PHI[i,j,NK-2,nbl,var]) +
            ab2*(PHI[i,j,2,nbl,var]    - PHI[i,j,NK-1,nbl,var])
        
        RHS[2] =
            bb4*(PHI[i,j,4,nbl,var]    - PHI[i,j,NK-1,nbl,var]) +
            ab2*(PHI[i,j,3,nbl,var]    - PHI[i,j,1,nbl,var])
        
        RHS[NK-1] =
            bb4*(PHI[i,j,2,nbl,var]    - PHI[i,j,NK-3,nbl,var]) +
            ab2*(PHI[i,j,NK,nbl,var]   - PHI[i,j,NK-2,nbl,var])
        
        RHS[NK] =
            bb4*(PHI[i,j,3,nbl,var]    - PHI[i,j,NK-2,nbl,var]) +
            ab2*(PHI[i,j,2,nbl,var]    - PHI[i,j,NK-1,nbl,var])
        
        # solve tridiagonal system
        tdmap!(1, nm1, @view(G.APD[1:nm1]), @view(G.ACD[1:nm1]), @view(G.AMD[1:nm1]), @view(RHS[1:nm1]))
        
        PHID[i,j,1:nm1,nbl,var] .= RHS[1:nm1]
        PHID[i,j,NK,nbl,var]     = PHID[i,j,1,nbl,var]
        
    end
    end
    end
    end
    
end


# Explicit schemes for grid derivatives (apply chain-rule for periodic boundaries)

function discretization_i_exp_grid!(PHI, PHID)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for j   = 1:G.NJ[nbl]
        
        NI = G.NI[nbl]
        
        for i = 3:NI-2
            PHID[i,j,k,nbl] =
                bb4*(PHI[i+2,j,k,nbl] - PHI[i-2,j,k,nbl]) +
                ab2*(PHI[i+1,j,k,nbl] - PHI[i-1,j,k,nbl])
        end
        
        LLx = PHI[NI,j,k,nbl] - PHI[1,j,k,nbl] 

        PHID[1,j,k,nbl] =
        bb4*(PHI[3,j,k,nbl]    - PHI[NI-2,j,k,nbl] + LLx) +
        ab2*(PHI[2,j,k,nbl]    - PHI[NI-1,j,k,nbl] + LLx)

        PHID[2,j,k,nbl] =
        bb4*(PHI[4,j,k,nbl]    - PHI[NI-1,j,k,nbl] + LLx) +
        ab2*(PHI[3,j,k,nbl]    - PHI[1,j,k,nbl])

        PHID[NI-1,j,k,nbl] =
        bb4*(PHI[2,j,k,nbl]    - PHI[NI-3,j,k,nbl] + LLx) +
        ab2*(PHI[NI,j,k,nbl]   - PHI[NI-2,j,k,nbl])

        PHID[NI,j,k,nbl] =
        bb4*(PHI[3,j,k,nbl]    - PHI[NI-2,j,k,nbl] + LLx) +
        ab2*(PHI[2,j,k,nbl]    - PHI[NI-1,j,k,nbl] + LLx)
    end
    end
    end
    
end

function discretization_j_exp_grid!(PHI, PHID)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for i   = 1:G.NI[nbl]
        
        NJ = G.NJ[nbl]
        
        for j = 3:NJ-2
            PHID[i,j,k,nbl] =
                bb4*(PHI[i,j+2,k,nbl] - PHI[i,j-2,k,nbl]) +
                ab2*(PHI[i,j+1,k,nbl] - PHI[i,j-1,k,nbl])
        end
        
        LLy = PHI[i,NJ,k,nbl] - PHI[i,1,k,nbl]

        PHID[i,1,k,nbl] =
        bb4*(PHI[i,3,k,nbl]    - PHI[i,NJ-2,k,nbl] + LLy) +
        ab2*(PHI[i,2,k,nbl]    - PHI[i,NJ-1,k,nbl] + LLy)

        PHID[i,2,k,nbl] =
        bb4*(PHI[i,4,k,nbl]    - PHI[i,NJ-1,k,nbl] + LLy) +
        ab2*(PHI[i,3,k,nbl]    - PHI[i,1,k,nbl])

        PHID[i,NJ-1,k,nbl] =
        bb4*(PHI[i,2,k,nbl]    - PHI[i,NJ-3,k,nbl] + LLy) +
        ab2*(PHI[i,NJ,k,nbl]   - PHI[i,NJ-2,k,nbl])

        PHID[i,NJ,k,nbl] =
        bb4*(PHI[i,3,k,nbl]    - PHI[i,NJ-2,k,nbl] + LLy) +
        ab2*(PHI[i,2,k,nbl]    - PHI[i,NJ-1,k,nbl] + LLy)
        
    end
    end
    end
    
end

function discretization_k_exp_grid!(PHI, PHID)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for nbl = 1:G.nblocks
    for j   = 1:G.NJ[nbl]
    for i   = 1:G.NI[nbl]
        
        NK = G.NK[nbl]
        
        for k = 3:NK-2
            PHID[i,j,k,nbl] =
                bb4*(PHI[i,j,k+2,nbl] - PHI[i,j,k-2,nbl]) +
                ab2*(PHI[i,j,k+1,nbl] - PHI[i,j,k-1,nbl])
        end
        
        LLz = PHI[i,j,NK,nbl] - PHI[i,j,1,nbl]

        PHID[i,j,1,nbl] =
        bb4*(PHI[i,j,3,nbl]    - PHI[i,j,NK-2,nbl] + LLz) +
        ab2*(PHI[i,j,2,nbl]    - PHI[i,j,NK-1,nbl] + LLz)

        PHID[i,j,2,nbl] =
        bb4*(PHI[i,j,4,nbl]    - PHI[i,j,NK-1,nbl] + LLz) +
        ab2*(PHI[i,j,3,nbl]    - PHI[i,j,1,nbl])

        PHID[i,j,NK-1,nbl] =
        bb4*(PHI[i,j,2,nbl]    - PHI[i,j,NK-3,nbl] + LLz) +
        ab2*(PHI[i,j,NK,nbl]   - PHI[i,j,NK-2,nbl])

        PHID[i,j,NK,nbl] =
        bb4*(PHI[i,j,3,nbl]    - PHI[i,j,NK-2,nbl] + LLz) +
        ab2*(PHI[i,j,2,nbl]    - PHI[i,j,NK-1,nbl] + LLz)
        
    end
    end
    end
    
end


# Compact schemes for grid derivatives

function discretization_i_comp_grid!(PHI, PHID)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for j   = 1:G.NJ[nbl]
        
        NI  = G.NI[nbl]
        nm1 = NI - 1
        RHS = zeros(NI)
        
        for i = 3:NI-2
            RHS[i] =
                bb4*(PHI[i+2,j,k,nbl] - PHI[i-2,j,k,nbl]) +
                ab2*(PHI[i+1,j,k,nbl] - PHI[i-1,j,k,nbl])
        end
        
        LLx = PHI[NI,j,k,nbl] - PHI[1,j,k,nbl] 

        RHS[1] =
        bb4*(PHI[3,j,k,nbl]    - PHI[NI-2,j,k,nbl] + LLx) +
        ab2*(PHI[2,j,k,nbl]    - PHI[NI-1,j,k,nbl] + LLx)

        RHS[2] =
        bb4*(PHI[4,j,k,nbl]    - PHI[NI-1,j,k,nbl] + LLx) +
        ab2*(PHI[3,j,k,nbl]    - PHI[1,j,k,nbl])

        RHS[NI-1] =
        bb4*(PHI[2,j,k,nbl]    - PHI[NI-3,j,k,nbl] + LLx) +
        ab2*(PHI[NI,j,k,nbl]   - PHI[NI-2,j,k,nbl])

        RHS[NI] =
        bb4*(PHI[3,j,k,nbl]    - PHI[NI-2,j,k,nbl] + LLx) +
        ab2*(PHI[2,j,k,nbl]    - PHI[NI-1,j,k,nbl] + LLx)
        
        tdmap!(1, nm1, @view(G.APD[1:nm1]), @view(G.ACD[1:nm1]), @view(G.AMD[1:nm1]), @view(RHS[1:nm1]))
        
        PHID[1:nm1,j,k,nbl] .= RHS[1:nm1]
        PHID[NI,j,k,nbl]     = PHID[1,j,k,nbl]
        
    end
    end
    end
    
end

function discretization_j_comp_grid!(PHI, PHID)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for i   = 1:G.NI[nbl]
        
        NJ  = G.NJ[nbl]
        nm1 = NJ - 1
        RHS = zeros(NJ)
        
        for j = 3:NJ-2
            RHS[j] =
                bb4*(PHI[i,j+2,k,nbl] - PHI[i,j-2,k,nbl]) +
                ab2*(PHI[i,j+1,k,nbl] - PHI[i,j-1,k,nbl])
        end
        
        LLy = PHI[i,NJ,k,nbl] - PHI[i,1,k,nbl]

        RHS[1] =
        bb4*(PHI[i,3,k,nbl]    - PHI[i,NJ-2,k,nbl] + LLy) +
        ab2*(PHI[i,2,k,nbl]    - PHI[i,NJ-1,k,nbl] + LLy)

        RHS[2] =
        bb4*(PHI[i,4,k,nbl]    - PHI[i,NJ-1,k,nbl] + LLy) +
        ab2*(PHI[i,3,k,nbl]    - PHI[i,1,k,nbl])

        RHS[NJ-1] =
        bb4*(PHI[i,2,k,nbl]    - PHI[i,NJ-3,k,nbl] + LLy) +
        ab2*(PHI[i,NJ,k,nbl]   - PHI[i,NJ-2,k,nbl])

        RHS[NJ] =
        bb4*(PHI[i,3,k,nbl]    - PHI[i,NJ-2,k,nbl] + LLy) +
        ab2*(PHI[i,2,k,nbl]    - PHI[i,NJ-1,k,nbl] + LLy)
        
        tdmap!(1, nm1, @view(G.APD[1:nm1]), @view(G.ACD[1:nm1]), @view(G.AMD[1:nm1]), @view(RHS[1:nm1]))
        
        PHID[i,1:nm1,k,nbl] .= RHS[1:nm1]
        PHID[i,NJ,k,nbl]     = PHID[i,1,k,nbl]
        
    end
    end
    end
    
end

function discretization_k_comp_grid!(PHI, PHID)
    
    bb4 = G.bdisc/4.0
    ab2 = G.adisc/2.0
    
    for nbl = 1:G.nblocks
    for j   = 1:G.NJ[nbl]
    for i   = 1:G.NI[nbl]
        
        NK  = G.NK[nbl]
        nm1 = NK - 1
        RHS = zeros(NK)
        
        for k = 3:NK-2
            RHS[k] =
                bb4*(PHI[i,j,k+2,nbl] - PHI[i,j,k-2,nbl]) +
                ab2*(PHI[i,j,k+1,nbl] - PHI[i,j,k-1,nbl])
        end
        
       LLz = PHI[i,j,NK,nbl] - PHI[i,j,1,nbl]

        RHS[1] =
        bb4*(PHI[i,j,3,nbl]    - PHI[i,j,NK-2,nbl] + LLz) +
        ab2*(PHI[i,j,2,nbl]    - PHI[i,j,NK-1,nbl] + LLz)

        RHS[2] =
        bb4*(PHI[i,j,4,nbl]    - PHI[i,j,NK-1,nbl] + LLz) +
        ab2*(PHI[i,j,3,nbl]    - PHI[i,j,1,nbl])

        RHS[NK-1] =
        bb4*(PHI[i,j,2,nbl]    - PHI[i,j,NK-3,nbl] + LLz) +
        ab2*(PHI[i,j,NK,nbl]   - PHI[i,j,NK-2,nbl])

        RHS[NK] =
        bb4*(PHI[i,j,3,nbl]    - PHI[i,j,NK-2,nbl] + LLz) +
        ab2*(PHI[i,j,2,nbl]    - PHI[i,j,NK-1,nbl] + LLz)
        
        tdmap!(1, nm1, @view(G.APD[1:nm1]), @view(G.ACD[1:nm1]), @view(G.AMD[1:nm1]), @view(RHS[1:nm1]))

        PHID[i,j,1:nm1,nbl] .= RHS[1:nm1]
        PHID[i,j,NK,nbl]     = PHID[i,j,1,nbl]

    end
    end
    end

end


# Periodic filtering routines (Pade-type filters, applied at final RK stage)

function filter_i!(PHI, nvars)

    M = div(G.fscheme, 2)              # Fortran:  fscheme/2 (integer)
    RHS = G.filter_rhs                 

    @inbounds for var = 1:nvars
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for j   = 1:G.NJ[nbl]

        Ni  = G.NI[nbl]
        nm1 = Ni - 1

        for i = 1:Ni
            acc = 0.0
            for c = 0:M
                ipc = i + c
                imc = i - c
                if ipc > Ni
                    ipc -= Ni - 1
                end
                if imc < 1
                    imc += Ni - 1
                end
                acc += 0.5 * G.fcoeff[c+1] *
                       (PHI[ipc,j,k,nbl,var] + PHI[imc,j,k,nbl,var])
            end
            RHS[i] = acc
        end

        tdmap!(1, nm1, @view(G.APF[1:nm1]), @view(G.ACF[1:nm1]), @view(G.AMF[1:nm1]), @view(RHS[1:nm1]))

        #element-wise copy without slice allocation
        for i = 1:nm1
            PHI[i,j,k,nbl,var] = RHS[i]
        end
        PHI[Ni,j,k,nbl,var] = PHI[1,j,k,nbl,var]

    end
    end
    end
    end
end

function filter_j!(PHI, nvars)

    M = div(G.fscheme, 2)
    RHS = G.filter_rhs                 

    @inbounds for var = 1:nvars
    for nbl = 1:G.nblocks
    for k   = 1:G.NK[nbl]
    for i   = 1:G.NI[nbl]

        Nj  = G.NJ[nbl]
        nm1 = Nj - 1

        for j = 1:Nj
            acc = 0.0
            for c = 0:M
                jpc = j + c
                jmc = j - c
                if jpc > Nj
                    jpc -= Nj - 1
                end
                if jmc < 1
                    jmc += Nj - 1
                end
                acc += 0.5 * G.fcoeff[c+1] *
                       (PHI[i,jpc,k,nbl,var] + PHI[i,jmc,k,nbl,var])
            end
            RHS[j] = acc
        end

        tdmap!(1, nm1, @view(G.APF[1:nm1]), @view(G.ACF[1:nm1]), @view(G.AMF[1:nm1]), @view(RHS[1:nm1]))

        #element-wise copy without slice allocation
        for j = 1:nm1
            PHI[i,j,k,nbl,var] = RHS[j]
        end
        PHI[i,Nj,k,nbl,var] = PHI[i,1,k,nbl,var]

    end
    end
    end
    end
end

function filter_k!(PHI, nvars)

    M = div(G.fscheme, 2)
    RHS = G.filter_rhs                 

    @inbounds for var = 1:nvars
    for nbl = 1:G.nblocks
    for j   = 1:G.NJ[nbl]
    for i   = 1:G.NI[nbl]

        Nk  = G.NK[nbl]
        nm1 = Nk - 1

        for k = 1:Nk
            acc = 0.0
            for c = 0:M
                kpc = k + c
                kmc = k - c
                if kpc > Nk
                    kpc -= Nk - 1
                end
                if kmc < 1
                    kmc += Nk - 1
                end
                acc += 0.5 * G.fcoeff[c+1] *
                       (PHI[i,j,kpc,nbl,var] + PHI[i,j,kmc,nbl,var])
            end
            RHS[k] = acc
        end

        tdmap!(1, nm1, @view(G.APF[1:nm1]), @view(G.ACF[1:nm1]), @view(G.AMF[1:nm1]), @view(RHS[1:nm1]))

        #element-wise copy without slice allocation
        for k = 1:nm1
            PHI[i,j,k,nbl,var] = RHS[k]
        end
        PHI[i,j,Nk,nbl,var] = PHI[i,j,1,nbl,var]

    end
    end
    end
    end
end

# Convert conservative variables Qc back to primitives Qp.

function set_primitives!()

    γ    = G.gamma
    Mach = G.Mach

    @inbounds for nbl = 1:G.nblocks
        for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]

            rhl = G.Qc[i,j,k,nbl,1]
            ul  = G.Qc[i,j,k,nbl,2] / rhl
            vl  = G.Qc[i,j,k,nbl,3] / rhl
            wl  = G.Qc[i,j,k,nbl,4] / rhl
            El  = G.Qc[i,j,k,nbl,5] / rhl
            Tl  = (El - 0.5*(ul^2 + vl^2 + wl^2)) * (γ*(γ-1.0)*Mach^2)
            pl  = rhl*Tl / (γ*Mach^2)

            G.Qp[i,j,k,nbl,1] = rhl
            G.Qp[i,j,k,nbl,2] = ul
            G.Qp[i,j,k,nbl,3] = vl
            G.Qp[i,j,k,nbl,4] = wl
            G.Qp[i,j,k,nbl,5] = pl
            G.Qp[i,j,k,nbl,6] = Tl
        end
    end
end

# Main solver step: build flux divergence and advance one RK sub-stage.

function unsteady!(stepl)

    γ    = G.gamma
    Mach = G.Mach
    Re   = G.Re
    Pr   = G.prandtl
    Tref = G.T_ref

    # ---- (1) Inviscid fluxes ---------------------------------------
    @inbounds for nbl = 1:G.nblocks
        for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]

            ixl = G.ix[i,j,k,nbl]; iyl = G.iy[i,j,k,nbl]; izl = G.iz[i,j,k,nbl]
            jxl = G.jx[i,j,k,nbl]; jyl = G.jy[i,j,k,nbl]; jzl = G.jz[i,j,k,nbl]
            kxl = G.kx[i,j,k,nbl]; kyl = G.ky[i,j,k,nbl]; kzl = G.kz[i,j,k,nbl]
            vol = 1.0 / G.Jac[i,j,k,nbl]

            rhl = G.Qp[i,j,k,nbl,1]
            ul  = G.Qp[i,j,k,nbl,2]
            vl  = G.Qp[i,j,k,nbl,3]
            wl  = G.Qp[i,j,k,nbl,4]
            pl  = G.Qp[i,j,k,nbl,5]
            Tl  = G.Qp[i,j,k,nbl,6]
            El  = Tl/(γ*(γ-1.0)*Mach^2) + 0.5*(ul^2 + vl^2 + wl^2)

            Ucont = ixl*ul + iyl*vl + izl*wl
            Vcont = jxl*ul + jyl*vl + jzl*wl
            Wcont = kxl*ul + kyl*vl + kzl*wl

            G.Fflux[i,j,k,nbl,1] = -rhl*Ucont*vol
            G.Fflux[i,j,k,nbl,2] = -(rhl*ul*Ucont + ixl*pl)*vol
            G.Fflux[i,j,k,nbl,3] = -(rhl*vl*Ucont + iyl*pl)*vol
            G.Fflux[i,j,k,nbl,4] = -(rhl*wl*Ucont + izl*pl)*vol
            G.Fflux[i,j,k,nbl,5] = -(rhl*El*Ucont + pl*Ucont)*vol

            G.Gflux[i,j,k,nbl,1] = -rhl*Vcont*vol
            G.Gflux[i,j,k,nbl,2] = -(rhl*ul*Vcont + jxl*pl)*vol
            G.Gflux[i,j,k,nbl,3] = -(rhl*vl*Vcont + jyl*pl)*vol
            G.Gflux[i,j,k,nbl,4] = -(rhl*wl*Vcont + jzl*pl)*vol
            G.Gflux[i,j,k,nbl,5] = -(rhl*El*Vcont + pl*Vcont)*vol

            G.Hflux[i,j,k,nbl,1] = -rhl*Wcont*vol
            G.Hflux[i,j,k,nbl,2] = -(rhl*ul*Wcont + kxl*pl)*vol
            G.Hflux[i,j,k,nbl,3] = -(rhl*vl*Wcont + kyl*pl)*vol
            G.Hflux[i,j,k,nbl,4] = -(rhl*wl*Wcont + kzl*pl)*vol
            G.Hflux[i,j,k,nbl,5] = -(rhl*El*Wcont + pl*Wcont)*vol
        end
    end

    # ---- (2) Viscous fluxes (only if viscous == 1) -----------------
    if G.viscous == 1

        # primitive variable derivatives in i,j,k via the chosen scheme
        if G.dscheme == 1 || G.dscheme == 2
            discretization_i_exp!(G.Qp, G.Qpi, G.nprims)
            discretization_j_exp!(G.Qp, G.Qpj, G.nprims)
            if G.grid2d != 1
                discretization_k_exp!(G.Qp, G.Qpk, G.nprims)
            else
                discretization_k2d_exp!(G.Qp, G.Qpk, G.nprims)
            end
        elseif G.dscheme == 3 || G.dscheme == 4
            discretization_i_comp!(G.Qp, G.Qpi, G.nprims)
            discretization_j_comp!(G.Qp, G.Qpj, G.nprims)
            if G.grid2d != 1
                discretization_k_comp!(G.Qp, G.Qpk, G.nprims)
            else
                discretization_k2d_exp!(G.Qp, G.Qpk, G.nprims)
            end
        end

        @inbounds for nbl = 1:G.nblocks
            for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]

                ixl = G.ix[i,j,k,nbl]; iyl = G.iy[i,j,k,nbl]; izl = G.iz[i,j,k,nbl]
                jxl = G.jx[i,j,k,nbl]; jyl = G.jy[i,j,k,nbl]; jzl = G.jz[i,j,k,nbl]
                kxl = G.kx[i,j,k,nbl]; kyl = G.ky[i,j,k,nbl]; kzl = G.kz[i,j,k,nbl]
                vol = 1.0 / G.Jac[i,j,k,nbl]

                rhl = G.Qp[i,j,k,nbl,1]
                ul  = G.Qp[i,j,k,nbl,2]
                vl  = G.Qp[i,j,k,nbl,3]
                wl  = G.Qp[i,j,k,nbl,4]
                Tl  = G.Qp[i,j,k,nbl,6]

                # primitive derivatives in computational coords
                uil = G.Qpi[i,j,k,nbl,2]; vil = G.Qpi[i,j,k,nbl,3]
                wil = G.Qpi[i,j,k,nbl,4]; Til = G.Qpi[i,j,k,nbl,6]

                ujl = G.Qpj[i,j,k,nbl,2]; vjl = G.Qpj[i,j,k,nbl,3]
                wjl = G.Qpj[i,j,k,nbl,4]; Tjl = G.Qpj[i,j,k,nbl,6]

                ukl = G.Qpk[i,j,k,nbl,2]; vkl = G.Qpk[i,j,k,nbl,3]
                wkl = G.Qpk[i,j,k,nbl,4]; Tkl = G.Qpk[i,j,k,nbl,6]

                # chain rule:  d/dx = ix*d/di + jx*d/dj + kx*d/dk  (etc.)
                u_x = ixl*uil + jxl*ujl + kxl*ukl
                v_x = ixl*vil + jxl*vjl + kxl*vkl
                w_x = ixl*wil + jxl*wjl + kxl*wkl
                T_x = ixl*Til + jxl*Tjl + kxl*Tkl

                u_y = iyl*uil + jyl*ujl + kyl*ukl
                v_y = iyl*vil + jyl*vjl + kyl*vkl
                w_y = iyl*wil + jyl*wjl + kyl*wkl
                T_y = iyl*Til + jyl*Tjl + kyl*Tkl

                u_z = izl*uil + jzl*ujl + kzl*ukl
                v_z = izl*vil + jzl*vjl + kzl*vkl
                w_z = izl*wil + jzl*wjl + kzl*wkl
                T_z = izl*Til + jzl*Tjl + kzl*Tkl

                div2b3 = (2.0/3.0)*(u_x + v_y + w_z)

                # Sutherland's law (non-dimensional)
                mul = (Tl^1.5)*(1.0 + 110.4/Tref) / (Tl + 110.4/Tref)

                Txx = (2.0*u_x - div2b3)*mul/Re
                Tyy = (2.0*v_y - div2b3)*mul/Re
                Tzz = (2.0*w_z - div2b3)*mul/Re
                Txy = (u_y + v_x)*mul/Re
                Txz = (w_x + u_z)*mul/Re
                Tyz = (v_z + w_y)*mul/Re

                facprM = 1.0/((γ - 1.0)*Pr*Mach^2)

                bx = ul*Txx + vl*Txy + wl*Txz + (mul/Re)*facprM*T_x
                by = ul*Txy + vl*Tyy + wl*Tyz + (mul/Re)*facprM*T_y
                bz = ul*Txz + vl*Tyz + wl*Tzz + (mul/Re)*facprM*T_z

                # add viscous part to the inviscid flux (continuity unchanged)
                G.Fflux[i,j,k,nbl,2] += (ixl*Txx + iyl*Txy + izl*Txz)*vol
                G.Fflux[i,j,k,nbl,3] += (ixl*Txy + iyl*Tyy + izl*Tyz)*vol
                G.Fflux[i,j,k,nbl,4] += (ixl*Txz + iyl*Tyz + izl*Tzz)*vol
                G.Fflux[i,j,k,nbl,5] += (ixl*bx  + iyl*by  + izl*bz )*vol

                G.Gflux[i,j,k,nbl,2] += (jxl*Txx + jyl*Txy + jzl*Txz)*vol
                G.Gflux[i,j,k,nbl,3] += (jxl*Txy + jyl*Tyy + jzl*Tyz)*vol
                G.Gflux[i,j,k,nbl,4] += (jxl*Txz + jyl*Tyz + jzl*Tzz)*vol
                G.Gflux[i,j,k,nbl,5] += (jxl*bx  + jyl*by  + jzl*bz )*vol

                G.Hflux[i,j,k,nbl,2] += (kxl*Txx + kyl*Txy + kzl*Txz)*vol
                G.Hflux[i,j,k,nbl,3] += (kxl*Txy + kyl*Tyy + kzl*Tyz)*vol
                G.Hflux[i,j,k,nbl,4] += (kxl*Txz + kyl*Tyz + kzl*Tzz)*vol
                G.Hflux[i,j,k,nbl,5] += (kxl*bx  + kyl*by  + kzl*bz )*vol

                G.tked[i,j,k,nbl] = 0.5*rhl*(ul^2 + vl^2 + wl^2)
                G.enst[i,j,k,nbl] = 0.5*(2.0*mul/Re)*rhl*
                                    ((w_y - v_z)^2 + (w_x - u_z)^2 + (v_x - u_y)^2)
            end
        end
    end

    # ---- (3) Net flux = divergence of (Fflux, Gflux, Hflux) --------
    fill!(G.net_flux, 0.0)
    if G.dscheme == 1 || G.dscheme == 2
        discretization_i_exp!(G.Fflux, G.fluxD, G.nconserv); G.net_flux .+= G.fluxD
        discretization_j_exp!(G.Gflux, G.fluxD, G.nconserv); G.net_flux .+= G.fluxD
        if G.grid2d != 1
            discretization_k_exp!(G.Hflux, G.fluxD, G.nconserv)
        else
            discretization_k2d_exp!(G.Hflux, G.fluxD, G.nconserv)
        end
        G.net_flux .+= G.fluxD
    elseif G.dscheme == 3 || G.dscheme == 4
        discretization_i_comp!(G.Fflux, G.fluxD, G.nconserv); G.net_flux .+= G.fluxD
        discretization_j_comp!(G.Gflux, G.fluxD, G.nconserv); G.net_flux .+= G.fluxD
        if G.grid2d != 1
            discretization_k_comp!(G.Hflux, G.fluxD, G.nconserv)
        else
            discretization_k2d_exp!(G.Hflux, G.fluxD, G.nconserv)
        end
        G.net_flux .+= G.fluxD
    end

    # ---- (4) Reduce  tke, enstpt  over interior 1..N-1 -------------
    
    tke_acc    = 0.0
    enstpt_acc = 0.0
    @inbounds for nbl = 1:G.nblocks
        for k = 1:G.NK[nbl]-1, j = 1:G.NJ[nbl]-1, i = 1:G.NI[nbl]-1
            tke_acc    += G.tked[i,j,k,nbl]
            enstpt_acc += G.enst[i,j,k,nbl]
        end
        denom       = (G.NI[nbl]-1)*(G.NJ[nbl]-1)*(G.NK[nbl]-1)
        tke_acc    /= denom
        enstpt_acc /= denom
    end
    G.tke    = tke_acc
    G.enstpt = enstpt_acc

   # ---- (5) RK update ---------------------------------------------
    
    @inbounds for var = 1:G.nconserv
        for nbl = 1:G.nblocks
            for k = 1:G.NK[nbl], j = 1:G.NJ[nbl], i = 1:G.NI[nbl]

                vol = 1.0 / G.Jac[i,j,k,nbl]

                G.Qcnew[i,j,k,nbl,var] += G.time_step * G.net_flux[i,j,k,nbl,var] *
                                          G.fac_RK[stepl] / vol

                if stepl <= G.rk_steps - 1
                    G.Qc[i,j,k,nbl,var] = G.Qcini[i,j,k,nbl,var] +
                                          G.time_step * G.net_flux[i,j,k,nbl,var] *
                                          G.fac_qini[stepl+1] / vol
                else
                    G.res[var] = max(G.res[var],
                                     abs(G.Qcnew[i,j,k,nbl,var] -
                                         G.Qcini[i,j,k,nbl,var]))
                    G.Qc[i,j,k,nbl,var] = G.Qcnew[i,j,k,nbl,var]
                end
            end
        end
    end
end