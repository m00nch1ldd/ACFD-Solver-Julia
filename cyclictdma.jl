# Periodic / cyclic tridiagonal solver.
# Buffers qq, ss, fei are reused from G to avoid per-call allocation.
function tdmap!(ji, jf, ap, ac, am, fi)

    ja = ji + 1
    jj = ji + jf

    qq  = G.tdma_qq
    ss  = G.tdma_ss
    fei = G.tdma_fei

    qq[ji] = -ap[ji] / ac[ji]
    ss[ji] = -am[ji] / ac[ji]
    fnn    =  fi[jf]
    fi[ji] =  fi[ji] / ac[ji]

    # forward elimination
    for j = ja:jf
        pp    = 1.0 / (ac[j] + am[j]*qq[j-1])
        qq[j] = -ap[j] * pp
        ss[j] = -am[j] * ss[j-1] * pp
        fi[j] = (fi[j] - am[j]*fi[j-1]) * pp
    end

    # backward pass
    ss[jf]  = 1.0
    fei[jf] = 0.0
    for i = ja:jf
        j      = jj - i
        ss[j]  = ss[j] + qq[j]*ss[j+1]
        fei[j] = fi[j] + qq[j]*fei[j+1]
    end

    fi[jf] = (fnn - ap[jf]*fei[ji] - am[jf]*fei[jf-1]) /
             (ap[jf]*ss[ji] + am[jf]*ss[jf-1] + ac[jf])

    # backward substitution
    for i = ja:jf
        j     = jj - i
        fi[j] = fi[jf]*ss[j] + fei[j]
    end
end
