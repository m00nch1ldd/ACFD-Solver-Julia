# Stand-alone verification driver for tdmap! (not used by the main solver).
# Mirrors the reference test in the Fortran cyclictdma-call.f90.

include("Module.jl")
include("cyclictdma.jl")
using .CFDVars
import .CFDVars: G

NMAX = 100

AD = zeros(NMAX)
BD = zeros(NMAX)
CD = zeros(NMAX)
DD = zeros(NMAX)

for i = 1:NMAX-1
    AD[i] = 0.2
    BD[i] = 1.0
    CD[i] = 0.5

    im1 = i - 1
    ip1 = i + 1
    if im1 < 1       im1 = NMAX-1  end
    if ip1 > NMAX-1  ip1 = 1       end

    DD[i] = CD[i]*im1 + BD[i]*i + AD[i]*ip1
end

tdmap!(1, NMAX-1, CD, BD, AD, DD)

for i = 1:NMAX-1
    println(DD[i])
end
