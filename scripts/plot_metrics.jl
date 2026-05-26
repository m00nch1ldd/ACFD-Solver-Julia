#!/usr/bin/env julia
using DelimitedFiles
using Plots

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTDIR = joinpath(ROOT, "benchmarks")

runtime = readdlm(joinpath(OUTDIR, "runtime_vs_grid.csv"), ',', String; skipstart=1)
speedup = readdlm(joinpath(OUTDIR, "speedup_vs_serial.csv"), ',', String; skipstart=1)

grid = sort(unique(parse.(Int, runtime[:,3])))

function mode_runtime(mode)
    [parse(Float64, runtime[(runtime[:,2].==mode) .& (runtime[:,3].==string(n)), 4][1]) for n in grid]
end

cpu = mode_runtime("cpu_serial")
async = mode_runtime("async_tasks")
threads = mode_runtime("multithreading")
dist = mode_runtime("distributed")
gpu = mode_runtime("gpu")

p1 = plot(grid, cpu, marker=:circle, lw=2, label="CPU serial", xlabel="Grid size (N in NxNxN)", ylabel="Runtime (s)", title="Runtime vs grid size")
plot!(p1, grid, async, marker=:utriangle, lw=2, label="Async tasks")
plot!(p1, grid, threads, marker=:square, lw=2, label="Multithreading")
plot!(p1, grid, dist, marker=:diamond, lw=2, label="Distributed")
plot!(p1, grid, gpu, marker=:star5, lw=2, label="GPU")
savefig(p1, joinpath(OUTDIR, "runtime_all_modes_vs_grid.png"))

for mode in ("async_tasks", "multithreading", "distributed", "gpu")
    rows = speedup[speedup[:,2].==mode, :]
    g = parse.(Int, rows[:,1])
    sp = parse.(Float64, rows[:,5])
    p = plot(g, sp, marker=:circle, lw=2, label="$(mode) speedup", xlabel="Grid size (N in NxNxN)", ylabel="T_serial / T_mode", title="Speedup vs serial: $(mode)")
    hline!(p, [1.0], ls=:dash, c=:black, label="Parity")
    savefig(p, joinpath(OUTDIR, "speedup_$(mode)_vs_serial.png"))
end

println("Saved plots under benchmarks/")
