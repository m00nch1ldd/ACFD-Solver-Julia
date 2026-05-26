#!/usr/bin/env julia
using DelimitedFiles
using Plots

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTDIR = joinpath(ROOT, "benchmarks")

runtime = readdlm(joinpath(OUTDIR, "runtime_vs_grid.csv"), ',', String; skipstart=1)
speedup = readdlm(joinpath(OUTDIR, "speedup_gpu_vs_cpu.csv"), ',', String; skipstart=1)

grid = parse.(Int, speedup[:,1])
cpu = parse.(Float64, speedup[:,2])
gpu = parse.(Float64, speedup[:,3])
sp  = parse.(Float64, speedup[:,4])

p1 = plot(grid, cpu, marker=:circle, lw=2, label="CPU serial", xlabel="Grid size (N in NxNxN)", ylabel="Runtime (s)", title="CPU runtime vs grid size")
savefig(p1, joinpath(OUTDIR, "cpu_runtime_vs_grid.png"))

p2 = plot(grid, gpu, marker=:diamond, lw=2, label="GPU mode", xlabel="Grid size (N in NxNxN)", ylabel="Runtime (s)", title="GPU runtime vs grid size")
savefig(p2, joinpath(OUTDIR, "gpu_runtime_vs_grid.png"))

p3 = plot(grid, sp, marker=:star5, lw=2, label="Speedup", xlabel="Grid size (N in NxNxN)", ylabel="T_CPU / T_GPU", title="GPU speedup vs grid size")
hline!(p3, [1.0], ls=:dash, c=:black, label="Parity")
savefig(p3, joinpath(OUTDIR, "gpu_speedup_vs_grid.png"))

println("Saved plots under benchmarks/")
