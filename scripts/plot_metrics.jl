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

function mode_speedup(mode)
    rows = speedup[speedup[:,2].==mode, :]
    g = parse.(Int, rows[:,1])
    sp = parse.(Float64, rows[:,5])
    return g, sp
end

cpu = mode_runtime("cpu_serial")
async = mode_runtime("async_tasks")
threads = mode_runtime("multithreading")
dist = mode_runtime("distributed")
gpu = mode_runtime("gpu")

# Panel 1: runtime comparison (all modes)
p_runtime = plot(grid, cpu, marker=:circle, lw=2, label="CPU serial",
    xlabel="Grid size (N in NxNxN)", ylabel="Runtime (s)", title="Runtime vs grid size")
plot!(p_runtime, grid, async, marker=:utriangle, lw=2, label="Async tasks")
plot!(p_runtime, grid, threads, marker=:square, lw=2, label="Multithreading")
plot!(p_runtime, grid, dist, marker=:diamond, lw=2, label="Distributed")
plot!(p_runtime, grid, gpu, marker=:star5, lw=2, label="GPU")

# Panels 2-5: speedup per mode vs serial
speed_modes = ["async_tasks", "multithreading", "distributed", "gpu"]

panels = Any[p_runtime]
for mode in speed_modes
    g, sp = mode_speedup(mode)
    p = plot(g, sp, marker=:circle, lw=2, label="$(mode) speedup",
        xlabel="Grid size (N in NxNxN)", ylabel="T_serial / T_mode",
        title="Speedup vs serial: $(mode)")
    hline!(p, [1.0], ls=:dash, c=:black, label="Parity")
    push!(panels, p)
end

# Single figure with subplots
p_all = plot(panels..., layout=(3,2), size=(1400,1200), legend=:best)
savefig(p_all, joinpath(OUTDIR, "all_metrics_subplots.png"))

# Keep the existing runtime figure output for backward compatibility
savefig(p_runtime, joinpath(OUTDIR, "runtime_all_modes_vs_grid.png"))

println("Saved plots under benchmarks/")
println("- ", joinpath(OUTDIR, "all_metrics_subplots.png"))
