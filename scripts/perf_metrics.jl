#!/usr/bin/env julia
using DelimitedFiles
using Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INPUT_PATH = joinpath(ROOT, "input.dat")
const MAIN_PATH = joinpath(ROOT, "Main.jl")
const OUTDIR = joinpath(ROOT, "benchmarks")
mkpath(OUTDIR)

# grid sizes used for report-ready curves
grid_sizes = [32, 48, 64]

# serial CPU and GPU backend (if available)
modes = Dict(
    "cpu_serial" => (0, 0),
    "async_tasks" => (1, 1),
    "multithreading" => (1, 2),
    "distributed" => (1, 3),
    "gpu" => (1, 4),
)

function rewrite_input(ni::Int, nj::Int, nk::Int, exec_mode::Int, parallel_mode::Int)
    lines = readlines(INPUT_PATH)
    lines[2] = "$(ni) $(nj) $(nk)                ! NImax,NJmax,NKmax"
    lines[8] = "$(exec_mode) $(parallel_mode)                     ! exec_mode parallel_mode"
    open(INPUT_PATH, "w") do io
        for l in lines
            println(io, l)
        end
    end
end

function run_case(ni::Int, mode_name::String, exec_mode::Int, parallel_mode::Int)
    rewrite_input(ni, ni, ni, exec_mode, parallel_mode)
    cmd = mode_name == "distributed" ? `julia -p 4 $(MAIN_PATH)` : `julia $(MAIN_PATH)`
    output = read(cmd, String)
    m = match(r"Time loop ONLY\s*:\s*([0-9.]+)", output)
    isnothing(m) && error("Could not parse loop time for mode=$mode_name, N=$ni")
    return parse(Float64, m.captures[1])
end

results = String[]
push!(results, "timestamp,mode,grid_n,runtime_s")

for n in grid_sizes
    for (mode, (em, pm)) in modes
        rt = run_case(n, mode, em, pm)
        push!(results, "$(Dates.now()),$(mode),$(n),$(rt)")
    end
end

runtime_csv = joinpath(OUTDIR, "runtime_vs_grid.csv")
open(runtime_csv, "w") do io
    println(io, join(results, "\n"))
end

# speedup CSV
raw = readdlm(runtime_csv, ',', String; skipstart=1)
open(joinpath(OUTDIR, "speedup_vs_serial.csv"), "w") do io
    println(io, "grid_n,mode,serial_runtime_s,mode_runtime_s,speedup")
    for n in grid_sizes
        cpu = parse(Float64, raw[(raw[:,2].=="cpu_serial") .& (raw[:,3].==string(n)), 4][1])
        for mode in ("async_tasks", "multithreading", "distributed", "gpu")
            rt = parse(Float64, raw[(raw[:,2].==mode) .& (raw[:,3].==string(n)), 4][1])
            println(io, "$(n),$(mode),$(cpu),$(rt),$(cpu/rt)")
        end
    end
end

println("Wrote: ", runtime_csv)
println("Wrote: ", joinpath(OUTDIR, "speedup_vs_serial.csv"))
println("Strong scaling (multi-GPU) requires a cluster + CUDA-aware setup and is documented in README.")
