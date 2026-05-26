#!/usr/bin/env julia
using Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INPUT_PATH = joinpath(ROOT, "input.dat")
const MAIN_PATH = joinpath(ROOT, "Main.jl")
const OUTDIR = joinpath(ROOT, "benchmarks")
mkpath(OUTDIR)

const MODES = [
    ("cpu_serial", 0, 0, `julia $(MAIN_PATH)`),
    ("async_tasks", 1, 1, `julia $(MAIN_PATH)`),
    ("multithreading", 1, 2, `julia $(MAIN_PATH)`),
    ("distributed", 1, 3, `julia -p 4 $(MAIN_PATH)`),
    ("gpu", 1, 4, `julia $(MAIN_PATH)`),
]

function rewrite_input(exec_mode::Int, parallel_mode::Int)
    lines = readlines(INPUT_PATH)
    lines[8] = "$(exec_mode) $(parallel_mode)                     ! exec_mode parallel_mode"
    open(INPUT_PATH, "w") do io
        for l in lines
            println(io, l)
        end
    end
end

function parse_loop_time(output::String)
    m = match(r"Time loop ONLY\s*:\s*([0-9.]+)", output)
    isnothing(m) && return missing
    return parse(Float64, m.captures[1])
end

rows = String["timestamp,mode,exec_mode,parallel_mode,status,runtime_s"]
println("Running all 5 modes sequentially...")

for (mode, em, pm, cmd) in MODES
    println("\n--- Mode: $mode (exec_mode=$em, parallel_mode=$pm) ---")
    rewrite_input(em, pm)
    t0 = time()
    status = "ok"
    runtime = missing
    try
        output = read(cmd, String)
        runtime = parse_loop_time(output)
        if ismissing(runtime)
            status = "no_time_parse"
        end
    catch err
        status = "failed"
        @warn "Run failed" mode exception=(err, catch_backtrace())
    end
    wall = round(time() - t0, digits=3)
    rt = ismissing(runtime) ? "" : string(runtime)
    push!(rows, "$(Dates.now()),$(mode),$(em),$(pm),$(status),$(rt)")
    println("status=$status, wall_s=$wall, loop_runtime_s=$(ismissing(runtime) ? "NA" : runtime)")
end

outcsv = joinpath(OUTDIR, "run_all_modes_summary.csv")
open(outcsv, "w") do io
    println(io, join(rows, "\n"))
end
println("\nWrote summary: $outcsv")
println("Tip: run `julia scripts/perf_metrics.jl` next for report metrics.")
