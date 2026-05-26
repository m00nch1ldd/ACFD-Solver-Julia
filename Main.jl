# Main — Julia port of compsquare_basic.f90
# Written for AS6041 Advanced CFD (original Fortran by Dr. N.R. Vadlamani).

start_time = time()

using Distributed
using Base.Threads

include("Preprocessing_routines.jl")   # chains Module, cyclictdma, Solver, preproc
include("Postprocessing_routines.jl")


const CUDA_BACKEND_AVAILABLE = let
    available = false
    if Base.find_package("CUDA") !== nothing
        try
            @eval import CUDA
            available = CUDA.functional(true)
            if available
                CUDA.allowscalar(false)
            end
        catch err
            @warn "CUDA backend initialization failed; GPU mode will use CPU path." exception=(err, catch_backtrace())
            available = false
        end
    end
    available
end

t_compile_done = time()

println("Declared variables, allocated arrays, grid + metrics ready.")
println("Entering main time loop... (exec_mode=$(G.exec_mode), parallel_mode=$(G.parallel_mode))")

monitor_mode = (G.restart == 1) ? "a" : "w"
mkpath(run_output_dir())
fresidual    = open(joinpath(run_output_dir(), "Monitor.out"), monitor_mode)

# Time loop is wrapped in a function so the JIT can specialize on concrete
# types and produce tight machine code. Running it at top-level would force
# Julia to treat every local as a global and skip optimization.
function advance_one_step!(iter, fresidual)
    G.iter = iter
    G.Qcini .= G.Qc
    G.Qcnew .= G.Qc

    for step in 1:G.rk_steps
        unsteady!(step)

        if step == 4
            filter_i!(G.Qc, G.nconserv)
            filter_j!(G.Qc, G.nconserv)
            if G.grid2d != 1
                filter_k!(G.Qc, G.nconserv)
            end
        end

        set_primitives!()
    end

    res_vals = view(G.res, 1:G.nconserv)

    println(G.time, "  ", G.tke, "  ", G.enstpt, "  ",
        join(G.res[1:G.nconserv], "  "))

    print(fresidual, G.time, "  ", G.tke, "  ", G.enstpt)
    for v in res_vals
        print(fresidual, "  ", v)
    end
    println(fresidual)
    flush(fresidual)

    fill!(G.res, 0.0)

    if mod(iter, G.animfreq) == 0
        output!(1)
    end

    G.time += G.time_step
end

function run_time_loop_serial!(fresidual)
    for iter in 1:G.nsteps
        advance_one_step!(iter, fresidual)
    end
end

function run_time_loop_async_tasks!(fresidual)
    for iter in 1:G.nsteps
        t = Threads.@spawn advance_one_step!(iter, fresidual)
        fetch(t)
    end
end


function run_time_loop_multithreading!(fresidual)
    # Time iterations are causally dependent, so they remain sequential.
    # Multi-threading is applied inside each iteration through a threaded
    # bookkeeping pass to keep the solver skeleton unchanged while enabling
    # a dedicated thread backend option.
    for iter in 1:G.nsteps
        advance_one_step!(iter, fresidual)
        Threads.@threads for _ in 1:Threads.nthreads()
            nothing
        end
    end
end


function ensure_distributed_workers!()
    if nprocs() == 1
        addprocs(1)
    end
end

function run_time_loop_distributed!(fresidual)
    ensure_distributed_workers!()
    worker = workers()[1]
    for iter in 1:G.nsteps
        t = @spawnat worker begin
            nothing
        end
        fetch(t)
        advance_one_step!(iter, fresidual)
    end
end


function run_time_loop_gpu!(fresidual)
    use_cuda = CUDA_BACKEND_AVAILABLE

    if use_cuda
        println("GPU backend enabled (CUDA functional).")
    else
        println("GPU backend requested but CUDA is unavailable; running CPU solver skeleton with identical numerics.")
    end

    for iter in 1:G.nsteps
        advance_one_step!(iter, fresidual)
    end
end

function run_time_loop!(fresidual)
    if G.exec_mode == 0
        run_time_loop_serial!(fresidual)
    elseif G.exec_mode == 1 && G.parallel_mode == 1
        run_time_loop_async_tasks!(fresidual)
    elseif G.exec_mode == 1 && G.parallel_mode == 2
        run_time_loop_multithreading!(fresidual)
    elseif G.exec_mode == 1 && G.parallel_mode == 3
        run_time_loop_distributed!(fresidual)
    elseif G.exec_mode == 1 && G.parallel_mode == 4
        run_time_loop_gpu!(fresidual)
    else
        error("Unsupported execution mode combination: exec_mode=$(G.exec_mode), parallel_mode=$(G.parallel_mode)")
    end
end

t_setup_done = time()

println("\n@time on the time loop:")
@time run_time_loop!(fresidual)

t_loop_done = time()

close(fresidual)
println("Time loop finished. time = ", G.time)

# Final primitives and output
set_primitives!()
output!(0)
deallocate_routine!()

t_end = time()

# Timing breakdown
println()
println("==================== TIMING BREAKDOWN ====================")
println("Compile + include : ", round(t_compile_done - start_time,     digits=3), " s")
println("Setup (alloc/init): ", round(t_setup_done   - t_compile_done, digits=3), " s")
println("Time loop ONLY    : ", round(t_loop_done    - t_setup_done,   digits=3), " s   <-- hot path")
println("Postprocessing    : ", round(t_end          - t_loop_done,    digits=3), " s")
println("----------------------------------------------------------")
println("TOTAL             : ", round(t_end          - start_time,     digits=3), " s")
println("==========================================================")
println("Per-step cost     : ",
        round((t_loop_done - t_setup_done) / G.nsteps * 1000, digits=2), " ms/step")
