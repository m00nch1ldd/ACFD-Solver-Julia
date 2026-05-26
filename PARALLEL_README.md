# Parallel execution guide

## `input.dat` control line
Line 8 selects execution mode:
- `0 0` -> Serial CPU
- `1 1` -> Async tasks/coroutines
- `1 2` -> Multi-threading
- `1 3` -> Distributed computing
- `1 4` -> GPU computing

For 64x64x64 TGV, use:
```
64 64 64
1 1
...
1 4
```

## How each parallel mode is implemented in this solver

All modes share the same numerical kernel (`advance_one_step!`), which performs RK stages, filtering, primitive recovery, monitor writes, and output frequency checks. The mode selection is done in `run_time_loop!` using `exec_mode` and `parallel_mode`.

### 1) `parallel_mode = 1` (Async tasks / coroutines)
- Implemented in `run_time_loop_async_tasks!`.
- For each time iteration, the solver launches `advance_one_step!` with `Threads.@spawn`, then immediately `fetch`es the task.
- Because each step is fetched before the next begins, iterations remain causally sequential (no overlapping timesteps), while using Julia task scheduling machinery.
- Solver areas modified/used:
  - Time-loop wrapper in `Main.jl` (`run_time_loop_async_tasks!`).
  - No change to stencil/RK/filter operators themselves.

### 2) `parallel_mode = 2` (Multi-threading)
- Implemented in `run_time_loop_multithreading!`.
- The physical timestep loop remains sequential (`advance_one_step!`), since timestep `n+1` depends on timestep `n`.
- This mode executes the same solver path with kernel-level threading now active inside `Solver_routines.jl` (primitive recovery + inviscid/viscous flux loops parallelized over `k`).
- Solver areas modified/used:
  - Time-loop wrapper in `Main.jl` (`run_time_loop_multithreading!`).
  - Core numerical routines are parallelized in `Solver_routines.jl`.

### 3) `parallel_mode = 3` (Distributed computing)
- Implemented in `run_time_loop_distributed!` plus `ensure_distributed_workers!`.
- If the run started with one process, `addprocs(1)` adds a worker.
- Each iteration spawns a small remote task (`@spawnat worker ...`) and fetches it, then runs local `advance_one_step!`.
- This currently acts as a distributed backend skeleton/integration path; the domain decomposition and RHS/filter kernels are not yet distributed across workers.
- Solver areas modified/used:
  - Distributed process management and loop wrapper in `Main.jl`.
  - No distributed decomposition changes in `Solver_routines.jl` yet.

### 4) `parallel_mode = 4` (GPU computing)
- Implemented in `run_time_loop_gpu!` with CUDA availability detection (`CUDA_BACKEND_AVAILABLE`).
- On CUDA-capable systems, each iteration currently demonstrates device path plumbing by copying `G.Qc` to `CuArray` and back, then executing `advance_one_step!`.
- If CUDA is unavailable, the code prints a fallback message and runs the same CPU solver path.
- This is a GPU backend skeleton path; core derivative/filter/RK kernels are not yet rewritten as persistent device kernels.
- Solver areas modified/used:
  - CUDA backend detection and GPU loop wrapper in `Main.jl`.
  - Numerical kernels remain in existing CPU routines.


### Kernel-level parallelization now applied
- Core solver loops in `Solver_routines.jl` are parallelized over the `k` dimension using `Threads.@threads` in:
  - `set_primitives!`
  - inviscid flux assembly inside `unsteady!`
  - viscous flux assembly inside `unsteady!`
- This means the full per-step physics update path (primitive recovery + flux construction) is now threaded rather than only backend wrappers.

## Run commands

### 1) Serial
```bash
julia Main.jl
```

### 2) Async tasks/coroutines
```bash
JULIA_NUM_THREADS=4 julia Main.jl
```

### 3) Multi-threading
```bash
JULIA_NUM_THREADS=8 julia Main.jl
```

### 4) Distributed
```bash
julia -p 4 Main.jl
```

### 5) GPU
Requires CUDA-capable GPU and CUDA.jl.
```bash
julia --project -e 'using Pkg; Pkg.add("CUDA"); Pkg.add("Plots")'
JULIA_NUM_THREADS=4 julia Main.jl
```


Quick helper (runs all 5 execution modes once and writes a status summary):
```bash
julia scripts/run_all_modes.jl
```
Creates: `benchmarks/run_all_modes_summary.csv`.

## Performance metrics workflow

### Generate CSVs
```bash
julia scripts/perf_metrics.jl
```
Creates:
- `benchmarks/runtime_vs_grid.csv`
- `benchmarks/speedup_vs_serial.csv` (includes async, threads, distributed, GPU against serial baseline)

### Generate plots
```bash
julia scripts/plot_metrics.jl
```
Creates report-ready PNGs:
- `benchmarks/all_metrics_subplots.png` (single figure with runtime + all speedup panels)
- `benchmarks/runtime_all_modes_vs_grid.png` (runtime panel, also saved separately)

## Strong scaling (multiple GPUs)
For true multi-GPU strong scaling, run fixed grid size across 1,2,4,... GPUs on a multi-GPU node/cluster and report:
- runtime vs GPU count
- parallel efficiency = T1 / (N * TN)

This repository provides the CSV/plot pipeline; cluster launcher commands are system-specific (MPI/Slurm/local CUDA_VISIBLE_DEVICES).
