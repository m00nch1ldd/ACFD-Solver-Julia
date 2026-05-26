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

## Performance metrics workflow

### Generate CSVs
```bash
julia scripts/perf_metrics.jl
```
Creates:
- `benchmarks/runtime_vs_grid.csv`
- `benchmarks/speedup_gpu_vs_cpu.csv`

### Generate plots
```bash
julia scripts/plot_metrics.jl
```
Creates report-ready PNGs:
- `benchmarks/cpu_runtime_vs_grid.png`
- `benchmarks/gpu_runtime_vs_grid.png`
- `benchmarks/gpu_speedup_vs_grid.png`

## Strong scaling (multiple GPUs)
For true multi-GPU strong scaling, run fixed grid size across 1,2,4,... GPUs on a multi-GPU node/cluster and report:
- runtime vs GPU count
- parallel efficiency = T1 / (N * TN)

This repository provides the CSV/plot pipeline; cluster launcher commands are system-specific (MPI/Slurm/local CUDA_VISIBLE_DEVICES).
