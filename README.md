# ACFD-Solver-Julia

Serial Julia finite-difference solver for periodic compressible flow on structured grids.
Current implemented cases:
- **Taylor–Green vortex (TGV, testcase=1, 3D)**
- **CoVo family (testcase=2/3/4, 2D)**

Numerics in this version:
- Spatial derivative schemes: **E2, E4, C4, C6**
- Explicit filtering: **F2/F4/F6/F8/F10**
- Time integration: **4-stage RK (classical RK4 coefficients)**
- Periodic treatment in all active directions

## 1) Repository structure

- `Main.jl`  
  Main driver. Includes preprocessing + postprocessing, runs RK loop, writes monitor/output files.
- `Module.jl`  
  Defines global simulation state `G` (all parameters + arrays).
- `Preprocessing_routines.jl`  
  Input parsing, memory allocation, grid generation, initial condition setup, discretization/filter coefficient setup, and metric computation.
- `Solver_routines.jl`  
  Derivative operators (explicit + compact), flux calculations, RK update (`unsteady!`), filtering, primitive recovery.
- `cyclictdma.jl`, `cyclictdma-call.jl`  
  Cyclic tridiagonal solver used by compact derivative/filter solves.
- `Postprocessing_routines.jl`  
  Binary writer for `grid.xyz`, `flowxxxxx.xyz`, `flow.xyz`, `restart.xyz`.
- `input.dat`  
  Runtime controls.

## 2) `input.dat` (exact 8-line format)

1. `restart  nblocks`  
2. `NImax  NJmax  NKmax`  
3. `testcase  viscous`  
4. `Re  Mach  gamma  prandtl  T_ref`  
5. `nprims  nconserv`  
6. `dscheme  fscheme  alpha_f`  
7. `rk_steps  nsteps  time_step  animfreq`
8. `exec_mode  parallel_mode`

### Field meaning

- `restart`: restart flag used by monitor mode handling (`0` fresh, `1` append monitor).
- `nblocks`: number of blocks (current setup is generally used as 1 block).
- `NImax,NJmax,NKmax`: structured grid size.
- `testcase`:  
  - `1` = TGV, uniform periodic 3D domain (`Lx=Ly=Lz=2π`)
  - `2` = CoVo, uniform Cartesian 2D grid
  - `3` = CoVo, random-perturbed 2D grid
  - `4` = CoVo, sinusoidally distorted 2D grid
- `viscous`: viscous terms flag.
- `Re, Mach, gamma, prandtl, T_ref`: physical/non-dimensional parameters.
- `nprims,nconserv`: expected primitive/conservative counts (normally `6` and `5`).
- `dscheme`: derivative scheme selector:  
  - `1=E2` (2nd-order explicit central)
  - `2=E4` (4th-order explicit central)
  - `3=C4` (4th-order compact)
  - `4=C6` (6th-order compact)
- `fscheme`: filter order selector (`2,4,6,8,10`).
- `alpha_f`: filter compactness parameter (for F10 commonly `0.495`).
- `rk_steps`: number of RK stages (current coefficients are RK4).
- `nsteps`: total physical time steps.
- `time_step`: Δt.
- `animfreq`: write `flowxxxxx.xyz` every `animfreq` iterations.
- `exec_mode`: execution selector (`0=serial`, `1=parallel`).
- `parallel_mode`: parallel backend selector (`1=asynchronous tasks/coroutines`, `2=multi-threading`).

For this phase:
- Use `0 0` on line 8 for serial execution.
- Use `1 1` on line 8 for async-task execution.
- Use `1 2` on line 8 for multi-threading execution.

## 3) How the current solver advances one case

1. **Read input + allocate arrays** (`read_input!`, `allocate_routine!`).
2. **Generate grid** from `testcase` (`generate_grid!`).
3. **Set initial flow state** (`initialize_non_dimensionalize!`) for TGV/CoVo.
4. **Compute grid metrics + Jacobian** (`metrics!`) using selected `dscheme`.
5. **Time loop** (`run_time_loop!` in `Main.jl`):
   - save current conservative state (`Qcini`, `Qcnew`)
   - RK stage loop: `unsteady!` computes/update RHS and state
   - apply filter at final RK stage: `filter_i!`, `filter_j!`, and `filter_k!` when 3D
   - recover primitives with `set_primitives!`
   - write monitor line (`time, tke, enstrophy, residuals`)
   - periodic snapshot output when `iter % animfreq == 0`
6. **Finalization**: final primitive update + `flow.xyz`, `restart.xyz`, and `grid.xyz` output.

## 4) Output files produced

- `Monitor.out`  
  Iteration history: time, TKE, enstrophy, and residual components.
- `grid.xyz`  
  Structured grid coordinates (Fortran unformatted record binary).
- `flowxxxxx.xyz`  
  Intermediate flow snapshots at output frequency (`xxxxx` = iteration index, 5-digit).
- `flow.xyz`  
  Final flow state.
- `restart.xyz`  
  Restart state.

### Variables written in flow files

For each block, solver writes:
1. primitive variables (`nprims`): typically `rho, u, v, w, p, T`
2. conservative variables (`nconserv`): typically `rho, rho*u, rho*v, rho*w, rho*E`
3. metric terms (`ix,iy,iz,jx,jy,jz,kx,ky,kz`) in `flow*.xyz` and `flow.xyz`

> Important: These `.xyz` files are **binary Fortran-record files**, not ASCII XYZ point-cloud format.

## 5) Running cases and schemes

Run (serial):

```bash
julia Main.jl
```

Run (parallel mode 1 = async tasks/coroutines):

```bash
JULIA_NUM_THREADS=4 julia Main.jl
```

Run (parallel mode 2 = multi-threading):

```bash
JULIA_NUM_THREADS=4 julia Main.jl
```

Edit `input.dat` line 8 to pick mode, then rerun. Suggested values:
- serial: `0 0`
- async tasks: `1 1`
- multi-threading: `1 2`

### Scheme examples with F10

Use line 6 as:
- **E2 + F10**: `1 10 0.495`
- **E4 + F10**: `2 10 0.495`
- **C4 + F10**: `3 10 0.495`
- **C6 + F10**: `4 10 0.495`

### Case selection examples

- **TGV 3D**: line 3 = `1 1`, and set `NKmax > 1`.
- **CoVo 2D**: line 3 = `2 1` or `3 1` or `4 1`, and set `NKmax = 1`.

## 6) ParaView visualization workflow

Since outputs are custom Fortran-record binary:

1. Convert `grid.xyz` + chosen `flowxxxxx.xyz` (or `flow.xyz`) into VTK structured-grid format (`.vts`/`.vtu`) using your converter script/tool.
2. Open converted VTK file in ParaView.
3. Apply `Slice`, `Contour`, `Glyph`, `Calculator` as needed for velocity, pressure, density, vortical structures.

For transient visualization, convert multiple `flowxxxxx.xyz` files and load them as a time series in ParaView.


## 7) Parallel mode (current stage)

The code now supports a serial/parallel switch through `input.dat` line 8.
This stage supports `parallel_mode=1` (asynchronous tasks/coroutines) and `parallel_mode=2` (multi-threading).
Future stages can extend this field to distributed and GPU backends without changing the input structure.
