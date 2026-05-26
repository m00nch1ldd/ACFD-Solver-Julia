# ACFD-Solver-Julia

Serial Julia CFD solver for periodic compressible test cases:
- **Taylor–Green vortex (TGV, 3D)**
- **Co-rotating vortex / CoVo variants (2D)**

The code uses explicit RK4 time integration, central finite-difference/compact spatial derivatives, and explicit high-order filtering.

## 1) Code structure (what each file does)

- `Main.jl`  
  Driver: runs preprocessing, time loop, filtering, monitoring, and final output.
- `Preprocessing_routines.jl`  
  Reads `input.dat`, allocates arrays, generates grid, initializes flow field, sets scheme/filter coefficients, computes metrics.
- `Solver_routines.jl`  
  Core numerics: derivatives (explicit/compact), fluxes, RK update, filtering (`F2/F4/F6/F8/F10`), primitive/conservative updates.
- `Postprocessing_routines.jl`  
  Writes solver outputs (`grid.xyz`, `flow*.xyz`, `restart.xyz`) in Fortran-record binary format.
- `Module.jl`  
  Global simulation state (`G`) and shared arrays/parameters.
- `cyclictdma.jl` + `cyclictdma-call.jl`  
  Cyclic TDMA solver used by compact schemes.
- `input.dat`  
  Runtime setup file (case, grid, numerics, timestep, etc.).

## 2) `input.dat` format

Current parser expects **7 lines**:

1. `restart  nblocks`  
2. `NImax  NJmax  NKmax`  
3. `testcase  viscous`  
   - `testcase = 1` → TGV (3D)
   - `testcase = 2` → CoVo uniform grid (2D)
   - `testcase = 3` → CoVo random-perturbed grid (2D)
   - `testcase = 4` → CoVo sinusoidal grid (2D)
4. `Re  Mach  gamma  prandtl  T_ref`  
5. `nprims  nconserv` (typically `6 5`)  
6. `dscheme  fscheme  alpha_f`  
   - `dscheme`: `1=E2`, `2=E4`, `3=C4`, `4=C6`
   - `fscheme`: `2,4,6,8,10` (filter order, e.g. `10` = F10)
   - `alpha_f`: filter compactness/control parameter
7. `rk_steps  nsteps  time_step  animfreq`

## 3) How solver works (short)

1. Read input and allocate all arrays.
2. Build grid from `testcase`.
3. Initialize primitive fields (`rho, u, v, w, p, T`) and conservative fields.
4. Compute grid metrics/Jacobian.
5. Time-march for `nsteps` with RK4:
   - evaluate RHS (`unsteady!`)
   - update conservative variables
   - at final RK stage, apply filter in i/j/(k for 3D)
   - recover primitives
6. Write monitor history each step; write periodic snapshots every `animfreq`; write final files.

## 4) Output files

- `Monitor.out`  
  Time history: time, TKE, enstrophy, residuals.
- `grid.xyz`  
  Grid coordinates (binary, Fortran-record format).
- `flow%05d.xyz`  
  Intermediate flow snapshots when `iter % animfreq == 0`.
- `flow.xyz`  
  Final flow solution.
- `restart.xyz`  
  Restart file.

> Notes:
> - Files are **binary**, not plain text.
> - `*.xyz` here is solver-specific Fortran-record data (not simple XYZ text).

## 5) Running different schemes/cases

Run command:

```bash
julia Main.jl
```

Edit `input.dat` line 6 for scheme/filter:

- **E2 + F10**: `1 10 0.495`
- **E4 + F10**: `2 10 0.495`
- **C4 + F10**: `3 10 0.495`
- **C6 + F10**: `4 10 0.495`

Change case via line 3:
- TGV: `1 1` (usually 3D grid, `NKmax > 1`)
- CoVo: `2 1` or `3 1` or `4 1` (typically 2D, use `NKmax = 1`)

## 6) ParaView visualization (grid + flow)

Because outputs are custom Fortran binary, use one of these paths:

1. **Preferred**: convert `grid.xyz` + `flow*.xyz` to VTK/VTU (script/tool), then open in ParaView.
2. **Alternative**: use a ParaView reader/plugin that supports Plot3D/Fortran-record format matching this writer.

Recommended workflow:
- Keep `grid.xyz` and selected `flowXXXXX.xyz` together.
- Convert to VTK (structured grid) with variables (rho,u,v,w,p,T, etc.).
- Open resulting `.vts/.vtu` in ParaView and apply contour/slice/vector filters.

---
If you want, next step I can add a **small Julia converter** in this repo to directly write `.vts` for ParaView from `grid.xyz` + `flow.xyz`.
