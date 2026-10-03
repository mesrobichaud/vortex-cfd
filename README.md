# CFD

OpenFOAM simulations of the CPB device. Meshes are built with gmsh and cases are
set up with PyFoam.

Each simulation run is a **case**, generated from a **case template** and a
set of **parameters**. The template is an ordinary OpenFOAM case in which a
few values are placeholders. The mesh, fields, physics and numerics are
separate files that you can edit directly or control through parameters.

## Layout

```
CFD/
  pyproject.toml              Python dependencies; installs src/cfdtools
  src/cfdtools/               Python code shared by all base cases
    mesh.py                   gmsh cell-count and sizing helpers
    units.py                  unit conversion to SI (rpm, mm, mPa.s, nu = mu/rho)
    monitor.py                live plots of a running case
    post.py                   read results and plot fields and profiles (matplotlib)
    pv.py                     open a case in ParaView with a publication view (python -m cfdtools.pv)
    paraview_view.py          the ParaView view used by pv.py (runs inside ParaView)

  geometry/                   gmsh mesh builders, shared by all base cases
    common.py                 device dimensions used by the builders
    annulus.py                2D r-theta annulus  (annulus_stream-plane)
    wedge.py                  axisymmetric r-z wedge  (annulus_span-plane)

  annulus_stream-plane/       base case: annular gap in the r-theta plane
    base/annulus-2D/          case template
    studies/                  parameter files for named runs and sweeps
    post/                     figure scripts and example notebook for this base case
    cases/                    generated cases (not tracked by git)

  annulus_span-plane/         base case: annular gap in the r-z plane
    base/wedge-2D/
    studies/
    cases/
```

Each `annulus_*` folder is a **base case**: one flow setup and its case
template. Usually `base/` holds a single template. Mesh builders are kept in
the top-level `geometry/` folder so that base cases modelling the same device
share one set of dimensions.

| Folder | Contents | Edit it to change |
|---|---|---|
| `src/cfdtools/` | Code that does not depend on a specific geometry | Shared tools |
| `geometry/` | Mesh builders and device dimensions | The geometry or mesh |
| `<base case>/base/<template>/` | Boundary conditions, models, schemes, solver settings, outputs | The physics or numerics |
| `<base case>/studies/` | Parameter files | The values to run |
| `<base case>/post/` | Figure scripts | The figures |
| `<base case>/cases/` | Generated cases | Nothing: regenerate instead |

Anything you want to keep goes in `base/` or `studies/`. Cases can be
deleted and regenerated at any time.

## Requirements

- **OpenFOAM v2606** (openfoam.com release), installed separately:
  - macOS: [openfoam-app](https://github.com/gerlero/openfoam-app), which
    provides the `openfoam` command used below.
  - Linux / WSL: [openfoam.com/download](https://www.openfoam.com/download).

  Other recent openfoam.com versions will probably work but are untested.
  The openfoam.org releases use different dictionary names and will not run
  these templates unchanged.
- **Python 3.10+**, in any environment (conda, venv, ...).

## Setup

From the repository root, in a new environment:

```bash
conda create -n cfd python=3.12 && conda activate cfd    # or: python -m venv .venv && source .venv/bin/activate
pip install -e .
```

This installs gmsh, PyFoam, matplotlib and `cfdtools`. `cfdtools` is
installed in editable mode, so changes under `src/` apply without
reinstalling. Run `pip install -e .` again if `pyproject.toml` changes.

PyFoam runs OpenFOAM programs, so OpenFOAM must be available in the same
shell. Start an OpenFOAM shell with the environment active:

```bash
conda activate cfd
openfoam            # macOS; on Linux, source OpenFOAM's etc/bashrc instead
```

Run all commands below in that shell, from inside a base case folder
(e.g. `annulus_stream-plane/`).

## Workflow

### 1. Generate a case

```bash
pyFoamPrepareCase.py cases/<name> --clone-case=base/annulus-2D
```

This copies the template to `cases/<name>` and then:

1. Reads `default.parameters` and applies any overrides (see step 2).
2. Runs `derivedParameters.py`, which computes SI values (`Ri`, `Ro`, `ecc`,
   `omega`, `nu`).
3. Fills in every `*.template` file: `foo.template` is written as `foo`.
4. Runs `meshCreate.sh`: gmsh builds the mesh, `gmshToFoam` converts it,
   patch types are set, and `checkMesh` runs.
5. Copies `0.org/` to `0/`.
6. Runs `caseSetup.sh`, if the template has one.

The parameter values used are saved in `cases/<name>/PyFoamPrepareCaseParameters`
(readable version: `.rst`). Mesh output is in `meshCreate.sh.log`.

### 2. Override parameters

For a single run:

```bash
pyFoamPrepareCase.py cases/ecc02 --clone-case=base/annulus-2D \
    --values-string="{'eccentricity_mm':0.2, 'speed_rpm':200}"
```

For runs you want to keep a record of, write a parameter file in
`studies/` (OpenFOAM dictionary syntax, no header):

```
// studies/ecc02.parameters
eccentricity_mm 0.2;
speed_rpm       200;
```

```bash
pyFoamPrepareCase.py cases/ecc02 --clone-case=base/annulus-2D \
    --parameter-file=studies/ecc02.parameters
```

`--parameter-file` can be given more than once; later files override earlier
ones. `--automatic-casename` names the case after the parameter files.

### 3. Run

Run the solver named by `application` in the case's `system/controlDict` (see
[Base cases](#base-cases)):

```bash
cd cases/<name>
pimpleFoam > log.pimpleFoam       # transient, e.g. annulus-2D
simpleFoam > log.simpleFoam       # steady, e.g. wedge-2D
```

Results:

- Time folders: full fields, saved every `writeInterval`. A transient solver
  names them by time (`0.01/`, `0.02/`, ...). A steady solver names them by
  iteration (`500/`, `1000/`, ...) and writes the last iteration when it stops.
- `postProcessing/`: time histories, e.g. `bobTorque/0/moment.dat` (torque
  on the bob) and `residuals/0/solverInfo.dat`.
- `log.<solver>`: the solver log.

Open `<name>.foam` in ParaView to view the fields.

### 4. Monitor a run

In a second terminal, with the same environment active:

```bash
python -m cfdtools.monitor cases/<name>              # refresh every 1 s
python -m cfdtools.monitor cases/<name> -r 5         # refresh every 5 s
python -m cfdtools.monitor cases/<name> -c 2         # two plots per row
python -m cfdtools.monitor cases/<name> -n 2000      # only the last 2000 iterations (or seconds)
python -m cfdtools.monitor cases/<name> --save m.png # save one image and exit
python -m cfdtools.monitor cases/rpm_sweep          # a folder of cases: follow the running one
```

The path can be relative to wherever you run the monitor. Given a folder of
cases (a sweep or mesh study), the monitor shows the case whose
`postProcessing` files changed most recently, and switches when the next case
starts. The case shown is in the title.

The monitor plots every `postProcessing/*/*.dat` file. Residuals are shown
on a log scale. Columns that are negligible compared with the largest one in
the same file (e.g. x and y torque in a 2D case) are left out.

For transient cases (`pimpleFoam`), the flow is fully developed when the
torque stops changing; residuals only show that each time step converged.
For steady cases (`simpleFoam`), check that the residuals have fallen and the
torque has stopped changing.

PyFoam's `pyFoamPlotRunner.py` and OpenFOAM's `foamMonitor` are not used.
On macOS, PyFoam's plotting does not work and `foamMonitor` needs X11,
gnuplot and GNU coreutils.

### 5. Re-run or rebuild a case

Running the solver again in the same case does not replace
`postProcessing/`. OpenFOAM writes new files next to the old ones
(`moment_0.dat`, ...). The monitor plots only the newest run.

To delete results but keep the mesh and setup:

```bash
pyFoamClearCase.py cases/<name>        # removes time folders (except 0) and postProcessing/
```

To regenerate a case with the case's own parameters and templates:

```bash
pyFoamPrepareCase.py cases/<name>
pyFoamPrepareCase.py cases/<name> --no-mesh-create      # keep the existing mesh
pyFoamPrepareCase.py cases/<name> --values-string="{'speed_rpm':200}"
```

This does **not** pick up changes made to `base/` after the case was
created. To use the current template, delete the case and generate it again:

```bash
mv cases/<name> ~/.Trash/
pyFoamPrepareCase.py cases/<name> --clone-case=base/annulus-2D
```

To rebuild only the mesh, run `./meshCreate.sh` inside the case.

### 6. Sweeps

Loop over values in the shell:

```bash
for e in 0.0 0.1 0.2 0.3 0.4; do
  pyFoamPrepareCase.py cases/ecc_$e --clone-case=base/annulus-2D \
      --values-string="{'eccentricity_mm':$e}"
  (cd cases/ecc_$e && pimpleFoam > log.pimpleFoam)
done
```

#### Sweep with continuation (`annulus_span-plane/studies/rpm_sweep.sh`)

Each speed starts from the converged solution of the previous speed, the way
the flow develops when the speed is increased slowly. This matters above
Taylor-vortex onset, where more than one steady solution can exist. The
sweep is defined by three files in `studies/`:

| File | Contents |
|---|---|
| `rpm_sweep.speeds` | Speeds in rpm, one per line, run in the order listed |
| `rpm_sweep.parameters` | Settings shared by every case, e.g. `maxIterations` and the mesh |
| `rpm_sweep.sh` | The script that builds and runs the cases |

Run from `annulus_span-plane/` in the OpenFOAM shell:

```bash
bash studies/rpm_sweep.sh
```

For each speed, the script:

1. builds `cases/rpm_sweep/rpm_<speed>` from `base/wedge-2D` with the
   `.parameters` file and that speed;
2. copies the previous case's solution into the new case's `0/` with
   `mapFields` (the first speed starts from rest). The new case keeps its own
   boundary conditions, so the rotating wall uses the new speed;
3. runs `simpleFoam`, and stops the sweep if it fails or does not converge
   within `maxIterations`.

- A case that already exists and converged is skipped and used as the
  starting point for the next speed. Running the script again continues the
  sweep, e.g. after adding speeds to the end of `rpm_sweep.speeds`.
- To redo a case, move it to the Trash and run the script again. Every case
  after it started from the old result, so move those to the Trash as well.
- Logs in each case: `log.prepareCase`, `log.mapFields` (names the source
  case), `log.simpleFoam`.

PyFoam's `pyFoamRunParameterVariation.py` can also build and run a sweep from
a `.variations` file, but every case then starts from rest.

#### Mesh study (`annulus_span-plane/studies/mesh_study.sh`)

Run before a sweep to choose its mesh. The study runs one speed (the top
speed of the sweep, the hardest case to resolve) on three meshes, each with
1.5x the cells of the previous one in every direction.

| File | Contents |
|---|---|
| `mesh_study.meshes` | One mesh per line: name, `nr_gap`, `nr_core`, `nz_ann`, `nz_bot`, run coarse to fine |
| `mesh_study.parameters` | Settings shared by every mesh: `speed_rpm`, `maxIterations` |
| `mesh_study.sh` | The script that builds and runs the cases |

```bash
bash studies/mesh_study.sh
```

It works like `rpm_sweep.sh`, with cases in `cases/mesh_study/<mesh name>`.
Each finer mesh starts from the coarser mesh's solution (`mapFields`
interpolates between the meshes), so all meshes stay on the same vortex
solution.

Compare on each mesh:

- the torque on the bob side and bottom;
- the number of sign changes of the radial velocity up the middle of the gap,
  which must be the same on every mesh (the same vortex pattern);
- the peak radial velocity (vortex strength);
- the peak wall shear stress on the bob (from the `wallShearStress` field).

```python
from pathlib import Path
import numpy as np
from cfdtools import post
from cfdtools.monitor import find_series, read_dat
from cfdtools.convergence import gci

names = ["coarse", "medium", "fine"]
cells, torque = [], []
for name in names:
    case = Path("cases/mesh_study") / name
    p = post.case_parameters(case)
    k, rho = float(p["wedge_to_full"]), float(p["rho_kg_m3"])
    mesh, _ = post.read_case(case)
    s = find_series(case)
    T = read_dat(s["rotorSideTorque/moment"])["total_z"][-1] * k

    r_mid = 0.5 * (float(p["Ri"]) + float(p["Ro"]))
    z0, z1 = mesh.bounds[4], mesh.bounds[5]
    _, ur = post.sample_line(mesh, (r_mid, 0, z0), (r_mid, 0, z1), "U", which="x", n=2000)
    sign_changes = np.count_nonzero(np.diff(np.sign(ur)))

    side, _ = post.read_patch(case, "rotorSide")
    tau_max = np.linalg.norm(side.cell_data["wallShearStress"], axis=1).max() * rho

    print(f"{name}: {mesh.n_cells} cells, T = {T:.5e} N m, sign changes = {sign_changes}, "
          f"max |u_r| = {abs(ur).max():.4e} m/s, max tau = {tau_max:.4g} Pa")
    cells.append(mesh.n_cells)
    torque.append(T)

print(gci(torque, cells))       # dim=2 for the wedge
```

`gci()` (in `cfdtools.convergence`) applies the Grid Convergence Index
procedure of Celik et al. (2008), *J. Fluids Eng.* 130, 078001. It returns
the observed order `p`, the `extrapolated` zero-cell-size value, and
`gci_fine`, the relative uncertainty of the fine-mesh value. Use the coarsest
mesh whose value is within your tolerance of the extrapolated value, and put
its counts in `rpm_sweep.parameters`. If `oscillatory` is `True`, the values
go up and down with refinement, and the GCI is unreliable: add a finer mesh.

### 7. Parallel runs

```bash
pyFoamDecompose.py . 4                        # write system/decomposeParDict, split the mesh into 4
mpirun -np 4 pimpleFoam -parallel > log.pimpleFoam
reconstructPar                                # merge processor*/ back into time folders
```

Parallel runs are only faster for large meshes: at least about 20,000 cells
per process. On Apple silicon, start with one process per performance core.

### 8. Figures

ParaView (open `<name>.foam`) is for looking at results. For figures in a
paper, use the scripts in `post/`, which use `cfdtools.post`:

```bash
python post/couette_figure.py cases/<name>                 # writes cases/<name>/figures/couette.pdf and .png
python post/couette_figure.py cases/<name> --time 0.1 -o fig.pdf
```

`couette_figure.py` (annulus-2D) draws:

- (a) velocity magnitude on the mesh cells, with the annulus unrolled
  (angle against radius) so the gap is visible;
- (b) tangential velocity across the gap. For a concentric case: cell values
  against the analytical Couette solution. For an eccentric case: profiles
  at the narrowest and widest gap.

To open a case in ParaView with a publication view already set up (white
background, parallel projection, no orientation axes, flat colours, the
Cool to Warm colour map, a colour bar, camera along the mesh's thin direction):

```bash
python -m cfdtools.pv cases/<name>                 # ParaView window
python -m cfdtools.pv                              # from inside a case folder
python -m cfdtools.pv cases/<name> --save          # image only, no window: cases/<name>/figures/paraview_U.png
python -m cfdtools.pv cases/<name> --save --field p --time 0.1 --range 0 0.1 --size 3000 2000 -o fig.png
```

It looks for ParaView in `$PARAVIEW_BIN` (a folder containing `paraview` and
`pvbatch`), then on the `PATH`, then in `/Applications/ParaView-*.app` on
macOS. The view is defined in `src/cfdtools/paraview_view.py`, which runs in
ParaView's own Python.

`annulus_stream-plane/post/post_examples.ipynb` shows each `cfdtools.post`
function on a case, including time histories from `postProcessing/`. To run
notebooks, install the kernel once with `pip install -e ".[notebook]"`, then
select the `cfd` environment's Python as the kernel.

`cfdtools.post` functions for writing other figure scripts:

| Function | Does |
|---|---|
| `read_case(case, time="latest")` | Reads the mesh and cell fields at one saved time (via pyvista) |
| `case_parameters(case)` | Returns the parameter values the case was built with |
| `use_style()` | Sets matplotlib defaults for print figures |
| `plot_field(ax, mesh, "U", which="mag")` | Draws a field on the cells of a 2D mesh; `transform=` maps coordinates, e.g. to unroll an annulus |
| `sample_line(mesh, start, end, "U")` | Samples a field along a line, linearly interpolated between cells |
| `cell_centres(mesh)`, `cell_vertices(mesh)` | Cell geometry |

`sample_line` interpolates between cell centres, so values within half a
cell of a wall are not exact (e.g. the velocity at the bob surface reads
slightly below the wall speed). Use cell values when that matters.

## Editing a case template

- Every file in the template is a normal OpenFOAM file and can be edited
  directly.
- Only files that contain parameters are templates (`*.template`). A
  placeholder is written `|-name-|`, where `name` is defined in
  `default.parameters` or `derivedParameters.py`. Python expressions work,
  e.g. `|-2*omega-|`. A line starting with `$$` is a Python assignment.
- To turn a fixed value into a parameter:
  1. Add it to `default.parameters` in lab units.
  2. Convert it in `derivedParameters.py`, if needed.
  3. Replace the value in the file with `|-name-|`.
  4. Rename the file to `<file>.template`.
- In `derivedParameters.py`, `del` any temporary names (e.g. imported
  modules) at the end. Otherwise they are saved as parameters.
- Name the initial-conditions folder `0.org`. `pyFoamPrepareCase.py
  --clone-case` does not copy `0.orig`.
- `meshCreate.sh.template` must be executable (`chmod +x`). The generated
  `meshCreate.sh` gets the same permissions.
- `meshCreate.sh` looks for a `geometry/` folder in the case folder and then
  in each parent folder, and uses the first one it finds. For the current
  base cases that is the top-level `CFD/geometry/`.

## Adding a base case

1. Create the folder:

   ```
   <new base case>/
     base/<template>/
     studies/
     cases/.gitkeep
   ```

2. Add the mesh builder to the top-level `geometry/` folder. It must write
   MSH 2.2 ASCII (the only format `gmshToFoam` reads). Physical group names
   become the OpenFOAM patch names. Put dimensions shared with other
   builders in `geometry/common.py`.
3. Copy the OpenFOAM tutorial closest to the new flow into
   `base/<template>/`:
   - `foamSearch $FOAM_TUTORIALS system/controlDict application` lists the
     solver each tutorial uses.
   - Rename `0` or `0.orig` to `0.org`.
   - Delete the tutorial's `Allrun`, `Allclean` and `blockMeshDict`.
   - Rename the boundary conditions in `0.org/*` to match your patch names.
4. Copy `meshCreate.sh.template` from `annulus_stream-plane/base/annulus-2D/`.
   Change the mesh builder command and the patch types (`empty`, `wedge`,
   `wall`, ...).
5. Get the case running with fixed values first, without parameters or
   templates.
6. Replace values that should vary with parameters (see
   [Editing a case template](#editing-a-case-template)).
7. Add a `solverInfo` function object (for the monitor) and function objects
   for the quantities you need (`forces`, `probes`, `fieldAverage`, ...).
8. Test it:

   ```bash
   pyFoamPrepareCase.py cases/test --clone-case=base/<template> --values-string="{'endTime':0.01}"
   cd cases/test && pimpleFoam > log.pimpleFoam
   ```

## Base cases

| Base case | Template | Mesh | Solver | Patches | Status |
|---|---|---|---|---|---|
| `annulus_stream-plane` | `annulus-2D` | `geometry/annulus.py`: 2D r-theta annulus, optional eccentricity | pimpleFoam, laminar, Newtonian | `innerWall` (rotating), `outerWall`, `frontAndBack` (empty) | Working |
| `annulus_span-plane` | `wedge-2D` | `geometry/wedge.py`: axisymmetric r-z wedge (5°), optional bottom gap | simpleFoam (SIMPLEC), laminar, Newtonian | `rotorSide`, `rotorBottom` (rotating), `cupWall`, `cupBottom`, `top` (slip), `front`, `back` (wedge) | Working. No eccentricity (axisymmetric). Torques are for the wedge; multiply by `wedge_to_full` |

Mesh files (`*.msh`) are not tracked by git. Each case builds its own mesh.
To build one directly, from the repository root, run e.g.
`python geometry/annulus.py --help`.

## License

MIT. See [LICENSE](LICENSE).
