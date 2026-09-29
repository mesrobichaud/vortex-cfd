# CFD

OpenFOAM cases for the CPB rig, meshed with gmsh and set up with PyFoam.

Every run is built from a **base case** (a normal OpenFOAM case with a few
templated values) plus a **parameter set**. Each step -- mesh, fields,
physics, numerics -- is its own file that can be edited by hand or driven
from parameters. There is no single script that builds everything.

## Layout

```
CFD/
  cfdtools/                 shared Python helpers (import cfdtools.*)
    mesh.py                 gmsh sizing helpers, element table, colour modes
    units.py                lab units -> SI (rpm, mm, mPa.s, nu = mu/rho)
  pyproject.toml            makes cfdtools pip-installable

  couette_ideal/            a case-group: one rig geometry, many runs
    geometry/               gmsh mesh builders (annulus.py, wedge.py) + rig dimensions
    base/                   base cases -- where the physics and numerics live
      annulus2D/            2D r-theta annulus, laminar pimpleFoam
    studies/                parameter files for named runs and sweeps
    cases/                  generated runs (gitignored, disposable)

  true-geometry_ideal/      another case-group, same structure
```

**What goes where**

- `cfdtools/`: code that doesn't depend on any particular geometry. Write
  helpers inside the case-group first, and move them here once a second
  group needs them.
- `<group>/geometry/`: the mesh builders and rig dimensions for that group.
- `<group>/base/<name>/`: the OpenFOAM set-up (BCs, schemes, solver
  settings, function objects). Change the physics here.
- `<group>/studies/`: which values to run. Change parameters here.
- `<group>/cases/`: output only. Anything worth keeping goes back into
  `base/` or `studies/`, so any case can be regenerated.

## One-time setup

```bash
conda activate cfd                  # has gmsh, PyFoam, numpy
pip install -e .                    # from CFD/, installs cfdtools in editable mode
```

OpenFOAM (v2606, macOS app) must be on the path whenever PyFoam runs,
because PyFoam calls the OpenFOAM utilities. Start the OpenFOAM shell from
the activated env:

```bash
conda activate cfd
openfoam                            # opens a shell with OpenFOAM + the cfd env
```

All commands below are run from inside that shell, in a case-group
directory (e.g. `couette_ideal/`).

## Workflow

### 1. Build a case

```bash
pyFoamPrepareCase.py cases/<name> --clone-case=base/annulus2D
```

This copies the base case to `cases/<name>` and then prepares it:

1. reads `default.parameters`, then applies any overrides (see below)
2. runs `derivedParameters.py`, which converts to SI (`Ri`, `Ro`, `ecc`, `omega`, `nu`)
3. fills every `*.template` file, writing `foo.template` out as `foo`
4. runs `meshCreate.sh`, which does gmsh, then `gmshToFoam`, sets the patch
   types, and runs `checkMesh`
5. copies `0.org/` to `0/`
6. runs `caseSetup.sh` if one exists (setFields, decomposePar, ...)

The values actually used are recorded in `cases/<name>/PyFoamPrepareCaseParameters`
(`.rst` has a readable version). Script output goes to `meshCreate.sh.log`, etc.

### 2. Override parameters

Inline, for one-offs:

```bash
pyFoamPrepareCase.py cases/ecc02 --clone-case=base/annulus2D \
    --values-string="{'eccentricity_mm':0.2, 'speed_rpm':200}"
```

From a file, for runs you want to keep track of:

```bash
# studies/ecc02.parameters  (OpenFOAM dictionary syntax, no header)
eccentricity_mm 0.2;
speed_rpm       200;
```
```bash
pyFoamPrepareCase.py cases/ecc02 --clone-case=base/annulus2D \
    --parameter-file=studies/ecc02.parameters
```

You can pass more than one `--parameter-file`; later files win.
`--automatic-casename` names the case after the parameter files used.

### 3. Run

```bash
cd cases/<name>
pimpleFoam > log.pimpleFoam                 # plain
pyFoamPlotRunner.py --progress pimpleFoam   # with live residual plots
```

Results: time directories, and `postProcessing/bobTorque/0/moment.dat` for
the torque on the bob. For ParaView, open `<name>.foam`.

### 4. Re-run one step

Each step can be run again on its own:

```bash
pyFoamPrepareCase.py cases/<name> --no-mesh-create     # new parameters, keep the mesh
pyFoamPrepareCase.py cases/<name> --stop-after-templates
```

You can also edit a generated dict directly in `cases/<name>` for a quick
test. Anything worth keeping belongs in `base/` or `studies/`.

### 5. Sweeps

A shell loop over parameter values:

```bash
for e in 0.0 0.1 0.2 0.3 0.4; do
  pyFoamPrepareCase.py cases/ecc_$e --clone-case=base/annulus2D \
      --values-string="{'eccentricity_mm':$e}"
  (cd cases/ecc_$e && pimpleFoam > log.pimpleFoam)
done
```

PyFoam also has a built-in sweep runner, `pyFoamRunParameterVariation.py
base/annulus2D studies/<sweep>.variations`. The `.variations` file is an
OpenFOAM dict with a `values` sub-dict listing each parameter's values; it
must include a `solver` entry, e.g. `solver (pimpleFoam);`. See `--help`
for the options.

## Editing a base case

- Base cases are ordinary OpenFOAM cases, so you can edit any dict directly.
- Only files containing a varying value are templates. A placeholder is
  `|-name-|`, where `name` is anything from `default.parameters` or
  `derivedParameters.py`. Python expressions also work, e.g.
  `|-2*omega-|`. A line starting with `$$` is a Python assignment.
- To make a fixed value adjustable: add it to `default.parameters` (in lab
  units), convert it in `derivedParameters.py` if needed, put `|-name-|` in the
  dict, and rename the dict to `foo.template`.
- Temporary names in `derivedParameters.py` (e.g. imported modules) must be
  `del`'d at the end, or they end up in the parameter record.
- Keep `0.org` as the folder name. PyFoam's cloning copies `*.org` folders but
  silently skips `0.orig`.
- `meshCreate.sh.template` must stay executable (`chmod +x`); the generated
  script inherits its mode.

## Adding a case-group

```
new_group/
  geometry/       mesh builder(s): write MSH 2.2 ASCII, name patches with physical groups
  base/<name>/    default.parameters, derivedParameters.py, meshCreate.sh.template,
                  0.org/, constant/, system/
  studies/
  cases/.gitkeep
```

Copy `couette_ideal/base/annulus2D` as a starting point. `meshCreate.sh`
walks up from the case directory to the nearest `geometry/`, so the builders
are found whether the case is in `cases/`, `base/` or a sweep directory.

## Existing base cases

| Group | Base | Mesh | Solver | Patches |
|---|---|---|---|---|
| couette_ideal | annulus2D | `geometry/annulus.py`, 2D r-theta, optional eccentricity | pimpleFoam, laminar, Newtonian | `innerWall` (rotating), `outerWall`, `frontAndBack` (empty) |
