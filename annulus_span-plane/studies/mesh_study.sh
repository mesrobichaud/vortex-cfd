#!/bin/bash
# Mesh-independence study: one case per mesh, run coarse to fine. Each finer
# mesh starts from the converged solution of the mesh before it (mapFields
# interpolates between the meshes), so every mesh stays on the same vortex
# solution. The first mesh starts from rest.
#
# Run in the OpenFOAM shell with the cfd environment active:
#   bash studies/mesh_study.sh
#
# Meshes:          studies/mesh_study.meshes (run in the order listed)
# Shared settings: studies/mesh_study.parameters (speed, maxIterations)
# Cases:           cases/mesh_study/<mesh name>
#
# A case that already exists and converged is skipped and used as the starting
# point for the next mesh, so running the script again continues the study.
# The script stops at the first case that fails or does not converge.

cd "$(dirname "$0")/.." || exit 1       # annulus_span-plane/

base=base/wedge-2D
parameters=studies/mesh_study.parameters
meshes=studies/mesh_study.meshes
out=cases/mesh_study

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

converged() {
    grep -q "SIMPLE solution converged" "$1/log.simpleFoam" 2>/dev/null
}

command -v simpleFoam > /dev/null || fail "simpleFoam not found: run inside the OpenFOAM shell"
command -v pyFoamPrepareCase.py > /dev/null || fail "pyFoamPrepareCase.py not found: activate the cfd environment"

mkdir -p "$out"
prev=""

# The list is read on file descriptor 3 so that commands in the loop can't
# read from it.
while read -r name nr_gap nr_core nz_ann nz_bot <&3; do
    case=$out/$name

    if [ -d "$case" ]; then
        converged "$case" || fail "$case exists but did not converge. Move it to the Trash and run again."
        echo "$case: already converged, skipping"
        prev=$case
        continue
    fi

    echo "$case: building (nr_gap $nr_gap, nr_core $nr_core, nz_ann $nz_ann, nz_bot $nz_bot)"
    log=$(pyFoamPrepareCase.py "$case" --clone-case="$base" \
              --parameter-file="$parameters" \
              --values-string="{'nr_gap':$nr_gap, 'nr_core':$nr_core, 'nz_ann':$nz_ann, 'nz_bot':$nz_bot}" 2>&1)
    [ -f "$case/constant/polyMesh/boundary" ] || { echo "$log"; fail "$case: case or mesh not created"; }
    echo "$log" > "$case/log.prepareCase"

    if [ -n "$prev" ]; then
        echo "$case: mapping fields from $prev"
        # interpolate: linear interpolation from the coarser mesh. This OpenFOAM
        # version accepts mapNearest, interpolate or cellPointInterpolate here.
        mapFields "$PWD/$prev" -case "$case" -consistent -sourceTime latestTime \
            -mapMethod interpolate \
            > "$case/log.mapFields" 2>&1 || fail "$case: mapFields failed, see $case/log.mapFields"
    fi

    echo "$case: running simpleFoam"
    simpleFoam -case "$case" > "$case/log.simpleFoam" 2>&1 || fail "$case: simpleFoam failed, see $case/log.simpleFoam"
    converged "$case" || fail "$case: did not converge within maxIterations, see $case/log.simpleFoam"
    echo "$case: $(grep "SIMPLE solution converged" "$case/log.simpleFoam")"

    prev=$case
done 3< <(sed 's|//.*||' "$meshes" | grep -v '^[[:space:]]*$')

echo "Mesh study complete: $out"
