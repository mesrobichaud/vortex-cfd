#!/bin/bash
# Rotation-rate sweep with continuation: each speed starts from the converged
# solution of the previous speed. The first speed starts from rest.
#
# Run in the OpenFOAM shell with the cfd environment active:
#   bash studies/rpm_sweep.sh
#
# Speeds:          studies/rpm_sweep.speeds (run in the order listed)
# Shared settings: studies/rpm_sweep.parameters
# Cases:           cases/rpm_sweep/rpm_<speed>
#
# A case that already exists and converged is skipped and used as the starting
# point for the next speed, so running the script again continues the sweep.
# The script stops at the first case that fails or does not converge.

cd "$(dirname "$0")/.." || exit 1       # annulus_span-plane/

base=base/wedge-2D
parameters=studies/rpm_sweep.parameters
speeds=studies/rpm_sweep.speeds
out=cases/rpm_sweep

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

for rpm in $(sed 's|//.*||' "$speeds"); do
    case=$out/rpm_$rpm

    if [ -d "$case" ]; then
        converged "$case" || fail "$case exists but did not converge. Move it to the Trash and run again."
        echo "$case: already converged, skipping"
        prev=$case
        continue
    fi

    echo "$case: building"
    log=$(pyFoamPrepareCase.py "$case" --clone-case="$base" \
              --parameter-file="$parameters" \
              --values-string="{'speed_rpm':$rpm}" 2>&1)
    [ -f "$case/constant/polyMesh/boundary" ] || { echo "$log"; fail "$case: case or mesh not created"; }
    echo "$log" > "$case/log.prepareCase"

    if [ -n "$prev" ]; then
        echo "$case: mapping fields from $prev"
        mapFields "$PWD/$prev" -case "$case" -consistent -sourceTime latestTime \
            > "$case/log.mapFields" 2>&1 || fail "$case: mapFields failed, see $case/log.mapFields"
    fi

    echo "$case: running simpleFoam"
    simpleFoam -case "$case" > "$case/log.simpleFoam" 2>&1 || fail "$case: simpleFoam failed, see $case/log.simpleFoam"
    converged "$case" || fail "$case: did not converge within maxIterations, see $case/log.simpleFoam"
    echo "$case: $(grep "SIMPLE solution converged" "$case/log.simpleFoam")"

    prev=$case
done

echo "Sweep complete: $out"
