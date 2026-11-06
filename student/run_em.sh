#!/usr/bin/env bash
# Energy-minimise every prepared system (folders made by "fes_prep.py prep").
#
#   bash run_em.sh [systems_dir]
#
# Environment variables (optional):
#   GMX=gmx            GROMACS executable (gmx_mpi on most clusters)
#   MDRUN_FLAGS=""     extra mdrun flags, e.g. "-ntmpi 1 -ntomp 8"
#   MAXWARN=0          grompp -maxwarn
#   FORCE=1            redo systems that already have em.gro
set -uo pipefail
ROOT=${1:-systems}
GMX=${GMX:-gmx}
MDRUN_FLAGS=${MDRUN_FLAGS:-}
MAXWARN=${MAXWARN:-0}
MDP=$(cd "$(dirname "$0")" && pwd)/minim.mdp

for d in "$ROOT"/*/; do
    d=${d%/}
    [ -f "$d/system.top" ] || { echo "skip $d (not prepared, see prep_report.json)"; continue; }
    if [ -f "$d/em.gro" ] && [ -z "${FORCE:-}" ]; then echo "skip $d (already minimised)"; continue; fi
    (
        cd "$d" &&
        $GMX grompp -f "$MDP" -c system.gro -p system.top -o em.tpr -maxwarn "$MAXWARN" > grompp.log 2>&1 &&
        $GMX mdrun -deffnm em $MDRUN_FLAGS > mdrun.log 2>&1
    ) && echo "ok   $d" || echo "FAIL $d (see grompp.log / mdrun.log)"
done
