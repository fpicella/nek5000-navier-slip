#!/bin/bash
# setup_couette.sh RUNDIR KN -- circular-Couette analytic check (annulus a=0.5,
# R=5, outer wall rotating at omega=1, viscosity 1, slip Kn on the cylinder).
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
RUN=$1; KN=$2
bash "$HERE/scripts/setup_run.sh" "$RUN" "$KN" --rinf 5 --nr 12 --dr1 0.08 > /dev/null
sed -i -e "s/^endTime *=.*/endTime = 60.0/" -e "s/^dt *=.*/dt = 1.0e-2/" -e "s/^viscosity *=.*/viscosity = 1.0/" \
       -e "s/^userParam02 *=.*/userParam02 = 1.0e-7/" -e "s/^userParam03 *=.*/userParam03 = 100/" \
       -e "/^userParam03/a userParam04 = 1\nuserParam05 = 1.0" "$RUN/cyl.par"
echo "ready: $RUN (Couette check, Kn=$KN)"
