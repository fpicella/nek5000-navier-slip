#!/bin/bash
# lane.sh NP MODE KN1 [KN2 ...] -- run the listed Kn cases one after another
# with NP ranks (MODE = re20 | couette).  Launch several lanes in parallel:
#   setsid nohup bash scripts/lane.sh 4 re20 0.0 0.5 5.0 > logs/lane1.log 2>&1 &
# Optional env: TAG (run-dir prefix, default MODE), MESHOPTS (gen_cyl_re2.py
# options, re20 mode only), e.g. TAG=fine MESHOPTS="--ntheta 56 --nr 34 --dr1 0.05",
# PAREDIT (sed script applied to cyl.par after setup, e.g. to set Re).
HERE=$(cd "$(dirname "$0")/.." && pwd)
NP=$1; MODE=$2; shift 2
TAG=${TAG:-$MODE}
unset DISPLAY
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 HWLOC_COMPONENTS=-gl
for KN in "$@"; do
  RUN=$HERE/runs/${TAG}_kn$KN
  if [ "$MODE" = couette ]; then bash "$HERE/scripts/setup_couette.sh" "$RUN" "$KN"
  else bash "$HERE/scripts/setup_run.sh" "$RUN" "$KN" $MESHOPTS; fi
  [ -n "$PAREDIT" ] && sed -i -e "$PAREDIT" "$RUN/cyl.par"
  cd "$RUN" && rm -f cyl_drag.dat
  echo "[$(date +%H:%M:%S)] start $TAG Kn=$KN np=$NP"
  mpirun --bind-to none -np "$NP" ./nek5000 > run.log 2>&1 < /dev/null
  echo "[$(date +%H:%M:%S)] done  $TAG Kn=$KN  $(grep cylmon run.log | tail -1)"
done
