#!/bin/bash
# timing_test.sh -- cost of the Navier slip wall, Re = 20 cylinder run to
# steady state (userParam02 = 1e-6), same mesh / N / dt / tolerances:
#   A  no-slip 'W', stressFormulation = no   (standard Nek set-up)
#   B  no-slip 'W', stressFormulation = yes  (the solver the slip wall needs)
#   C  Kn = 1e-3 ('shl' + implicit Robin), stressFormulation = yes
#   D  as C, wall diagnostics every step (userParam03 = 1) instead of 200
# Two rounds (r1, r2); the four cases of a round run concurrently, NP ranks
# each.  Results: runs/timing/<round>_<case>/run.log (Nek runstat block),
# summary with scripts/timing_parse.py.
P=$(cd "$(dirname "$0")/.." && pwd)
T=$P/runs/timing; NP=${NP:-4}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 HWLOC_COMPONENTS=-gl
unset DISPLAY
mkdir -p "$T"; LOG=$T/timing.log
for r in r1 r2; do for c in A B C D; do
  case $c in
    A) kn=0.0;   sf=no;  np3=200 ;;
    B) kn=0.0;   sf=yes; np3=200 ;;
    C) kn=0.001; sf=yes; np3=200 ;;
    D) kn=0.001; sf=yes; np3=1   ;;
  esac
  d=$T/${r}_$c
  bash "$P/scripts/setup_run.sh" "$d" $kn > /dev/null
  sed -i -e "s/^stressFormulation *=.*/stressFormulation = $sf/" \
         -e "s/^userParam03 *=.*/userParam03 = $np3/" \
         -e "s/^endTime *=.*/endTime = 300.0/" "$d/cyl.par"
done; done
echo "$(date +%T) $(hostname -s) setup done" >> "$LOG"
for r in r1 r2; do
  for c in A B C D; do
    ( cd "$T/${r}_$c" && mpirun --bind-to none -np $NP ./nek5000 > run.log 2>&1 < /dev/null ) &
  done
  wait
  echo "$(date +%T) round $r done" >> "$LOG"
done
echo ALLDONE >> "$LOG"
