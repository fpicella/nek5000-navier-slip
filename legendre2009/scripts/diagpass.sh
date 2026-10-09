#!/bin/bash
# diagpass.sh NP RUNDIR... -- restart each finished run from its last field
# for one step with the current executable and userParam08 = 1, in
# RUNDIR/diag, to write the wall distribution (wall_<rank>.dat).
# NEK_EXE = executable (default ../../cylinder2d/nek5000).  Needs a .usr with the
# wall dump (userParam08), which the cylinder2d/cyl.usr of this repository does not have.
P=$(cd "$(dirname "$0")/.." && pwd); EXE=${NEK_EXE:-$(dirname "$P")/cylinder2d/nek5000}; NP=$1; shift; ORIG=$(pwd)
for D in "$@"; do
  cd "$ORIG"; D=$(cd "$D" && pwd); last=$(ls "$D"/cyl0.f0* 2>/dev/null | tail -1)
  [ -n "$last" ] || { echo "no field file in $D"; continue; }
  rm -rf "$D/diag"; mkdir -p "$D/diag"; cd "$D/diag"
  cp "$last" restart.f00001; cp ../cyl.re2 ../cyl.ma2 ../cyl.par .; cp "$EXE" .
  sed -i -e 's/^stopAt *=.*/stopAt = numSteps\nnumSteps = 1\nstartFrom = restart.f00001/' -e '/^endTime/d' \
         -e 's/^userParam02 *=.*/userParam02 = 0.0/' -e 's/^userParam03 *=.*/userParam03 = 1/' \
         -e '/^userParam0[78]/d' cyl.par
  sed -i -e 's/^\(userParam03 *=.*\)$/\1\nuserParam08 = 1/' cyl.par
  printf "cyl\n%s/\n" "$PWD" > SESSION.NAME
  OMP_NUM_THREADS=1 HWLOC_COMPONENTS=-gl mpirun --bind-to none -np "$NP" ./nek5000 > run.log 2>&1 < /dev/null
  grep -q "run successful" run.log && echo "ok   $(basename "$D") ($(cat wall_*.dat | wc -l) wall nodes)" || echo "FAIL $(basename "$D")"
done
