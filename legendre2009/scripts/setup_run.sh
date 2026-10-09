#!/bin/bash
# setup_run.sh RUNDIR KN [gen_cyl_re2.py options...]
#   Creates RUNDIR with the mesh (.re2/.ma2/.co2), par (userParam01=KN),
#   and a copy of the compiled executable.
#   Env: NEK5000 = path to the (patched) Nek5000 tree, for gencon/genmap (required);
#        NEK_EXE = executable to copy (default ../../cylinder2d/nek5000, built by cylinder2d/build.sh).
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
REPO=$(dirname "$HERE")
NEK=${NEK5000:?set NEK5000=/path/to/Nek5000}
EXE=${NEK_EXE:-$REPO/cylinder2d/nek5000}
[ -x "$EXE" ] || { echo "no executable $EXE: run cylinder2d/build.sh first (or set NEK_EXE)"; exit 1; }
RUN=$1; KN=$2; shift 2
mkdir -p "$RUN" && cd "$RUN"
cp "$REPO/cylinder2d/cyl.par" .
sed -i -e "s/^userParam01 *=.*/userParam01 = $KN/" cyl.par
python3 "$REPO/cylinder2d/gen_cyl_re2.py" -o cyl "$@" > mesh.log
head -2 mesh.log
printf "cyl\n0.2\n" | "$NEK/bin/gencon" > gencon.log 2>&1 || true
printf "cyl\n0.2\n" | "$NEK/bin/genmap" > genmap.log 2>&1
ls cyl.re2 cyl.ma2 > /dev/null
cp "$EXE" .
printf "cyl\n%s/\n" "$(pwd)" > SESSION.NAME
echo "ready: $(pwd)  Kn=$KN"
