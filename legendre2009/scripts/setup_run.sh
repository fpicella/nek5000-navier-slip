#!/bin/bash
# setup_run.sh RUNDIR KN [gen_cyl_re2.py options...]
#   Creates RUNDIR with the mesh (.re2/.ma2/.co2), par (userParam01=KN),
#   and a copy of the compiled executable from $PROJECT/build_<host>.
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
NEK=$HERE/Nek5000
RUN=$1; KN=$2; shift 2
mkdir -p "$RUN" && cd "$RUN"
cp "$HERE/case/cyl.par" .
sed -i -e "s/^userParam01 *=.*/userParam01 = $KN/" cyl.par
python3 "$HERE/mesh/gen_cyl_re2.py" -o cyl "$@" > mesh.log
head -2 mesh.log
printf "cyl\n0.2\n" | "$NEK/bin/gencon" > gencon.log 2>&1 || true
printf "cyl\n0.2\n" | "$NEK/bin/genmap" > genmap.log 2>&1
ls cyl.re2 cyl.ma2 > /dev/null
cp "$HERE/build_$(hostname -s)/nek5000" .   # executable built on this host
printf "cyl\n%s/\n" "$(pwd)" > SESSION.NAME
echo "ready: $(pwd)  Kn=$KN"
