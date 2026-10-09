#!/bin/bash
# build_verif.sh N [LELG] [DIM] -- compile verif.usr at polynomial order N into
# verif/build_N<N> (DIM = 3: verif/build_N<N>_3d, ldim = 3).  Needs a host with MPI headers.
# (Tested with Ubuntu 22.04, gfortran, OpenMPI.)
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
NEK=$HERE/Nek5000
N=$1; LELG=${2:-64}; DIM=${3:-2}
LX1=$((N+1)); LXD=$(( (3*LX1+1)/2 ))
python3 "$HERE/patch/apply_robin_patch.py" "$NEK/core/subs1.f" > /dev/null
B=$HERE/verif/build_N$N; [ "$DIM" = 3 ] && B=${B}_3d
mkdir -p "$B" && cd "$B"
cp "$HERE/verif/verif.usr" "$HERE/case/navier_slip.f" .
sed -e "s/parameter (lx1=[0-9]*)/parameter (lx1=$LX1)/" \
    -e "s/parameter (ldim=[0-9]*)/parameter (ldim=$DIM)/" \
    -e "s/parameter (lxd=[0-9]*)/parameter (lxd=$LXD)/" \
    -e "s/parameter (lelg=[0-9]*)/parameter (lelg=$LELG)/" \
    -e "s/parameter (lpmin=[0-9]*)/parameter (lpmin=1)/" \
    "$HERE/case/SIZE" > SIZE
export NEK_SOURCE_ROOT=$NEK FC=mpif90 CC=mpicc
"$NEK/bin/makenek" verif > build.log 2>&1 || {
  grep -n -i -B2 -A4 "error" build.log | head -40; exit 1; }
echo "N=$N lx1=$LX1 lxd=$LXD lelg=$LELG ldim=$DIM built in $B"
