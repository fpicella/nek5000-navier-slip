#!/bin/bash
# run.sh [NPROCS] -- run the case on NPROCS MPI ranks (default 4, at least 2 with this SIZE).
# Nek's log goes to logfile; forces and wall diagnostics to cyl_drag.dat (see cyl.usr).
cd "$(dirname "$0")" || exit 1
NP=${1:-4}
printf "cyl\n%s/\n" "$(pwd)" > SESSION.NAME
mpirun -np "$NP" ./nek5000 > logfile 2>&1
grep -a -E "steady state|run successful" logfile
tail -1 cyl_drag.dat
