#!/bin/bash
# build.sh /path/to/Nek5000 -- patch Nek5000 (one core routine, see README.md) and compile this case.
# The patch edits core/subs1.f of the given Nek5000 tree; it is idempotent (does nothing if the tree
# is already patched).  Use a Nek5000 clone dedicated to slip cases.
set -e
NEK=${1:?usage: ./build.sh /path/to/Nek5000}
NEK=$(cd "$NEK" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
python3 apply_robin_patch.py "$NEK/core/subs1.f"
export NEK_SOURCE_ROOT=$NEK FC=${FC:-mpif90} CC=${CC:-mpicc}
"$NEK/bin/makenek" cyl > build.log 2>&1 || { tail -30 build.log; exit 1; }
echo "built $HERE/nek5000 (log in build.log); run it with ./run.sh NPROCS"
