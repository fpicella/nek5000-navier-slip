# Verification against exact solutions

Annulus a = 0.5 < r < R = 1.5, with the Navier slip condition on the inner cylinder and
Dirichlet conditions on the outer one. `verif.usr` holds the exact solutions:

| case | flow | what it tests |
|---|---|---|
| 1 | Stokes, outer cylinder translating along x | slip velocity that varies along the wall; force on the cylinder |
| 2 | decaying azimuthal mode (Navier–Stokes) | time accuracy: the decay rate σ = ν k² |
| 3 | circular Couette (steady Navier–Stokes) | the curvature term of the slip condition; torque |
| 4 | axial Couette, 3D only | the Robin term on the spanwise velocity |

## Results in this folder

- `data/<run>/verif_summary.dat`: final errors of every run (`data/c3_*` also has `profile.dat`,
  the velocity at all grid points). Run names are explained at the top of `make_runs.py`.
- `verif_results.txt`, `verif_convergence.png`: p-, h- and dt-refinement (cases 1, 2, 3).
- `couette_validation.png`: the circular Couette flow, Kn = 1e-3 to 1e3.

The 3D runs (`v3*`) can be set up with `make_runs.py v3`, but their results are not in `data/`.

## Remaking the figures from the shipped data

```bash
python3 plot_verif.py       # reads data/, writes verif_results.txt and verif_convergence.png
python3 plot_couette.py     # reads data/, writes couette_validation.png
```

Needs Python 3 with NumPy and Matplotlib. `data/` also holds six `p1tight_*` runs (case 1 with
solver tolerances 1e-14); `plot_verif.py` adds them to the table and to panels (a)–(b).

## Running the suite again

```bash
export NEK5000=/path/to/Nek5000             # a clone dedicated to slip cases: core/subs1.f is patched
for N in 3 5 6 7 8 9 10 11 12 13; do ./build_verif.sh $N; done     # build_N<N>/nek5000
./build_verif.sh 4 1024                      # N = 4 also runs the h-refinement, up to 1024 elements
./build_verif.sh 7 64 3                      # only for the 3D runs (build_N7_3d)
python3 make_runs.py p1 p3 c3 h1 t2          # runs/<name>/ and jobs.txt; ALWAYS give name prefixes
xargs -P 16 -I{} bash -c '{}' < jobs.txt     # one MPI job per line
python3 plot_verif.py runs; python3 plot_couette.py runs
```

`build_verif.sh N [LELG] [DIM]` starts from `../cylinder2d/SIZE` and `../src/navier_slip.f`.
`make_runs.py` deletes the old results of every run it sets up again, so always pass the
prefixes of the runs you want.
