# Comparison with Legendre, Lauga & Magnaudet (2009)

D. Legendre, E. Lauga and J. Magnaudet, *Influence of slip on the dynamics of two-dimensional
wakes*, J. Fluid Mech. 633, 437–447 (2009), arXiv:0905.0648.

122 runs of the 2D cylinder case (`../cylinder2d`, 1040 elements, N = 7) from Re = 5.5 to 800,
plus two runs on a refined mesh (1904 elements) at Re = 800. The main results are summarised
in the [top-level README](../README.md#validation-in-brief).

## Contents

| Path | Content |
|---|---|
| `data/legendre_campaign.txt` | one line per run (Cd, CL, St, lift growth rate, wall vorticity slope and separation angle), then the derived separation and shedding boundaries |
| `data/legendre_fig*.dat`, `data/legendre_re100_*.dat` | curves and markers digitised from the paper (`scripts/digitize_fig3.py`, `digitize_fig4.py` read them from the arXiv PDF with PyMuPDF) |
| `data/re20_sweep*.dat`, `data/re100_sweep.dat` | Re = 20 and Re = 100 sweeps in Kn |
| `data/legendre_fig134_summary.txt` | figures 1, 3, 4 against this work, in numbers |
| `data/timing_re20.txt` | cost of the slip wall: no slip vs Navier slip, same mesh and time step |
| `figures/` | the comparison figures |
| `scripts/` | run set-up, job queue and plotting scripts (below) |
| `rerun_2026-10-08/` | independent re-run of the Re = 20 sweep and of the Couette check, on a fresh Nek5000 (commit in `nek5000_commit.txt`, patch in `subs1_robin_patch.diff`), with its analysis |

## What is not here

- The raw time series (`cyl_drag.dat` of every run, about 120 MB) and the flow fields. The
  plotting scripts read them from `runs/<family><Re>_kn<Kn>/`, so they only run after the
  campaign has been run again.
- The linear stability (nekStab) results used by `scripts/plot_legendre_lsa.py` for figures 1
  (`legendre_fig1_lsa.png`), 3 and 4.
- The campaign `.usr` file wrote extra wall diagnostics (vorticity slope at the rear stagnation
  point, separation angle, wall and vorticity dumps). `../cylinder2d/cyl.usr` is a cleaned-up
  version without them: with it you get Cd, CL, St and the shedding threshold, but not the
  separation boundary (figure 1), figure 4 or the vorticity insets. `scripts/diagpass.sh` and the
  Couette mode of `scripts/setup_couette.sh` / `lane.sh couette` need the same extended file.

## Running the campaign again

```bash
export NEK5000=/path/to/Nek5000                  # patched by cylinder2d/build.sh
../cylinder2d/build.sh $NEK5000                  # executable ../cylinder2d/nek5000 (or set NEK_EXE)
python3 scripts/campaign.py                       # list the runs and their settings
python3 scripts/campaign.py setup                 # runs/L<Re>_kn<Kn>/ and jobs_<host>.txt
xargs -P 8 -I{} bash -c '{}' < jobs_host1.txt     # 4 MPI ranks per job
```

Edit `HOSTS` in `campaign.py` to split the jobs over your machines (name: number of 4-rank job
slots). `scripts/setup_run.sh RUNDIR KN [mesh options]` sets up a single run; `scripts/lane.sh`
runs a list of Kn one after the other (the Re = 20 sweep: `runs/re20_kn*`).

## Plotting

All plotting scripts write to `data/` and `figures/` and take the runs folder as first argument
(default `runs/`). Python 3 with NumPy, SciPy and Matplotlib.

| Script | Reads | Writes |
|---|---|---|
| `plot_fig2a.py` | `runs/re20_kn*`, `runs/fine_kn*` | `re20_sweep.dat`, `fig2a_re20.png` |
| `plot_re100.py` | `runs/re100_kn*` | `re100_sweep.dat`, `re100_vs_legendre.png` |
| `plot_legendre.py` | `runs/{re20,re100,L,Lfine}*_kn*` | `legendre_campaign.txt`, `legendre_fig1.png`, `legendre_fig2.png`, `legendre_fig1_insets.png` |
| `plot_legendre_lsa.py NEKSTAB_RESULTS` | `legendre_campaign.txt`, nekStab results | `legendre_fig1_lsa.png`, `legendre_fig3.png`, `legendre_fig4.png`, `legendre_fig134_summary.txt` |
| `timing_parse.py` | `runs/timing` (from `timing_test.sh`) | table on screen |
| `rerun_2026-10-08/rerun_analysis.py` | `rerun_2026-10-08/runs`, `data/` | files in `rerun_2026-10-08/` |
