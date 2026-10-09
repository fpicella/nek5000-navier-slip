# Navier slip on curved walls in Nek5000

An implicit Navier (partial-slip) boundary condition for the spectral-element code
[Nek5000](https://github.com/Nek5000/Nek5000), exact on curved walls, with:

- a ready-to-run test case: 2D flow past a circular cylinder;
- a verification suite against exact solutions;
- a full comparison with Legendre, Lauga & Magnaudet (2009), *Influence of slip on the dynamics of
  two-dimensional wakes*, J. Fluid Mech. 633, 437–447.

On the wall, with n the unit normal pointing into the fluid and S the rate-of-strain tensor:

```
u · n = 0 ,        u_t = 2 λ (S · n)_t
```

λ is the slip length (λ = 0: no slip; λ → ∞: shear free). The curvature term is included, and
the slip friction is treated implicitly, so any slip length runs with the usual time step.

## How to cite

If you use this code, please cite:

> F. Picella, J.-Ch. Robinet and S. Cherubini,
> *Laminar–turbulent transition in channel flow with superhydrophobic surfaces modelled as a partial
> slip wall*, J. Fluid Mech. **881**, 462–497 (2019).
> [doi:10.1017/jfm.2019.740](https://doi.org/10.1017/jfm.2019.740)

GitHub's **"Cite this repository"** button (from [`CITATION.cff`](CITATION.cff)) gives the BibTeX entry.

## Contents

| Folder | Content |
|---|---|
| [`src/`](src) | `navier_slip.f` (the slip module, included by the `.usr` file) and `apply_robin_patch.py` (the only change to Nek5000: about 30 lines added to `core/subs1.f`) |
| [`cylinder2d/`](cylinder2d) | complete 2D cylinder case (mesh, `SIZE`, `.par`, `.usr`, build and run scripts). **Its [README](cylinder2d/README.md) explains the method in detail.** |
| [`verification/`](verification) | exact-solution tests on an annulus (Stokes, decaying Navier–Stokes mode, slip Couette, in 2D and 3D): case file, run generator, plotting scripts, summary data and figures |
| [`legendre2009/`](legendre2009) | the comparison with Legendre *et al.* (2009): 122 runs from Re = 5.5 to 800, summary tables, digitised reference curves, figures and scripts; plus an independent re-run of the Re = 20 sweep (October 2026) |

## Quick start

```bash
cd cylinder2d
./build.sh /path/to/Nek5000     # patches core/subs1.f once (idempotent), then makenek
./run.sh 4                      # Re = 20, Kn = λ/a = 0.1; about 10 minutes on 4 cores
```

The run stops by itself at steady state. The last line of `cyl_drag.dat` should give Cd = 1.9201.

Requirements: Nek5000 (tested with master, commit ff73775, August 2026), gfortran and an MPI library
(OpenMPI), Python 3 with NumPy for the mesh generator and Matplotlib/SciPy for the plots. Use a
Nek5000 clone dedicated to slip cases, since one core routine is patched.

## Validation in brief

**Exact solutions** (`verification/`):
- spectral (exponential) convergence down to about 1e−10 at polynomial order N = 13 for no-slip, slip
  and shear-free walls;
- second-order time accuracy with BDF2, third-order with BDF3;
- the slip-Couette annulus, whose solution contains the curvature term, is matched to about 1e−6.

**Legendre, Lauga & Magnaudet (2009)** (`legendre2009/`):
- Re = 20: Cd(0) = 2.0314 and Cd(∞) = 1.3286; the normalised drag Cd\*(Kn) agrees with their
  figure 2a to within 0.002 (0.005 at Kn = 1).
- Separation boundary (their figure 1): within 2–4 % at all Re, including the bend at high Re.
- Shear-free drag: equal to 3 digits at all Re.
- Onset of vortex shedding: our critical slip is **lower** than theirs.

  | Re | 50 | 100 | 200 | 500 | 800 |
  |---|---|---|---|---|---|
  | Kn_c, this code | 0.052 | 0.33 | 0.42 | 0.36 | 0.31 |
  | Kn_c, Legendre *et al.* | 0.08 | 0.36 | 0.47 | 0.50 | 0.52 |

  We found no hysteresis, and our Re = 800 result is mesh-converged. Their no-slip Cd, CL and St at
  Re ≥ 500 are 4–18 % below ours, which may point to resolution on their side, but the discrepancy
  is not settled. Details in `legendre2009/data/legendre_campaign.txt`.

## What is not here

The raw time series of the Legendre campaign (about 120 MB of `cyl_drag.dat` files) and the flow
fields are not in the repository; the summary tables and figures are.

## License

BSD 3-Clause, see [`LICENSE`](LICENSE). Nek5000 itself is distributed under its own license.

## References

- D. Legendre, E. Lauga and J. Magnaudet, Influence of slip on the dynamics of two-dimensional wakes,
  J. Fluid Mech. 633, 437–447 (2009), arXiv:0905.0648.
- S. C. R. Dennis and G.-Z. Chang, Numerical solutions for steady flow past a circular cylinder at
  Reynolds numbers up to 100, J. Fluid Mech. 42, 471–489 (1970).
