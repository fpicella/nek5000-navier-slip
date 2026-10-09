# Navier slip on a curved wall in Nek5000: the 2D circular cylinder

This folder is a complete, small Nek5000 case: two-dimensional flow past a circular cylinder whose
wall obeys a Navier (partial-slip) condition with a prescribed slip length λ. The slip is exact on
the curved wall (curvature term included) and implicit in time, so any slip length can be used with
the usual time step.

The price is a small addition to one Nek5000 core routine: about 30 lines, applied by a script.
Everything else lives in the user files.

## Quick start

```bash
./build.sh /path/to/Nek5000     # patches core/subs1.f once, then makenek
./run.sh 4                      # 4 MPI ranks; about 10 minutes on a desktop
```

The default case is Re = 20 with Kn = λ/a = 0.1. It stops by itself when the flow is steady. The
last line of `cyl_drag.dat` should give Cd = 1.9201. See the reference values below.

Tested with Nek5000 master (commit ff73775, August 2026), gfortran and OpenMPI.

## 1. The boundary condition

The cylinder has radius a = D/2 and sits in a uniform stream U; Re = UD/ν. On the wall, with n the
unit normal pointing into the fluid, the Navier condition reads

```
u · n = 0                      (no penetration)
u_t   = 2 λ (S · n)_t          (slip proportional to the tangential shear)
```

Here S = (∇u + ∇uᵀ)/2 is the rate-of-strain tensor and the index t denotes the tangential part.

- λ = 0 is the no-slip wall and λ → ∞ the shear-free wall.
- The slip is measured by Kn = λ/a.
- On a flat wall y = 0 the condition reduces to the familiar u = λ ∂u/∂y.

It is useful to write the same condition as a tangential traction. With σ = −p I + 2μS and n_out
the normal pointing out of the fluid (Nek's convention), the traction exerted by the wall on the
fluid is

```
t_t = (σ · n_out)_t = −(μ/λ) u_t
```

This is a friction proportional to the slip velocity: the larger λ, the weaker the friction.

On the cylinder, in polar coordinates with u_r = 0 at r = a:

```
2 (S · n)_t = ∂u_θ/∂r − u_θ/a      so      u_θ = λ (∂u_θ/∂r − u_θ/a)
```

The term −u_θ/a is the curvature term. It is part of the frame-invariant condition (Legendre,
Lauga & Magnaudet 2009, eq. 2.2). A flat-wall formula u_t = λ ∂u_t/∂n would miss it.

## 2. How it is imposed in Nek5000

### 2.1 Where the slip enters the weak form

Nek5000 solves the momentum equation in weak form. With the stress formulation, the viscous term
tested against a velocity test function v is

```
∫_Ω 2μ S(u) : S(v) dΩ  −  ∫_Γ (σ · n_out) · v dΓ
```

On the slip wall both u and v satisfy u·n = v·n = 0, so only the tangential traction survives in the
boundary integral. Replacing it with the Navier traction gives

```
− ∫_wall t_t · v dΓ  =  + (μ/λ) ∫_wall u · v dΓ
```

This is a symmetric, positive "boundary mass" term: a Robin condition. It belongs on the left-hand
side, next to the viscous operator, and makes the discrete problem better conditioned, not worse.

### 2.2 Two ingredients of Nek5000 that are used as they are

- **The stress formulation** (`stressFormulation = yes`). It is needed for two reasons:
  - Nek can impose u·n = 0 on a wall that is not aligned with the axes only in this formulation. It
    removes the normal velocity component with the local face normal, node by node.
  - The natural boundary term of this formulation is the full traction 2μS·n. So the curvature term
    of section 1 comes out of the discretisation automatically, with nothing to add by hand.

  The formulation requires PN–PN−2 (`lx2 = lx1−2`) and `lx1m = lx1` in `SIZE`.
- **The `shl` boundary condition** on the cylinder faces: zero normal velocity plus a natural
  (traction) condition on the tangential components. `userbc` returns a zero traction.

### 2.3 The trick: the slip friction goes into the operator, implicitly

In the discrete weak form, the boundary term of section 2.1 is a diagonal matrix R that acts on the
velocity at the wall nodes:

```
(R u)_i = rbc_i u_i ,      rbc_i = (μ/λ) w_i
```

Here w_i is the quadrature weight of wall node i, that is the surface Jacobian times the GLL
weights. Nek stores these weights in the `area` array. The coefficient `rbc` is built once, in
`usrdat3`, by `nslip_coef`.

- A node shared by two wall faces gets the sum of both weights, which is what the assembled integral
  requires.
- The `shl` mask already removes the normal component, so the same isotropic coefficient acts only
  on the tangential velocity.

Each time step, Nek solves a Helmholtz problem for the velocity, H u = f, with a conjugate-gradient
solver. In the stress formulation the matrix–vector product with H is the core routine `axhmsf`
(file `core/subs1.f`). The patch adds a routine, `axhm_robin`, that does

```
au ← au + rbc · u         (each velocity component)
```

and calls it at the end of `axhmsf`. The friction is then treated implicitly: it is part of the
matrix, at the new time level, exactly like the viscous term. There is no time-step restriction
from the slip length, and λ → 0 recovers the no-slip wall smoothly.

**Why not impose it explicitly, without touching the core?** The obvious alternative sets the
traction t_t = −(μ/λ)u_t in `userbc`, from the velocity of the previous step. That treats a very stiff
friction explicitly. It is stable only for a time step of order

```
dt ≲ λ h / (μ N²)
```

where h is the size of the wall elements and N the polynomial order. For small slip lengths that
is orders of magnitude below the time step the flow itself needs. The implicit term has no such
limit, at the cost of a few lines in the core.

The patch is written by `apply_robin_patch.py`:
- It appends `axhm_robin` to `core/subs1.f` and inserts the call to it in `axhmsf`.
- It is idempotent, and it does nothing unless the user switches the term on (`ifrobin`, set by
  `nslip_coef`). An unmodified case runs exactly as before.

The coefficient lives in a common block, `/robinbc/`, filled by the user files. The CG solver's
Jacobi preconditioner does not include the new term. At very small slip lengths (Kn ≈ 1e−4) this
costs a few percent more iterations, and nothing at moderate Kn.

### 2.4 Checking that the condition holds

`nslip_diag` evaluates, at every wall node, the tangential velocity and the tangential shear from
the computed field. It reports the largest residual of the Navier condition,
max |u_t − 2λ(S·n)_t| (the `res_navier` column of `cyl_drag.dat`). In the default case it is about
1e−7, against slip velocities of about 0.3.

## 3. Files

| File | Content |
|---|---|
| `cyl.usr` | the case: boundary and initial conditions, forces, diagnostics, steady-state stop |
| `navier_slip.f` | the slip condition (included by `cyl.usr`): faces, Robin coefficient, diagnostics |
| `apply_robin_patch.py` | the only change to Nek5000: adds `axhm_robin` to `core/subs1.f` |
| `SIZE` | 2D, polynomial order 7 (`lx1 = 8`), PN–PN−2, `lx1m = lx1` |
| `cyl.par` | Re = 20, Kn = 0.1, BDF2 with OIFS, dt = 5e−3 |
| `cyl.re2`, `cyl.ma2`, `cyl.co2` | the mesh, ready to use (1040 elements) |
| `gen_cyl_re2.py` | the mesh generator (writes `cyl.re2` directly) |
| `build.sh`, `run.sh` | patch and compile; run |

The mesh is a polar O-grid. It has 40 sectors and 26 rings, stretched geometrically from a first
ring 0.08 D thick to an outer circle at r = 40 a. The rings are exact circular arcs. The outer
circle is an inflow (`v`) for x < 0 and an outflow (`O`) for x > 0; the cylinder faces are `W` in
the file and are switched to `shl` at run time.

## 4. Parameters and output

In `cyl.par`:
- `viscosity = −Re`.
- `userParam01` = Kn = λ/a: > 0 slip, 0 no slip, < 0 shear free.
- `userParam02`: stop when max |∂u/∂t| drops below this value; 0 runs to `endTime`.
- `userParam03`: interval of the diagnostics, in steps.
- `userParam06`: amplitude of a small transverse-velocity blob behind the cylinder at t = 0, to
  start the shedding at Re > 47 (e.g. 0.1).

`cyl_drag.dat` has one line per diagnostic, with the columns
`istep time Cd Cd_p Cd_v Cl us_max us_min res_navier dudt_max`:
- Cd = 2Fx, split into its pressure (Cd_p) and viscous (Cd_v) parts; Cl = 2Fy.
- us is the tangential wall velocity, oriented from the front to the rear of the cylinder. us_min < 0
  means a separated flow.
- `res_navier` is the residual described in section 2.4; `dudt_max` is the steady-state monitor.

## 5. Reference results

**Steady flow at Re = 20** (this mesh, N = 7), drag coefficient against Kn:

| Kn | Cd | Cd_p | Cd_v |
|---|---|---|---|
| 0 (no slip) | 2.0314 | 1.2246 | 0.8068 |
| 0.01 | 2.0197 | 1.2200 | 0.7997 |
| 0.05 | 1.9736 | 1.1976 | 0.7759 |
| 0.1 | 1.9201 | 1.1668 | 0.7533 |
| 0.2 | 1.8308 | 1.1098 | 0.7210 |
| 0.5 | 1.6668 | 0.9985 | 0.6683 |
| 1 | 1.5447 | 0.9144 | 0.6302 |
| shear free | 1.3286 | 0.7680 | 0.5605 |

- The no-slip value compares with 2.045 of Dennis & Chang (1970).
- The ratio Cd(Kn)/Cd(0) agrees with figure 2a of Legendre, Lauga & Magnaudet (2009) to within 0.002
  (0.005 at Kn = 1).
- Refining the mesh to 1904 elements changes Cd by less than 1e−6.

**Checks against exact solutions.** A circular-Couette annulus, with slip on the inner cylinder and
the outer wall rotating, has an exact solution u_θ = A r + B/r. Its slip coefficient contains the
curvature term. The computed wall velocity and torque match it to about 1e−6. With exact
solutions, the error decreases exponentially with the polynomial order, and the time-stepping error
is of second order with BDF2.

**Unsteady flow at Re = 100** (`viscosity = −100`, `userParam06 = 0.1`, `userParam02 = 0`, run to
t ≈ 300): without slip Cd = 1.331, lift amplitude 0.323, Strouhal number 0.1649. With slip the
shedding weakens; it stops for Kn > 0.33 (linear stability of this code). Legendre et al. report
0.36–0.38.

## 6. Changing the case

- **Reynolds number:** `viscosity = −Re`. Above Re ≈ 47 the flow sheds vortices; set
  `userParam06 = 0.1`, `userParam02 = 0` and a longer `endTime`.
- **Slip length:** `userParam01` = Kn. Any value works with the same time step; `Kn < 0` gives a
  shear-free wall.
- **Mesh:** `python3 gen_cyl_re2.py --ntheta 40 --nr 26 --rinf 40 --dr1 0.08` writes `cyl.re2`
  (options shown with their defaults). Then run Nek's `genmap` (answer `cyl`, then `0.2`) to
  rebuild `cyl.ma2`; `gencon` likewise for `cyl.co2`. Check `lelg` in `SIZE` against the number of
  elements.
- **Another geometry:** `navier_slip.f` is independent of the cylinder. Give the slip faces a
  boundary ID with `nslip_faces`, fill the coefficient with `nslip_coef`, and keep the stress
  formulation. The same files work in 3D.

## 7. Limitations

- **Viscosity and slip length:** constant viscosity and an isotropic slip length (the same in every
  tangential direction). Both are easy to lift: `rbc` can vary from node to node.
- **Formulation:** stress formulation only. In this case it costs no more per time step than Nek's
  default Laplacian formulation (15.0 against 16.6 ms per step on 4 cores).
- **Core patch:** one core routine is patched. Keep a Nek5000 clone dedicated to these cases.

## References

- D. Legendre, E. Lauga and J. Magnaudet, Influence of slip on the dynamics of two-dimensional wakes,
  J. Fluid Mech. 633, 437–447 (2009), arXiv:0905.0648.
- S. C. R. Dennis and G.-Z. Chang, Numerical solutions for steady flow past a circular cylinder at
  Reynolds numbers up to 100, J. Fluid Mech. 42, 471–489 (1970).
- Nek5000 documentation: boundary conditions and the stress formulation,
  https://nek5000.github.io/NekDoc/.
