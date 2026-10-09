#!/usr/bin/env python3
"""Patch Nek5000 core/subs1.f: add an implicit Robin (Navier-slip) boundary
term to the stress-formulation velocity Helmholtz operator (axhmsf).

    au_i <- au_i + rbc * u_i      for i = 1..ldim

rbc(lx1,ly1,lz1,lelv) (common /robinbc/) is filled by the user with
(mu/lambda)*area at the GLL nodes of the slip ('shl') faces and is zero
elsewhere.  It is the discrete form of  (mu/lambda) int_Gamma u.v dS,
i.e. the natural-BC contribution of the Navier traction  t = -(mu/lambda) u_t.
The 'shl' mask removes the normal component (u.n = 0), so the isotropic
coefficient acts on the tangential velocity only, and the full stress
tensor of the formulation carries the curvature terms of (S.n)_t.
Usage: apply_robin_patch.py path/to/Nek5000/core/subs1.f   (idempotent)
"""
import sys

ROBIN_SUB = r"""
c-----------------------------------------------------------------------
      subroutine axhm_robin (au1,au2,au3,u1,u2,u3,nel)
c
c     Implicit Robin (Navier-slip) boundary term for the stress-form
c     velocity Helmholtz operator:   au_i <- au_i + rbc * u_i .
c     rbc (common /robinbc/) = (mu/lambda)*area on the 'shl' wall nodes,
c     zero elsewhere (set in usrdat3);  ifrobin (common /robinbl/)
c     switches the term on.  See patch/apply_robin_patch.py.
c
      include 'SIZE'
      common /robinbc/ rbc(lx1,ly1,lz1,lelv)
      common /robinbl/ ifrobin
      logical ifrobin
      real au1(1),au2(1),au3(1),u1(1),u2(1),u3(1)

      if (.not.ifrobin) return
      n = lx1*ly1*lz1*nel
      call addcol3 (au1,rbc,u1,n)
      call addcol3 (au2,rbc,u2,n)
      if (ldim.eq.3) call addcol3 (au3,rbc,u3,n)

      return
      end
c-----------------------------------------------------------------------
"""

CALL = ("      if (ifield.eq.1 .and. matmod.ge.0)               ! Navier slip\n"
        "     $   call axhm_robin (au1,au2,au3,u1,u2,u3,nel)  ! (Robin term)\n\n")

def main():
    path = sys.argv[1]
    s = open(path).read()
    if 'axhm_robin' in s:
        print('subs1.f already patched'); return
    anchor = "      taxhm=taxhm+(dnekclock()-etime1)\n\n      return\n      end\n"
    assert s.count(anchor) == 1, 'anchor not found exactly once in axhmsf'
    s = s.replace(anchor, CALL + anchor) + ROBIN_SUB
    open(path, 'w').write(s)
    print('patched', path)

if __name__ == '__main__':
    main()
