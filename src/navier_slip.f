c=======================================================================
c  navier_slip.f -- Navier (partial) slip on curved walls in Nek5000
c=======================================================================
c  On the wall faces with boundary ID bid, with n the unit normal
c  pointing into the fluid and S the rate-of-strain tensor:
c        u.n = 0 ,        u_t = 2*lambda*(S.n)_t ,
c  i.e. the tangential traction exerted on the fluid is
c        t_t = -(mu/lambda) u_t .
c     lambda > 0 : Navier slip  ('shl' faces + implicit Robin term)
c     lambda = 0 : no slip      ('W  ' faces)
c     lambda < 0 : shear free   ('shl' faces, no Robin term)
c
c  How it works (details in README.md): the 'shl' condition gives
c  u.n = 0 and a natural (traction) condition on the wall.  The Navier
c  traction enters the weak form as the boundary term
c        (mu/lambda) int_wall u.v dS ,
c  which is added implicitly to the velocity Helmholtz operator: it is
c  stored as rbc = (mu/lambda)*area on the wall GLL nodes (nslip_coef)
c  and applied as au <- au + rbc*u by axhm_robin, the routine that
c  apply_robin_patch.py adds to the Nek5000 core.
c
c  Requirements
c   - stressFormulation = yes (PN/PN-2: lx2 = lx1-2, and lx1m = lx1 in
c     SIZE).  Only this formulation imposes u.n = 0 on faces that are
c     not aligned with the axes, and its natural boundary term is the
c     full traction 2 mu S.n, so the curvature term of (S.n)_t is
c     included automatically.
c   - The core patch (apply_robin_patch.py) applied to core/subs1.f.
c   - userbc must return trn = tr1 = tr2 = 0 on the 'shl' faces.
c   - Constant viscosity and an isotropic slip length.
c
c  Usage: put the line
c        include 'navier_slip.f'
c  after the last routine of the .usr file, and call
c     usrdat  : call nslip_init
c     usrdat2 : call nslip_faces('W  ',bid,lambda,nface)
c     usrdat3 : call nslip_coef(bid,lambda)
c     userchk : call nslip_diag(bid,lambda,utmax,resmax)   (optional)
c-----------------------------------------------------------------------
      subroutine nslip_init
c     Switch the Robin term off and clear its coefficient (call in usrdat)
      include 'SIZE'
      common /robinbc/ rbc(lx1,ly1,lz1,lelv)
      common /robinbl/ ifrobin
      logical ifrobin

      ifrobin = .false.
      call rzero(rbc,lx1*ly1*lz1*lelv)

      return
      end
c-----------------------------------------------------------------------
      subroutine nslip_faces(cbwall,bid,slip,nface)
c     Tag the faces whose velocity BC is cbwall with boundaryID = bid and
c     give them the BC that matches the slip length.  Call in usrdat2,
c     i.e. before Nek builds the velocity masks.
      include 'SIZE'
      include 'TOTAL'
      character*3 cbwall
      integer bid,nface,e,f
      real slip

      nface = 0
      do e=1,nelv
      do f=1,2*ldim
         if (cbc(f,e,1).eq.cbwall) then
            boundaryID(f,e) = bid
            nface = nface + 1
            if (slip.eq.0.0) then
               cbc(f,e,1) = 'W  '
            else
               cbc(f,e,1) = 'shl'
            endif
         endif
      enddo
      enddo
      nface = iglsum(nface,1)
      if (nio.eq.0) write(6,'(a,i3,a,i7,a,1pe13.5)') ' nslip: bid',bid,
     $   '  faces',nface,'  slip length',slip

      return
      end
c-----------------------------------------------------------------------
      subroutine nslip_coef(bid,slip)
c     Robin coefficient rbc = (mu/lambda)*area on the nodes of the 'shl'
c     faces of boundary ID bid.  Call in usrdat3 (geometry final).  AREA
c     holds the face quadrature weights (surface Jacobian times GLL
c     weights), so that sum rbc*u*v = (mu/lambda) int u.v dS; a node
c     shared by two wall faces receives both contributions, as the
c     assembled (dssum-ed) boundary integral requires.  Face points are
c     visited in Nek's face ordering (dsset/skpdat, as in drgtrq).
      include 'SIZE'
      include 'TOTAL'
      common /robinbc/ rbc(lx1,ly1,lz1,lelv)
      common /robinbl/ ifrobin
      logical ifrobin
      integer bid,e,f,pf,i
      real slip,visc,coef,wmeas

      if (slip.le.0.0) return
      visc = cpfld(1,1)                       ! dynamic viscosity
      if (visc.lt.0.0) visc = -1.0/visc
      coef  = visc/slip
      wmeas = 0.0
      call dsset(lx1,ly1,lz1)
      do e=1,nelv
      do f=1,2*ldim
         if (boundaryID(f,e).eq.bid .and. cbc(f,e,1).eq.'shl') then
            pf     = eface1(f)
            js1    = skpdat(1,pf)
            jf1    = skpdat(2,pf)
            jskip1 = skpdat(3,pf)
            js2    = skpdat(4,pf)
            jf2    = skpdat(5,pf)
            jskip2 = skpdat(6,pf)
            i = 0
            do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
               i = i + 1
               rbc(j1,j2,1,e) = rbc(j1,j2,1,e) + coef*area(i,1,f,e)
               wmeas = wmeas + area(i,1,f,e)
            enddo
            enddo
         endif
      enddo
      enddo
      wmeas   = glsum(wmeas,1)
      ifrobin = .true.
      if (nio.eq.0) write(6,'(a,i3,a,1p2e15.7)') ' nslip: bid',bid,
     $   '  Robin term on; mu/lambda, wall measure =',coef,wmeas

      return
      end
c-----------------------------------------------------------------------
      subroutine nslip_diag(bid,slip,utmax,resmax)
c     Wall diagnostics on the faces of boundary ID bid.  With n the Nek
c     normal (out of the fluid) and tau = 2 S.n:
c        utmax  = max |u_t|
c        resmax = max |u_t + lambda tau_t|, the residual of the Navier
c                 condition (max |u_t| for no slip, max |tau_t| for a
c                 shear-free wall)
c     Also fills common /nslipw/: nslut = tangential wall velocity
c     (2D: u.t with t = (-n_y, n_x); 3D: |u_t|) and nslwm = 1 on the
c     wall nodes of bid, 0 elsewhere.
      include 'SIZE'
      include 'TOTAL'
      common /nslipw/ nslut(lx1*ly1*lz1*lelv),nslwm(lx1*ly1*lz1*lelv)
      real nslut,nslwm
      common /nslips/ nssij(lx1*ly1*lz1*6*lelv)
      real nssij
      real nsur(lx1*ly1*lz1),nsus(lx1*ly1*lz1),nsut(lx1*ly1*lz1)
     $   , nsvr(lx1*ly1*lz1),nsvs(lx1*ly1*lz1),nsvt(lx1*ly1*lz1)
     $   , nswr(lx1*ly1*lz1),nsws(lx1*ly1*lz1),nswt(lx1*ly1*lz1)
      integer bid,e,f,pf,i,k,ijk,iv,ib,nxyz,nij
      real slip,utmax,resmax
      real nsn(3),nsu(3),nst(3),nss(3,3),nsr(3),nsut3(3),nstt(3)
      real nsun,nstn,nsrm,nsum

      nxyz = lx1*ly1*lz1
      nij  = 3
      if (if3d) nij = 6
      call comp_sij(nssij,nij,vx,vy,vz,nsur,nsus,nsut,nsvr,nsvs,nsvt,
     $              nswr,nsws,nswt)
      call rzero(nslut,nxyz*nelv)
      call rzero(nslwm,nxyz*nelv)
      utmax  = 0.0
      resmax = 0.0
      call dsset(lx1,ly1,lz1)
      do e=1,nelv
      do f=1,2*ldim
         if (boundaryID(f,e).eq.bid) then
            pf     = eface1(f)
            js1    = skpdat(1,pf)
            jf1    = skpdat(2,pf)
            jskip1 = skpdat(3,pf)
            js2    = skpdat(4,pf)
            jf2    = skpdat(5,pf)
            jskip2 = skpdat(6,pf)
            i = 0
            do j2=js2,jf2,jskip2
            do j1=js1,jf1,jskip1
               i   = i + 1
               ijk = j1 + lx1*(j2-1)
               iv  = ijk + nxyz*(e-1)
               ib  = ijk + nxyz*nij*(e-1)
               nsn(1)   = unx(i,1,f,e)
               nsn(2)   = uny(i,1,f,e)
               nsu(1)   = vx(ijk,1,1,e)
               nsu(2)   = vy(ijk,1,1,e)
               nss(1,1) = nssij(ib)
               nss(2,2) = nssij(ib+nxyz)
               if (if3d) then
                  nsn(3)   = unz(i,1,f,e)
                  nsu(3)   = vz(ijk,1,1,e)
                  nss(3,3) = nssij(ib+2*nxyz)
                  nss(1,2) = nssij(ib+3*nxyz)
                  nss(2,3) = nssij(ib+4*nxyz)
                  nss(1,3) = nssij(ib+5*nxyz)
               else
                  nsn(3)   = 0.0
                  nsu(3)   = 0.0
                  nss(3,3) = 0.0
                  nss(1,2) = nssij(ib+2*nxyz)
                  nss(2,3) = 0.0
                  nss(1,3) = 0.0
               endif
               nss(2,1) = nss(1,2)
               nss(3,2) = nss(2,3)
               nss(3,1) = nss(1,3)
               nsun = 0.0
               do k=1,3
                  nst(k) = nss(k,1)*nsn(1) + nss(k,2)*nsn(2)
     $                   + nss(k,3)*nsn(3)
                  nsun   = nsun + nsu(k)*nsn(k)
               enddo
               nstn = nst(1)*nsn(1) + nst(2)*nsn(2) + nst(3)*nsn(3)
               nsrm = 0.0
               nsum = 0.0
               do k=1,3
                  nsut3(k) = nsu(k) - nsun*nsn(k)
                  nstt(k)  = nst(k) - nstn*nsn(k)
                  if (slip.gt.0.0) then
                     nsr(k) = nsut3(k) + slip*nstt(k)
                  elseif (slip.lt.0.0) then
                     nsr(k) = nstt(k)
                  else
                     nsr(k) = nsut3(k)
                  endif
                  nsrm = nsrm + nsr(k)**2
                  nsum = nsum + nsut3(k)**2
               enddo
               resmax = max(resmax,sqrt(nsrm))
               utmax  = max(utmax ,sqrt(nsum))
               if (if3d) then
                  nslut(iv) = sqrt(nsum)
               else
                  nslut(iv) = -nsu(1)*nsn(2) + nsu(2)*nsn(1)
               endif
               nslwm(iv) = 1.0
            enddo
            enddo
         endif
      enddo
      enddo
      utmax  = glmax(utmax,1)
      resmax = glmax(resmax,1)

      return
      end
c-----------------------------------------------------------------------
