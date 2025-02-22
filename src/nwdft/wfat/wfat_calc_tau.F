c=======================================================
c=======================================================
      subroutine wfat_calc_tau0(n, m, npt, r, th, ph, erg, zq, tau0)

      implicit none
#include "errquit.fh"
      integer n                     ! [input], 
      integer m                     ! [input], 
      integer npt                   ! [input], 
      double precision r(npt)       ! [input], r coordinate in a.u.
      double precision th(npt)      ! [input], theta coordinate in radians.
      double precision ph(npt)      ! [input], phi coordinate in radians.
      double precision erg          ! [input], 
      double precision zq           ! [input], 
      double complex tau0(npt)      ! [output],

      integer i
      double precision kappa, a, b, g, d, eta, r0, xi, phi0, z
      double precision beta0, pi_const, factorial, mm
      double complex zz, chgm_, imnum, czero
      external factorial
      parameter ( pi_const = 3.14159265359d0 )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      

      mm = dble(m)
      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))
      

c     ==== Calculate tau_\nu^{(0)} ====
      do i = 1, npt
c        == R_\nu^{(0)} ===
         eta = r(i) * (1 - cos(th(i)))
         z = kappa * eta
         zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
         call wfat_cchg(a, b, zz, chgm_)
         if (m /= 0) then
            r0 = d * kappa**(b/2) * eta**(abs(mm)/2) * exp(-z/2)
     &           * real(chgm_)
         else
            r0 = d * kappa**(b/2) * exp(-z/2) * real(chgm_)
         endif

c        == phi_\nu^{(0)} ==
         xi = r(i) * (1 + cos(th(i)))         
         call wfat_calc_phi0xi(xi, n, m, erg, phi0)

         tau0(i) = r0 * phi0 * exp(imnum*m*ph(i)) / sqrt(2*pi_const)
         if ( isnan(abs( tau0(i) )) ) tau0(i) = czero
      enddo

      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
      tau0 = -2*kappa**(beta0/kappa-1) * tau0
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_tau0_lexp(n, m, npt, r, th, ph, erg, zq,
     &           lmax, tau0)

      implicit none
#include "stdio.fh"      
#include "errquit.fh"
#include "mafdecls.fh"
      integer n                     ! [input],
      integer m                     ! [input],
      integer npt                   ! [input],
      double precision r(npt)       ! [input], r coordinate in a.u.
      double precision th(npt)      ! [input], theta coordinate in radians.
      double precision ph(npt)      ! [input], phi coordinate in radians.
      double precision erg          ! [input],
      double precision zq           ! [input],
      integer lmax                  ! [input],
      double complex tau0(npt)      ! [output],

      integer i, nsh, nl, l, i_rf, i_sh
      integer l_sh, k_sh, l_hf, k_hf, l_rf, k_rf
      double precision kappa, lmb, a, b
      double complex xx

      
      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0

      
c     ==== Allocations ====
      nsh = (lmax + 1)**2
      if (.not. MA_alloc_get(MT_DCPL, nsh, 'sphm', l_sh, k_sh))
     &     call errquit('WFAT: sphm allocation failed', 59, MA_ERR)
      nl = lmax - abs(m) + 1
      if (.not. MA_alloc_get(MT_DBL, nl, 'chgm', l_hf, k_hf))
     &     call errquit('WFAT: chgm allocation failed', 29, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl, 'radf', l_rf, k_rf))
     &     call errquit('WFAT: radf allocation failed', 29, MA_ERR)
      

c     ==== Calculate tau_\nu^{(0)} ====
      do i = 1, npt
         call wfat_chgm_1pt(r(i), abs(m), lmax, a, b, lmb, dbl_mb(k_hf))
         call wfat_radf_1pt(n, m, int(zq), lmax, kappa, r(i), 
     &        dbl_mb(k_hf), dbl_mb(k_rf))
         call wfat_sphm_1pt(lmax, th(i), ph(i), dcpl_mb(k_sh))
         
         tau0(i) = 0.0d0
         do l = abs(m), lmax
            i_rf = k_rf + l - abs(m)
            i_sh = k_sh + (l**2 + l + m + 1) - 1
            xx = dbl_mb(i_rf) * dcpl_mb(i_sh)
            tau0(i) = tau0(i) + xx
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sh))
     &     call errquit('WFAT: sphm deallocation failed', 71, MA_ERR)
      if (.not. MA_free_heap(l_hf))
     &     call errquit('WFAT: chgm deallocation failed', 71, MA_ERR)
      if (.not. MA_free_heap(l_rf))
     &     call errquit('WFAT: radf deallocation failed', 71, MA_ERR)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_calc_tau1(n, m, npt, r, th, ph, erg, zq, bt, gm,
     &           chmmt, tau1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"      
      integer n                     ! [input], 
      integer m                     ! [input], 
      integer npt                   ! [input], 
      double precision r(npt)       ! [input], r coordinate in a.u.
      double precision th(npt)      ! [input], theta coordinate in radians.
      double precision ph(npt)      ! [input], phi coordinate in radians.
      double precision erg          ! [input], 
      double precision zq           ! [input], 
      double precision bt           ! [input], beta orientation angle in degrees. 
      double precision gm           ! [input], gamma orientation angle in degrees.
      double precision chmmt(3)     ! [input], 
      double complex tau1(npt)      ! [output], 

      integer i, j, nlgl, l_r1, k_r1, l_et, k_et, l_es, k_es, l_id, k_id
      double precision kappa, a, b, g, d, eta, r0, r1, xi, phi0, phi1, z
      double precision beta0, pi_const, factorial, mm, sumphi1, sumr1
      double complex zz, chgm_, imnum, czero
      external factorial
      parameter ( pi_const=3.14159265359d0, nlgl=8 )
      parameter ( imnum=cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


c      if (npt /= 1)
c     &     call errquit('WFAT: For the moment, the argument npt of '//
c     &     'wfat_calc_tau1 must be equal to one', 27, UERR)
      mm = dble(m)
      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))


c     ==== Obtain eta ====
      if (.not. MA_alloc_get(MT_DBL, npt, 'r1', l_r1, k_r1))
     &     call errquit('WFAT: r1 allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'eta', l_et, k_et))
     &     call errquit('WFAT: eta allocation failed', 84, MA_ERR)
      do i = 1, npt
         dbl_mb(k_et+i-1) = r(i) * (1 - cos(th(i)))
      enddo
      if (.not. MA_alloc_get(MT_DBL, npt, 'e_sorted', l_es, k_es))
     &     call errquit('WFAT: e_sorted allocation failed', 29, MA_ERR)
      if (.not. MA_alloc_get(MT_INT, npt, 'indx', l_id, k_id))
     &     call errquit('WFAT: indx allocation failed', 29, MA_ERR)
      dbl_mb(k_es:k_es+npt-1) = dbl_mb(k_et:k_et+npt-1)
      do i = 1, npt
         int_mb(k_id+i-1) = i
      enddo
      call wfat_dsortx(dbl_mb(k_es), 1, npt, int_mb(k_id))
cc      if (ga_nodeid()==0)  then
cc         do i = 1, npt
cc            write(luout, '("eta, id : ", es21.12, i10)')
cc     &           dbl_mb(k_es+i-1), int_mb(k_id+i-1)
cc         enddo
cc      endif


c     ==== Calculate R_\nu^{(1)} ====
c      call wfat_r1(npt, dbl_mb(k_et), n, m, erg, zq, bt, gm, chmmt,
c     &     dbl_mb(k_r1))
      call wfat_r1(npt, dbl_mb(k_es), n, m, erg, zq, bt, gm, chmmt,
     &     nlgl, dbl_mb(k_r1))
      dbl_mb(k_es:k_es+npt-1) = dbl_mb(k_r1:k_r1+npt-1)
      do i = 1, npt
         j = int_mb(k_id+i-1)
         dbl_mb(k_r1+j-1) = dbl_mb(k_es+i-1)
      enddo
cfordebug      write(luout, '("node, sum r1: ", i3, es20.10)')
cfordebug     &     ga_nodeid(), sum(dbl_mb(k_r1:k_r1+npt-1))


c     ==== Calculate tau_\nu^{(1)} ====
cfordebug      sumr1 = 0.0d0
cfordebug      sumphi1 = 0.0d0
      do i = 1, npt
c        == R_\nu^{(0)} / \sqrt(\eta) ==
         eta = dbl_mb(k_et+i-1)
         z = kappa * eta
         zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
         call wfat_cchg(a, b, zz, chgm_)
         if (m /= 0) then
            r0 = d * kappa**(b/2) * eta**(abs(mm)/2) * exp(-z/2)
     &           * real(chgm_)
         else
            r0 = d * kappa**(b/2) * exp(-z/2) * real(chgm_)
         endif
         if (isnan(r0) .or. abs(r0) > 1.0d20) r0 = 0.0d0

c        == R_\nu^{(1)} / \sqrt(\eta) ==
         if (eta > 1.0d-12) then
            r1 = dbl_mb(k_r1+i-1) / sqrt(eta)
         else
            r1 = 0.0d0         ! 2)
         endif
         if (isnan(r1) .or. abs(r1) > 1.0d20) r1 = 0.0d0
c        2) r1 is R^{(1)}/\sqrt{\eta} where R^{(1)} is defined in eq. (A1a) of Ref. [1], 
c        and satisfies the asymptotic property of eq. (A2a). Thus, when \eta approaches
c        0, the division by \sqrt{\eta} would still make r1 approaches 0.
         
cfordebug         if (.not. isnan(r1)) sumr1 = sumr1 + r1

         
c        == phi_\nu^{(0)} and phi_\nu^{(1)} ==
         xi = r(i) * (1 + cos(th(i)))
         call wfat_calc_phi0xi(xi, n, m, erg, phi0)
         call wfat_calc_phi1xi(xi, n, m, erg, bt, gm, chmmt, phi1)
         
cfordebug         sumphi1 = sumphi1 + phi1

         tau1(i) = (r0 * phi1 + r1 * phi0) * exp(imnum*m*ph(i)) /
     &             sqrt(2*pi_const)
         if ( isnan(abs( tau1(i) )) ) tau1(i) = czero
      enddo
      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
      tau1 = -2*kappa**(beta0/kappa-1) * tau1

cfordebug      write(luout, '(" sum Phi1 = ", i3, f10.4, 4es20.10)')
cfordebug     &     ga_nodeid(), bt, sumphi1, sumr1, sum(tau1)
      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_es))
     &     call errquit('WFAT: e_sorted deallocation failed',
     &     46, MA_ERR)
      if (.not. MA_free_heap(l_id))
     &     call errquit('WFAT: indx deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_r1))
     &     call errquit('WFAT: r1 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_et))
     &     call errquit('WFAT: eta deallocation failed', 46, MA_ERR)

      end
c     References:
c     [1] PHYSICAL REVIEW A 93, 033412 (2016).
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_tau1_lexp(n, m, npt, r, th, ph, erg, zq,
     &           bt, gm, chmmt, lmax, tau1)

      implicit none
#include "stdio.fh"
#include "global.fh"      
#include "errquit.fh"
#include "mafdecls.fh"
      integer n                     ! [input], 
      integer m                     ! [input], 
      integer npt                   ! [input], 
      double precision r(npt)       ! [input], r coordinate in a.u.
      double precision th(npt)      ! [input], theta coordinate in radians.
      double precision ph(npt)      ! [input], phi coordinate in radians.
      double precision erg          ! [input], 
      double precision zq           ! [input], 
      double precision bt           ! [input], beta orientation angle in degrees. 
      double precision gm           ! [input], gamma orientation angle in degrees.
      double precision chmmt(3)     ! [input], 
      integer lmax                  ! [input], 
      double complex tau1(npt)      ! [output],

      integer nlgl, nsh, nl, i, j, l, i_q1, i_q2, i_sh
      integer l_rs, l_id, l_sh, l_q1, l_q2
      integer k_rs, k_id, k_sh, k_q1, k_q2
      double precision mu(3), protmat(3,3)
      double complex xx
      parameter(nlgl=8)


cc      if (npt /= 1)
cc     &     call errquit('WFAT: For the moment, the argument npt of '//
cc     &     'wfat_calc_tau1_lexp must be equal to one', 27, UERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'r_sorted', l_rs, k_rs))
     &     call errquit('WFAT: r_sorted allocation failed', 29, MA_ERR)
      if (.not. MA_alloc_get(MT_INT, npt, 'indx', l_id, k_id))
     &     call errquit('WFAT: indx allocation failed', 29, MA_ERR)
      dbl_mb(k_rs:k_rs+npt-1) = r
      do i = 1, npt
         int_mb(k_id+i-1) = i
      enddo
      call wfat_dsortx(dbl_mb(k_rs), 1, npt, int_mb(k_id))
cc      if (ga_nodeid()==0)  then
cc         do i = 1, npt
cc            write(luout, '("r, id : ", es21.12, i10)')
cc     &           dbl_mb(k_rs+i-1), int_mb(k_id+i-1)
cc         enddo
cc      endif

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      mu = matmul(transpose(protmat), -chmmt)


c     ==== Calculate Q_l's ====
      nsh = (lmax + 1)**2
      if (.not. MA_alloc_get(MT_DCPL, nsh, 'sphm', l_sh, k_sh))
     &     call errquit('WFAT: sphm allocation failed', 59, MA_ERR)
      nl = lmax - abs(m) + 1
      if (.not. MA_alloc_get(MT_DBL, nl*npt, 'q1', l_q1, k_q1))
     &     call errquit('WFAT: q1 allocation failed', 29, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl*npt, 'q2', l_q2, k_q2))
     &     call errquit('WFAT: q2 allocation failed', 29, MA_ERR)
c      call wfat_calc_q(npt, r, n, m, lmax, nlgl, erg, zq,
c     &     dbl_mb(k_q1), dbl_mb(k_q2))
      call wfat_calc_q(npt, dbl_mb(k_rs), n, m, lmax, nlgl, erg, zq,
     &     dbl_mb(k_q1), dbl_mb(k_q2))
      call wfat_reorder_q(lmax, m, npt, int_mb(k_id),
     &     dbl_mb(k_q1), dbl_mb(k_q2))
      
      
c     ==== Calculate tau_\nu^{(1)} ====
      do i = 1, npt
         call wfat_sphm_1pt(lmax, th(i), ph(i), dcpl_mb(k_sh))

         tau1(i) = 0.0d0
         do l = abs(m), lmax
            i_q1 = k_q1 + (i-1)*nl + l - abs(m)
            i_q2 = k_q2 + (i-1)*nl + l - abs(m)
            i_sh = k_sh + (l**2 + l + m + 1) - 1
            xx = (dbl_mb(i_q1) + mu(3)*dbl_mb(i_q2)) * dcpl_mb(i_sh)
            tau1(i) = tau1(i) + xx
         enddo
      enddo
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_rs))
     &     call errquit('WFAT: r_sorted deallocation failed',
     &     46, MA_ERR)
      if (.not. MA_free_heap(l_id))
     &     call errquit('WFAT: indx deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_sh))
     &     call errquit('WFAT: sphm deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_q1))
     &     call errquit('WFAT: q1 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_q2))
     &     call errquit('WFAT: q2 deallocation failed', 46, MA_ERR)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_collect_tau0(geom, nch, nqn, chpos, erg, nid, 
     &     mid, bt, gm, nptp, zi, v_m, tau0)

      implicit none
#include "stdio.fh"
      integer geom
      integer nch
      integer nqn
      double precision chpos(3,nch)
      double precision erg(nch)     ! [input],
      integer nid(nqn)              ! [input], 
      integer mid(nqn)              ! [input], 
      double precision bt           ! [input], beta orientation angle in degrees.
      double precision gm           ! [input], gamma orientation angle in degrees.
      integer nptp                  ! [input], the number of quadrature points for the current process.
      double precision zi           ! [input],
      double precision v_m(3,nptp)  ! in molecular frame.
      double complex tau0(nptp, nqn, nch)

      logical t1
      integer i, j, k, n, m
      double precision r, t, p, v_f(3), v_r(3), protmat(3,3)

      
      call wfat_passrot(bt, gm, protmat)
      
      do i = 1, nch
         do j = 1, nqn
            n = nid(j)
            m = mid(j)

            do k = 1, nptp
c              v_m = grid point in molecular frame
c              v_f = grid point in unrotated channel-center frame (the axes of the
c                    unrotated frame are parallel to those of the molecular frame)
c              v_r = grid point in rotated channel-center frame
               call geom_cent2arb(geom, 0, v_m(:,k), chpos(:,i), v_f,
     &              t1)
               v_r = matmul(transpose(protmat), v_f)
               r = sqrt(sum(v_r*v_r))
               t = acos(v_r(3)/r)
               p = atan(v_r(2), v_r(1))
         
               call wfat_calc_tau0(n, m, 1, r, t, p, erg(i), zi, 
     &              tau0(k,j,i))
            enddo
            
         enddo
      enddo
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_collect_tau1(geom, nch, nqn, chmmt, chpos, erg,
     &     nid, mid, bt, gm, nptp, zi, v_m, wfat0, tau1)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                  ! [input],
      integer nch                   ! [input],
      integer nqn                   ! [input],
      double precision chmmt(3,nch) ! [input],
      double precision chpos(3,nch) ! [input],
      double precision erg(nch)     ! [input],
      integer nid(nqn)              ! [input],
      integer mid(nqn)              ! [input],
      double precision bt           ! [input], beta orientation angle in degrees.
      double precision gm           ! [input], gamma orientation angle in degrees.
      integer nptp                  ! [input], the number of quadrature points for the current process.
      double precision zi           ! [input],
      double precision v_m(3,nptp)  ! [input], in molecular frame.
      logical wfat0(nqn,nch)        ! [input],
      double complex tau1(nptp, nqn, nch)    ! [output],

      logical t1
      integer l_r, l_t, l_p
      integer k_r, k_t, k_p
      integer i, j, k, n, m
      double precision v_f(3), v_r(3), protmat(3,3)
      double complex czero

      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


      call wfat_passrot(bt, gm, protmat)

c     ==== Coordinates in molecular frame ====
      if (.not. MA_alloc_get(MT_DBL, nptp, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 't', l_t, k_t))
     &     call errquit('WFAT: t allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 'p', l_p, k_p))
     &     call errquit('WFAT: p allocation failed', 59, MA_ERR)

      do i = 1, nch
         do k = 1, nptp
c           v_m = grid point in molecular frame
c           v_f = grid point in unrotated channel-center frame (the axes of the
c                 unrotated frame are parallel to those of the molecular frame)
c           v_r = grid point in rotated channel-center frame
            call geom_cent2arb(geom, 0, v_m(:,k), chpos(:,i), v_f, t1)
            v_r = matmul(transpose(protmat), v_f)
            dbl_mb(k_r+k-1) = sqrt(sum(v_r*v_r))
            dbl_mb(k_t+k-1) = acos(v_r(3)/dbl_mb(k_r+k-1))
            dbl_mb(k_p+k-1) = atan(v_r(2), v_r(1))
         enddo
         
         do j = 1, nqn
            n = nid(j)
            m = mid(j)
            if (wfat0(j,i)) then
               tau1(:,j,i) = czero
            else
               call wfat_calc_tau1(n, m, nptp, dbl_mb(k_r), dbl_mb(k_t),
     &              dbl_mb(k_p), erg(i), zi, bt, gm, chmmt(:,i),
     &              tau1(:,j,i))
            endif
         enddo
      enddo
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_t))
     &     call errquit('WFAT: t deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_p))
     &     call errquit('WFAT: p deallocation failed', 61, MA_ERR)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_reorder_q(lmax, m, nr, indx, q_1, q_2)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer lmax              ! [input], 
      integer m                                   ! [input], 
      integer nr                                  ! [input], 
      integer indx(nr)                            ! [input], 
      double precision q_1(abs(m):lmax, nr)       ! [in/output], orientation-independent
      double precision q_2(abs(m):lmax, nr)       ! [in/output], orientation-dependent

      integer nl, nn, l_tm, k_tm, i, j, i_tm1, i_tm2
      

      nl = (lmax - abs(m) + 1)
      nn = nl * nr
      if (.not. MA_alloc_get(MT_DBL, nn, 'temp', l_tm, k_tm))
     &     call errquit('WFAT: temp allocation failed', 29, MA_ERR)


      do i = 1, nr
         i_tm1 = k_tm + (i-1)*nl
         i_tm2 = k_tm + (i-1)*nl + lmax - abs(m)
         dbl_mb(i_tm1:i_tm2) = q_1(:,i)
      enddo
         
      do i = 1, nr
         j = indx(i)
         i_tm1 = k_tm + (i-1)*nl
         i_tm2 = k_tm + (i-1)*nl + lmax - abs(m)
         q_1(:,j) = dbl_mb(i_tm1:i_tm2)
      enddo
      

      do i = 1, nr
         i_tm1 = k_tm + (i-1)*nl
         i_tm2 = k_tm + (i-1)*nl + lmax - abs(m)
         dbl_mb(i_tm1:i_tm2) = q_2(:,i)
      enddo
      
      do i = 1, nr
         j = indx(i)
         i_tm1 = k_tm + (i-1)*nl
         i_tm2 = k_tm + (i-1)*nl + lmax - abs(m)
         q_2(:,j) = dbl_mb(i_tm1:i_tm2)
      enddo
      

      if (.not. MA_free_heap(l_tm))
     &     call errquit('WFAT: temp deallocation failed', 46, MA_ERR)

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_collect_tau0_OLD(geom, nch, nqn, erg, chpos, nid, 
     &     mid, bt, gm, natp, nptp, centid, ptlim, gq, zi, tau0)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer nch
      integer nqn
      double precision erg(nch)     ! [input],
      double precision chpos(3,nch)
      integer nid(nqn)              ! [input], 
      integer mid(nqn)              ! [input], 
      double precision bt           ! [input], beta orientation angle in degrees.
      double precision gm           ! [input], gamma orientation angle in degrees.
      integer natp                  ! [input], the number of atoms for the current process.
      integer nptp                  ! [input], the number of quadrature points for the current process.
      integer centid(natp)          ! [input],
      integer ptlim(2,natp)         ! [input],
      integer gq                    ! [input],
      double precision zi           ! [input],
      double complex tau0(nqn, nch, nptp)

      integer l_r, l_t, l_p
      integer k_r, k_t, k_p
      integer i, j, n, m


c     ==== Coordinates in molecular frame ====
      if (.not. MA_alloc_get(MT_DBL, nptp, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 't', l_t, k_t))
     &     call errquit('WFAT: t allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 'p', l_p, k_p))
     &     call errquit('WFAT: p allocation failed', 59, MA_ERR)

      
      do i = 1, nch
c        ==== Get the grid points with the current channel center ====
         call wfat_get_gridpt_rot(geom, natp, centid, ptlim, chpos(:,i),
     &        gq, bt, gm, dbl_mb(k_r), dbl_mb(k_t), dbl_mb(k_p))

         do j = 1, nqn
            n = nid(j)
            m = mid(j)
            call wfat_calc_tau0(n, m, nptp, dbl_mb(k_r), dbl_mb(k_t),
     &           dbl_mb(k_p), erg(i), zi, tau0(j,i,:))
         enddo
      enddo
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_t))
     &     call errquit('WFAT: t deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_p))
     &     call errquit('WFAT: p deallocation failed', 61, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_collect_tau1_OLD(geom, nch, nqn, erg, chpos, nid,
     &     mid, bt, gm, natp, nptp, centid, ptlim, gq, zi, chmmt, tau1)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer nch
      integer nqn
      double precision erg(nch)     ! [input],
      double precision chpos(3,nch)
      integer nid(nqn)              ! [input], 
      integer mid(nqn)              ! [input], 
      double precision bt           ! [input], beta orientation angle in degrees.
      double precision gm           ! [input], gamma orientation angle in degrees.
      integer natp                  ! [input], the number of atoms for the current process.
      integer nptp                  ! [input], the number of quadrature points for the current process.
      integer centid(natp)          ! [input],
      integer ptlim(2,natp)         ! [input],
      integer gq                    ! [input],
      double precision zi           ! [input],
      double precision chmmt(3,nch) ! [input],
      double complex tau1(nqn, nch, nptp)

      integer l_r, l_t, l_p
      integer k_r, k_t, k_p
      integer i, j, n, m
      

c     ==== Coordinates in molecular frame ====
      if (.not. MA_alloc_get(MT_DBL, nptp, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 't', l_t, k_t))
     &     call errquit('WFAT: t allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 'p', l_p, k_p))
     &     call errquit('WFAT: p allocation failed', 59, MA_ERR)

      
      do i = 1, nch
c        ==== Get the grid points with the current channel center ====
         call wfat_get_gridpt_rot(geom, natp, centid, ptlim, chpos(:,i),
     &        gq, bt, gm, dbl_mb(k_r), dbl_mb(k_t), dbl_mb(k_p))

         do j = 1, nqn
            n = nid(j)
            m = mid(j)
c            call wfat_calc_tau0(n, m, nptp, dbl_mb(k_r), dbl_mb(k_t),
c     &           dbl_mb(k_p), erg(i), zi, tau0(j,i,:))

            call wfat_calc_tau1(n, m, nptp, dbl_mb(k_r), dbl_mb(k_t),
     &           dbl_mb(k_p), erg(i), zi, bt, gm, chmmt(:,i),
     &           tau1(j,i,:))
         enddo
      enddo
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_t))
     &     call errquit('WFAT: t deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_p))
     &     call errquit('WFAT: p deallocation failed', 61, MA_ERR)

      end
c=======================================================
      
