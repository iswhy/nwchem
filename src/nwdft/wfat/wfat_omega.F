c=======================================================
c=======================================================
      subroutine wfat_omega_1pt_c(icent, geom, nch, v_m, erg, chpos, 
     &     cgrid, nqn, nid, mid, lmax_, zq, omega)

      implicit none
      integer icent                                    ! [input]
      integer geom                                     ! [input]
      integer nch                                      ! [input]
      double precision v_m(3,nch)                      ! [input]
      double precision erg(nch)                        ! [input]
      double precision chpos(3,nch)                    ! [input]
      logical cgrid(nch)                               ! [input]
      integer nqn                                      ! [input]
      integer nid(nqn)                                 ! [input]
      integer mid(nqn)                                 ! [input]
      integer lmax_                                    ! [input]
      integer zq                                       ! [input]
      double complex omega(1:(1+lmax_)**2, nqn, nch)   ! [output]

      integer i
      double complex czero

      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      

      if (icent /= -1) then
         call wfat_omega_1pt(geom, v_m(:,1), nch, erg, chpos, nqn, 
     &        nid, mid, lmax_, zq, omega)
      else
         do i = 1, nch
            if (cgrid(i)) then
               call wfat_omega_1pt(geom, v_m(:,i), 1, erg(i), 
     &              chpos(:,i), nqn, nid, mid, lmax_, zq,
     &              omega(:,:,i))
            else
               omega(:,:,i) = czero
            endif
         enddo
      endif

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_omega_1pt(geom, v_m, nch, erg, chpos, nqn,
     &                          nid, mid, lmax_, zq, omega)

      implicit none
#include "errquit.fh"
      integer geom                                     ! [input]
      double precision v_m(3)                          ! [input]
      integer nch                                      ! [input]
      double precision erg(nch)                        ! [input]
      double precision chpos(3,nch)                    ! [input]
      integer nqn                                      ! [input]
      integer nid(nqn)                                 ! [input]
      integer mid(nqn)                                 ! [input]
      integer lmax_                                    ! [input]
      integer zq                                       ! [input]
      double complex omega(1:(1+lmax_)**2, nqn, nch)   ! [output]

      logical t1
      integer imf, i
      parameter (imf = 0)
      double precision v_f(3)


      if (nch <= 0)
     &     call errquit('WFAT: Argument nch of wfat_omega_1pt '//
     &                  'must be larger than or equal to 1', 1, UERR)
      if (nqn <= 0)
     &     call errquit('WFAT: Argument nqn of wfat_omega_1pt '//
     &                  'must be larger than or equal to 1', 2, UERR)
      if (lmax_ < 0)
     &     call errquit('WFAT: Argument lmax_ of wfat_omega_1pt '//
     &                  'cannot be negative', 3, UERR)
      if (minval(nid) < 0)
     &     call errquit('WFAT: Element of argument array nid '//
     &                  'of wfat_omega_1pt cannot be negative', 4, UERR)
      

      do i = 1, nch
c        ==== Calculate Omega^(n,m)_(l,m') functions ====
         call geom_cent2arb(geom, imf, v_m, chpos(:,i), v_f, t1)
         call wfat_omega_1ptx(v_f, erg(i), nqn, nid, mid, lmax_, zq,
     &                        omega(:,:,i))
      enddo
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_omega_1ptx(vp, erg, nqn, nid, mid,
     &                          lmax_, zq, omega)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"      
      integer nqn, nid(1:nqn), mid(1:nqn), lmax_, zq   ! [input]
      double precision vp(3), erg                      ! [input]
      double complex omega(1:(1+lmax_)**2, 1:nqn)      ! [output]

      integer iq, ih0, ihn, ir0, irn, ish, i_af
      integer nsh, nl
      integer lmin_, m_, l, m
      integer l_af, l_hf, l_rf
      integer k_af, k_hf, k_rf
      double precision r, th, ph, kappa, lmb, a, b
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


c      write(luout, '("ZERO here = ", 2f23.14)') real(czero), imag(czero)

      if (nqn <= 0)
     &     call errquit('WFAT: Argument nqn of wfat_omega_1ptx '//
     &                  'must be larger than or equal to 1', 2, UERR)
      if (lmax_ < 0)
     &     call errquit('WFAT: Argument lmax_ of wfat_omega_1ptx '//
     &                  'cannot be negative', 3, UERR)
      if (minval(nid) < 0)
     &     call errquit('WFAT: Element of argument array nid '//
     &                  'of wfat_omega_1ptx cannot be negative', 4,
     &                  UERR)
      if (maxval(abs(mid)) > lmax_) then
         call errquit('WFAT: Argument lmax_ of wfat_omega_1ptx'//
     &                'must be larger than maxval(abs(mid)).', 20, UERR)
      end if
      omega = czero
      
      
c     ==== Convert Cartesian inputs to spherical ones ====
      r = sqrt(sum(vp*vp))
      th = acos(vp(3)/r)
      if (th >= 0.0d0 .and. th < 1.0d-12) then
         ph = 0.0d0
      else
         ph = atan2(vp(2), vp(1))
      endif
      
c     ==== Calculate angular functions (spherical harmonics) ====
      nsh = (lmax_ + 1)**2
      if (.not. MA_alloc_get(MT_DCPL, nsh, 'angf', l_af, k_af))
     &     call errquit('WFAT: angf allocation failed', 59, MA_ERR)
      call wfat_sphm_1pt(lmax_, th, ph, dcpl_mb(k_af:k_af+nsh-1))
c      write(*, '("      th, ph, ang = ", e14.5, e14.5, 3x, 2e14.5)')
c     &     th, ph, sum(dcpl_mb(k_af:k_af+nsh-1))
c      write(*, '("      th, ph, ang = ", e14.5, e14.5, 3x, e14.5)')
c     &     th, ph, sum(abs(imag(dcpl_mb(k_af:k_af+nsh-1))))
      
c     ==== Calculate confluent hypergeometric functions ====
      lmin_ = minval(abs(mid))
      nl = lmax_ - lmin_ + 1
      if (.not. MA_alloc_get(MT_DBL, nl, 'chgf', l_hf, k_hf))
     &     call errquit('WFAT: chgf allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl, 'radf', l_rf, k_rf))
     &     call errquit('WFAT: radf allocation failed', 59, MA_ERR)
      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0
      call wfat_chgm_1pt(r, lmin_, lmax_, a, b, lmb, dbl_mb(k_hf))
      
      do iq = 1, nqn
c        ==== Calculate radial functions ====
         m_ = mid(iq)
         ih0 = k_hf - lmin_ + abs(m_)
         ihn = k_hf - lmin_ + lmax_
         ir0 = k_rf - lmin_ + abs(m_)
         irn = k_rf - lmin_ + lmax_         
         dbl_mb(k_rf : k_rf+nl-1) = 0.0d0
         call wfat_radf_1pt(nid(iq), m_, zq, lmax_, kappa, r,
     &                      dbl_mb(ih0:ihn), dbl_mb(ir0:irn))
      
c        ==== Calculate channel functions ====
         do l = abs(m_), lmax_
            do m = -l, l
               ish = l**2 + l + m + 1
               i_af = k_af - 1 + ish
               omega(ish, iq) = dbl_mb(k_rf - lmin_ + l) *
     &                          dcpl_mb(i_af)

               if ( isnan(abs( omega(ish,iq) )) )
     &              omega(ish,iq) = czero
            enddo
         enddo
         
      enddo
c      write(*, '("      r, chg, rad = ", e13.5, 3x, e13.5, e13.5)')
c     &     r, sum(dbl_mb(k_hf : k_hf + nl - 1)),
c     &     sum(dbl_mb(k_rf : k_rf + nl - 1))
            
      
      if (.not. MA_free_heap(l_af))
     &     call errquit('WFAT: angf deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_hf))
     &     call errquit('WFAT: chgf deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_rf))
     &     call errquit('WFAT: radf deallocation failed', 61, MA_ERR)

      end
c=======================================================

      
c=======================================================      
c=======================================================      
      subroutine wfat_sphm_1pt(lmax, th, ph, sphm)
      
      implicit none      
      integer lmax              ! [input]
      double precision th, ph   ! [input], theta and phi angles in radians.
      double complex sphm(1:*)  ! [output]
      
      integer nl, inow, ipast1, ipast2, imin, l, m, i
      double precision pi_const, rpart, ipart
      double complex imnum, scterm

      parameter (pi_const = 3.14159265359d0)
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )
      
      
      nl = (lmax + 1)**2
      
c     ==== (l,m) = (0,0) ====
      inow = 1
      sphm(inow) = 0.5d0/sqrt(pi_const)
      if (lmax == 0) return
      
c     ==== (l,m) = (1,-1) ====  
      inow = 2
      sphm(inow) = 0.5d0 * sqrt(1.5d0/pi_const) * sin(th) *
     &             exp(-imnum * ph)
      
c     ==== (l,m) = (1,0) ====  
      inow = 3
      sphm(inow) = 0.5d0 * sqrt(3/pi_const) * cos(th)
      
c     ==== (l,m) = (1,+1) ====  
      inow = 4
      sphm(inow) = -1 * conjg(sphm(2))
      if (lmax == 1) return
      
c     ==== l >= 2 ====
      do l = 2, lmax
         do m = 0, l
            inow = l**2 + l + m + 1
            imin = l**2 + l - m + 1        
            if (m == l) then
               ipast1 = (l-1)**2 + (l-1) + (l-1) + 1
               sphm(inow) = - sqrt((2*l+1) / dble(2*l)) * 
     &                      sin(th) * exp(imnum * ph) * 
     &                      sphm(ipast1)
               sphm(imin) = (-1)**l * conjg(sphm(inow))
            else if (m > 0 .and. m < l) then
               ipast1 = (l-1)**2 + (l-1) + m + 1
               ipast2 = (l-2)**2 + (l-2) + m + 1
               if (l-2 < abs(m)) scterm = 0.0d0
               if (l-2 >= abs(m)) scterm = sqrt((l+m-1) * dble(l-m-1)) *
     &                                     sphm(ipast2)
               sphm(inow) = sqrt((2*l+1)/dble(2*l-3)) *
     &                      ( sqrt(dble(4*(l-1)**2 - 1)) * cos(th) *
     &                      sphm(ipast1) - scterm ) / sqrt((l+m) *
     &                      dble(l-m))
               sphm(imin) = (-1)**m * conjg(sphm(inow))
            else if (m == 0) then
               ipast1 = (l-1)**2 + (l-1) + 0 + 1
               ipast2 = (l-2)**2 + (l-2) + 0 + 1
               sphm(inow) = sqrt((2*l+1)/dble(2*l-3)) *
     &                      ( sqrt(dble(4*(l-1)**2 - 1)) * cos(th) *
     &                      sphm(ipast1) - (l-1) * sphm(ipast2) ) / l
            endif
         enddo
      enddo

      i = 0
      do l = 0, lmax
         do m = -l, l
            i = i + 1
            rpart = real(sphm(i))
            ipart = imag(sphm(i))
            if (abs(rpart) <= 1.0d-12) rpart = 0.0d0
            if (abs(ipart) <= 1.0d-12) ipart = 0.0d0
            sphm(i) = cmplx(rpart, ipart, kind=kind(1.0d0))
         enddo
      enddo
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_radf_1pt(n, m, znuc, ll, kp, r, chgf, radf)
      
      !!    z = kappa * r
      implicit none
      integer n, m, ll, znuc                    ! [input]
      double precision kp, r, chgf(abs(m):ll)   ! [input]
      double precision radf(abs(m):ll)          ! [output]

      integer l
      double precision cc

      
      do l = abs(m), ll
         call wfat_calc_radnorm(l, n, m, znuc, kp, cc)
         radf(l) = cc * (kp*r)**l * exp(-kp*r) * chgf(l)


cdebug         radf(l) = exp(kp*r)
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_radnorm(l, n, m, znuc, kp, cc)
      
      implicit none
#include "errquit.fh"      
#include "stdio.fh"
      integer l, n, m, znuc     ! [input]  
      double precision kp      ! [input]
      double precision cc      ! [output]
      
      integer a, pa, pb, i
      double precision b, c, pc, d, e, f, ff, xx
      double precision factorial
      external factorial
      
      
      if (abs(m) > l) then
         call errquit('WFAT: The argument m and l of '//
     &                'wfat_calc_radnorm must be such that abs(m)<=l',
     &                61, UERR)
      end if
      
      
      pa = l + (abs(m)-m)/2 + 1
      a = (-1)**pa
      
      pb = 2*l + 3
      b = sqrt(2.0d0)**pb
      
      pc = znuc/kp - (abs(m)+1)/2.0d0 - n
      c = kp**pc
      
      d = sqrt( dble(2*l+1) * factorial(l+m) * factorial(l-m) *
     &    factorial(abs(m)+n) * factorial(n) )
      
      e = factorial(l)/factorial(2*l+1)
      
      f = 0.0d0
      do i = 0, min(n, l-abs(m))
         xx = l + 1 - znuc/kp + n - i
         call wfat_gammaf(xx, ff)

         ff = ff / ( factorial(i) * factorial(l-i) * 
     &        factorial(abs(m)+i) * factorial(l-abs(m)-i) * 
     &        factorial(n-i) )
         f = f + ff
      end do
      
      cc = a * b * c * d * e * f

c      write(luout, '("cc = ", es15.6)') cc
      
      end
c=======================================================
