c=======================================================
c=======================================================
      subroutine wfat_calc_d(n, m, l, erg, zq, d1, d2)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n
      integer m
      integer l
      double precision erg
      double precision zq
      double precision d1          ! [Output], orientation-independent
      double precision d2          ! [Output], orientation-dependent

      integer i, n_, ii
      double precision kappa, beta0, gm, gn, gnn, vv, bt1_1, e1, d0, e1_
      double precision hl1, hl3, hlt1, hlt3, f1, bt1_2, e2, e2_, hl2
      double precision hlt2, f2
      double precision cc(-2:1), c(-2:2)
      double precision factorial
      external factorial


cOLD      if (m < 0)
cOLD     &     call errquit('WFAT: For the moment, the argument m of '//
cOLD     &     'wfat_calc_d cannot be negative.', 51, UERR)
      if (l < abs(m))
     &     call errquit('WFAT: Argument l and m of wfat_calc_d '//
     &     'should be such that l >= |m|.', 84, UERR)
      

      kappa = sqrt(2*abs(erg))
      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
      call wfat_gammaf(1+abs(m)+n-zq/kappa, gm)
      call wfat_gammaf(1+n-zq/kappa, gn)
      call wfat_gammaf(-n+zq/kappa, gnn)
      vv = (-1)**(l-abs(m)+m) * sqrt(2.0d0)*2**l * kappa *
     &     sqrt((2*l+1) * factorial(l+m) * factorial(l-m)) *
     &     factorial(l) / factorial(2*l+1) * (-1)**( (abs(m)-m)/2 )

      cc = 0.0d0
      do i = -2, 1
         n_ = n + i
         if (n_ >= 0) then
            cc(i) = sqrt(dble(n_+1) * dble(n_+abs(m)+1))
         endif
      enddo
      

c     ==== Orientation-independent part ====
      bt1_1 = -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      c(-2) = cc(-2) * cc(-1) / 8 / kappa**3
      c(-1) = -cc(-1)/2 * (2*n + abs(m))/kappa**3
      c(0) = 0.0d0
      c(1) = cc(0)/2 * (2*n + abs(m) + 2)/kappa**3
      c(2) = -cc(0) * cc(1) / 8 / kappa**3

      e1 = 0.0d0
      do i = n-2, n+2
         if (i == n) cycle

         ii = i - n
         call wfat_calc_d0(n, i, m, l, erg, zq, d0)
         e1_ = sqrt(factorial(i+abs(m)) * factorial(i)) * c(ii) * d0
         e1 = e1 + e1_
      enddo
      
      call wfat_calc_hl(n, m, l, 1, erg, zq, hl1)
      call wfat_calc_hl(n, m, l, 3, erg, zq, hl3)
      call wfat_calc_hlt(n, m, l, 1, erg, zq, hlt1)
      call wfat_calc_hlt(n, m, l, 3, erg, zq, hlt3)
      f1 = sqrt( factorial(n+abs(m)) * factorial(n) ) / gm *
     &     (
     &     (-1)**(abs(m)+1)/gn *
     &     (bt1_1/kappa*hl1 + 1.0d0/4/kappa**3*hl3)
     &     - gnn * (bt1_1/kappa*hlt1 + 1.0d0/4/kappa**3*hlt3)
     &     )
      d1 = (e1 + f1) * vv
      d1 = -2*kappa**(beta0/kappa-1) * d1

      
c     ==== Orientation-dependent part ====
      bt1_2 = -1.0d0/2/kappa * (2*n+1+abs(m))
      c(-2) = 0.0d0
      c(-1) = -cc(-1)/2 / kappa**2
      c(0) = 0.0d0
      c(1) = cc(0)/2 / kappa**2
      c(2) = 0.0d0

      e2 = 0.0d0
      do i = n-1, n+1
         if (i == n) cycle

         ii = i - n
         call wfat_calc_d0(n, i, m, l, erg, zq, d0)
         e2_ = sqrt(factorial(i+abs(m)) * factorial(i)) * c(ii) * d0
         e2 = e2 + e2_
      enddo

      call wfat_calc_hl(n, m, l, 2, erg, zq, hl2)
      call wfat_calc_hlt(n, m, l, 2, erg, zq, hlt2)
      f2 = sqrt( factorial(n+abs(m)) * factorial(n) ) / gm *
     &     (
     &     (-1)**(abs(m)+1)/gn *
     &     (bt1_2/kappa*hl1 - 1.0d0/2/kappa**2*hl2)
     &     - gnn * (bt1_2/kappa*hlt1 - 1.0d0/2/kappa**2*hlt2)
     &     )
      d2 = (e2 + f2) * vv
      d2 = -2*kappa**(beta0/kappa-1) * d2

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_calc_d0(n, nn, m, l, erg, zq, d0)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n
      integer nn
      integer m
      integer l
      double precision erg
      double precision zq
      double precision d0

      integer k
      double precision kappa, dd, gk, factorial
      external factorial


      if (l < abs(m))
     &     call errquit('WFAT: Argument l and m of wfat_calc_d0 '//
     &     'should be such that l >= |m|.', 84, MA_ERR)
      

      kappa = sqrt(2*abs(erg))

      d0 = 0.0d0
      do k = 0, min(nn, l-abs(m))
         call wfat_gammaf(l+1-zq/kappa+n-k, gk)
         dd = gk / factorial(k) / factorial(l-k) / factorial(abs(m)+k) /
     &        factorial(l-abs(m)-k) / factorial(nn-k)

         d0 = d0 + dd
      enddo

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_calc_hl(n, m, l, ll, erg, zq, hl)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n
      integer m
      integer l
      integer ll
      double precision erg
      double precision zq
      double precision hl

      integer imax, i, j, j_, k, k_, s1
      double precision kappa, a, b, ai, aj, aj_, ak, ak_, hl0
      double precision gj, gk, gk_, ps1, ps2, ps3, factorial
      external factorial


      if (ll /= 1 .and. ll /= 2 .and. ll /= 3)
     &     call errquit('WFAT: Argument ll of wfat_calc_hl should be '//
     &     'either 1, 2, or 3.', 84, MA_ERR)
      if (l < abs(m))
     &     call errquit('WFAT: Argument l and m of wfat_calc_hl '//
     &     'should be such that l >= |m|.', 84, MA_ERR)
      

      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = abs(m) + 1
      

      imax = l - 2*abs(m) - ll
      hl = 0.0d0
      do i = 0, imax
         do j = 0, imax
            do j_ = 0, n
               do k = 0, imax
                  do k_ = 0, imax

                     s1 = i + j + j_ + k + k_
                     if (s1 /= l - 2*abs(m) - ll) cycle
                     
                     ai = (-1)**i / factorial(i)

                     call wfat_gammaf(a+j, gj)
                     aj = gj / factorial(abs(m)+j) / factorial(j)
                     
                     aj_ = 1.0d0 / factorial(abs(m)+j_) / 
     &                    factorial(j_) / factorial(n-j_)

                     call wfat_gammaf(a+k, gk)
                     ak = gk / factorial(abs(m)+k) / factorial(k)

                     call wfat_gammaf(a+k_, gk_)
                     ak_ = gk_ / factorial(abs(m)+k_) / factorial(k_)

                     call wfat_psif(a+k, ps1)
                     call wfat_psif(dble(k+1), ps2)
                     call wfat_psif(dble(k+b), ps3)

                     hl0 = ai * aj * aj_ * ak * ak_ *
     &                    ((k-k_)*(ps1-ps2-ps3)/dble(l-abs(m)-k-j_)/
     &                    dble(l-abs(m)-k_-j_) +
     &                    1.0d0/(dble(l-abs(m)-k_-j_))**2)
                     hl = hl + hl0
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_hlt(n, m, l, ll, erg, zq, hlt)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n
      integer m
      integer l
      integer ll
      double precision erg
      double precision zq
      double precision hlt

      integer imax, i, j, j_, k, k_, s2
      double precision kappa, a, b, ai, aj, aj_, ak, ak_, hlt0
      double precision gj, gk, gk_, factorial
      external factorial
      

      if (ll /= 1 .and. ll /= 2 .and. ll /= 3)
     &     call errquit('WFAT: Argument ll of wfat_calc_hlt should '//
     &     'be either 1, 2, or 3.', 84, MA_ERR)
      if (l < abs(m))
     &     call errquit('WFAT: Argument l and m of wfat_calc_hlt '//
     &     'should be such that l >= |m|.', 84, MA_ERR)

      
      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = abs(m) + 1
      

      imax = l - abs(m) - ll
      hlt = 0.0d0
      do i = 0, imax
         do j = 0, imax
            do j_ = 0, n
               do k = 1, abs(m)
                  do k_ = 0, imax

                     s2 = i + j + j_ - k + k_
                     if (s2 /= l - 2*abs(m) - ll) cycle

                     ai = (-1)**i / factorial(i)

                     call wfat_gammaf(a+j, gj)
                     aj = gj / factorial(abs(m)+j) / factorial(j)

                     aj_ = 1.0d0 / factorial(abs(m)+j_) / 
     &                    factorial(j_) / factorial(n-j_)

                     call wfat_gammaf(zq/kappa-n-abs(m)+k, gk)
                     ak = factorial(k-1) / factorial(abs(m)-k) / gk

                     call wfat_gammaf(a+k_, gk_)
                     ak_ = gk_ / factorial(abs(m)+k_) / factorial(k_)

                     hlt0 = ai * aj * aj_ * ak * ak_ *
     &                    dble(k+k_) / dble(l-abs(m)+k-j_) /
     &                    dble(l-abs(m)-k_-j_)
                     hlt = hlt + hlt0

                  enddo
               enddo
            enddo
         enddo
      enddo

      end
c=======================================================
      
