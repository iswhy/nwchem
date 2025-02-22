c=======================================================
c=======================================================
      subroutine wfat_calc_phi0xi(xi, n, m, erg, phi0)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      double precision xi
      integer n
      integer m
      double precision erg
      double precision phi0

      integer nlag, l_la, k_la
      double precision kappa, s(1), factorial
      external factorial
      

      kappa = sqrt(2 * abs(erg))
      
      nlag = n + 1
      if (.not. MA_alloc_get(MT_DBL, nlag, 'lag', l_la, k_la))
     &     call errquit('WFAT: lag allocation failed', 17, MA_ERR)
      s = kappa * xi
      call wfat_lagfun(1, n, dble(abs(m)), s, dbl_mb(k_la))

      phi0 = sqrt(kappa) * sqrt( factorial(n)/factorial(n+abs(m)) ) *
     &       s(1)**(dble(abs(m))/2) * exp(-s(1)/2) * dbl_mb(k_la + n)

      if (.not. MA_free_heap(l_la))
     &     call errquit('WFAT: lag deallocation failed', 46, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_phi1xi(xi, n, m, erg, bt, gm, chmmt, phi1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      double precision xi
      integer n
      integer m
      double precision erg
      double precision bt           ! [input], beta orientation angle in degree.
      double precision gm           ! [input], gamma orientation angle in degree.
      double precision chmmt(3)
      double precision phi1

      integer nlag, i, n_, ii, l_la, k_la
      double precision protmat(3,3), mu(3), kappa, s(1), ff
      double precision phi0_i, cc(-2:1), c(-2:2), factorial
      external factorial
      

c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      mu = matmul(transpose(protmat), -chmmt)


c     ==== Obtain the required Laguerre functions ====
      kappa = sqrt(2 * abs(erg))
      nlag = (n+2) + 1
      if (.not. MA_alloc_get(MT_DBL, nlag, 'lag', l_la, k_la))
     &     call errquit('WFAT: lag allocation failed', 17, MA_ERR)
      s = kappa * xi
      call wfat_lagfun(1, n+2, dble(abs(m)), s, dbl_mb(k_la))
      

c     ==== Calculate the coefficients for phi0 ====
      cc = 0.0d0
      do i = -2, 1
         n_ = n + i
         if (n_ >= 0) then
            cc(i) = sqrt(dble(n_+1) * dble(n_+abs(m)+1))
         endif
      enddo
      c(-2) = cc(-2) * cc(-1) / 8 / kappa**3
      c(-1) = -cc(-1)/2 * (kappa*mu(3) + 2*n + abs(m))/kappa**3
      c(0) = 0.0d0
      c(1) = cc(0)/2 * (kappa*mu(3) + 2*n + abs(m) + 2)/kappa**3
      c(2) = -cc(0) * cc(1) / 8 / kappa**3


c     ==== Calculate phi1 ====
      ff = sqrt(kappa) * s(1)**(dble(abs(m))/2) * exp(-s(1)/2)
      phi1 = 0.0d0
      do i = n-2, n+2
         if (i < 0 .or. i == n) cycle

         ii = i - n           ! the values of ii are -2, -1, 0, 1, 2.
         phi0_i = ff * sqrt( factorial(i)/factorial(i+abs(m)) ) *
     &            dbl_mb(k_la + i)
         phi1 = phi1 + c(ii) * phi0_i
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_la))
     &     call errquit('WFAT: lag deallocation failed', 46, MA_ERR)

      end
c=======================================================
