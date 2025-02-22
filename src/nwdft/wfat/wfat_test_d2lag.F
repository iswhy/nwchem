c================================================
c================================================
      subroutine wfat_test_d2lag(n, m)
c     DESCRIPTION:
c     Calculate and prints the second derivative of r\phi_n^{|m|}(r)
c     where
c     \phi_n^{|m|}(x) = (kr)^{|m|/2} e^{-kr/2} L_{n-1}^{|m|}(kr)
c     numerically and exactly.
      implicit none
#include "mafdecls.fh"
#include "stdio.fh"      
#include "errquit.fh"
#include "wfat.fh"
      integer n
      integer m

      character*255 fname
      integer l_pt, l_ln, l_ln1, l_ln2
      integer k_pt, k_ln, k_ln1, k_ln2
      integer npt, i, p0, p1, p2
      double precision mm, ke, r_start, r_last, dr, x, ln, ln1, ln2
      double precision d2, d2_, x2h, x1h, x0h, x1h_, x2h_, f2h, f1h
      double precision f0h, f1h_, f2h_

      
      if (n < 1) 
     &     call errquit('WFAT: The input argument n of '//
     &     'wfat_test_d2lag cannot be smaller than 1', 40, UERR)

      npt = 401
      if (.not. MA_alloc_get(MT_DBL, npt, 'xpoints', l_pt, k_pt))
     &     call errquit('WFAT: xpoints allocation failed', 74, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt*(n+1), 'L_n', l_ln, k_ln))
     &     call errquit('WFAT: L_n allocation failed', 74, MA_ERR)

c     ==== Pointers to L_n^m, L_{n-1}^m, and L_{n-2}^m, respectively ====
      p0 = k_ln + n*npt
      p1 = k_ln + (n-1)*npt
      p2 = k_ln + (n-2)*npt

      mm = dble(m)
      ke = 4.2d0
      r_start = 0.d0
      r_last = 10.0d0
      dr = (r_last - r_start) / (npt - 1)
      do i = 1, npt
         dbl_mb(k_pt+i-1) = ke*r_start + (i-1)*ke*dr
      enddo


c     ==== Get all generalized Laguerre polynomials from 0 to n ====
      call wfat_lagfun(npt, n, mm, dbl_mb(k_pt), dbl_mb(k_ln))

      
      call util_file_name('test.wd2l', .false., .false., fname)
      open(WF_D2LAGFUNTEST_UNIT, file=fname, action='write')
      write(WF_D2LAGFUNTEST_UNIT, '("# n = ", i3)') n
      write(WF_D2LAGFUNTEST_UNIT, '("# m = ", i3)') m
      write(WF_D2LAGFUNTEST_UNIT, '("# k = ", f12.6)') ke
      write(WF_D2LAGFUNTEST_UNIT, '("#", a4, 2a15, 2a22, a15)')
     &        'No.', 'r', 'x = kr', 'Exact', 'Numerical', 'Rel. error'
      do i = 1, npt
         x = dbl_mb(k_pt+i-1)
         ln = dbl_mb(p0+i-1)

c        ==== Depending on the value of n, some polynomials might ====
c        ====   not contribute to the exact derivative formula    ====
         if (n <= 0) then
            ln1 = 0.0d0
         else
            ln1 = dbl_mb(p1+i-1)
         endif
         if (n <= 1) then
            ln2 = 0.0d0
         else
            ln2 = dbl_mb(p2+i-1)
         endif

         
c        ==== Exact 2nd derivative ====
         d2 = ke * exp(-x/2) *
     &        (0.25d0*x**(1+mm/2)*ln1 +
     &        (n+mm/2)*n*x**(mm/2-1)*ln -
     &        (n+mm/2)**2*x**(mm/2-1)*ln1 -
     &        (n-1+mm)*(n-1)*x**(mm/2-1)*ln1 +
     &        (n-1+mm)*(n-1+mm/2)*x**(mm/2-1)*ln2)


c        ==== Numerical 2nd derivative ====
         if (i <= 2 .or. i >= npt-1) then
            x0h = dbl_mb(k_pt+(i+0)-1)
            d2_ = 0.0d0
         else
            x2h = dbl_mb(k_pt+(i+2)-1)
            x1h = dbl_mb(k_pt+(i+1)-1)
            x0h = dbl_mb(k_pt+(i+0)-1)
            x1h_ = dbl_mb(k_pt+(i-1)-1)
            x2h_ = dbl_mb(k_pt+(i-2)-1)

            f2h = x2h**(1+mm/2) * exp(-x2h/2) * dbl_mb(p1+(i+2)-1)
            f1h = x1h**(1+mm/2) * exp(-x1h/2) * dbl_mb(p1+(i+1)-1)
            f0h = x0h**(1+mm/2) * exp(-x0h/2) * dbl_mb(p1+(i+0)-1)
            f1h_ = x1h_**(1+mm/2) * exp(-x1h_/2) * dbl_mb(p1+(i-1)-1)
            f2h_ = x2h_**(1+mm/2) * exp(-x2h_/2) * dbl_mb(p1+(i-2)-1)
            d2_ = (-f2h + 16*f1h - 30*f0h + 16*f1h_ - f2h_)
     &           / (12*dr**2) / ke
         endif

         
c        ==== Print to file ====
         write(WF_D2LAGFUNTEST_UNIT, '(i5, 2es15.6, 2es22.12, f15.8)')
     &        i, x0h/ke, x0h, d2, d2_, (d2-d2_)/abs(d2)
      enddo
      close(WF_D2LAGFUNTEST_UNIT)

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate xpoints', 32, MA_ERR)
      if (.not. ma_free_heap(l_ln))
     &     call errquit('WFAT: Cannot deallocate L_n', 32, MA_ERR)

      end
c================================================
      
