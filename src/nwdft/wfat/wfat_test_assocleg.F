c=======================================================
c=======================================================
      subroutine wfat_test_assocleg(m, n)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer n

      integer l_pt, k_pt, l_wt, k_wt, i
      double precision xi, intg, pi_const, d1, d2, intg_exact, a, b1, b2
      double precision yi
      parameter (pi_const=3.14159265359d0)


      if (.not. MA_alloc_get(MT_DBL, n, 'node', l_pt, k_pt))
     &     call errquit('WFAT: node allocation failed', 74, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, n, 'weight', l_wt, k_wt))
     &     call errquit('WFAT: weight allocation failed', 74, MA_ERR)
     
      
c     ==== 1st integral check ====
c     I = \int_{-1}^1 dx (1-x^2)^m exp(-x)
      write(luout, '(/, a)') "== The first integral check =="
      call wfat_assocleg_quad(-1.0d0, 1.0d0, m, n, dbl_mb(k_pt),
     &     dbl_mb(k_wt))
      write(luout, '(a)') 'Associated Legendre quadrature nodes and '//
     &     'weights'
      do i = 1, n
         write(luout, '(i5, 2es18.6)')
     &        i, dbl_mb(k_pt+i-1), dbl_mb(k_wt+i-1)
      enddo

      intg = 0.0d0
      do i = 1, n
         xi = dbl_mb(k_pt+i-1)
         intg = intg + dbl_mb(k_wt+i-1) * exp(-xi)
      enddo
      write(luout, '("1st integral:")')
      write(luout, '(2x, a12, f19.10)') 'quadrature:', intg


c     ==== 2nd integral check ====
c     I = \int_{-1}^1 dx cos(b1 x) sin(b2 x) exp(ax)
      a = -2.0d0
      b1 = 2*pi_const/4
      b2 = 5.0d0
      write(luout, '(/, a)') "== The second integral check =="
      call wfat_assocleg_quad(-1.0d0, 1.0d0, m, n, dbl_mb(k_pt),
     &     dbl_mb(k_wt))
      write(luout, '(a)') 'Associated Legendre quadrature nodes and '//
     &     'weights'
      do i = 1, n
         write(luout, '(i5, 2es18.6)')
     &        i, dbl_mb(k_pt+i-1), dbl_mb(k_wt+i-1)
      enddo

      intg = 0.0d0
      do i = 1, n
         xi = dbl_mb(k_pt+i-1)
         intg = intg + dbl_mb(k_wt+i-1) *
     &        cos(b1*xi) * sin(b2*xi) * exp(a*xi) / (1-xi**2)**m
      enddo
      
      call wfat_int_expsin(-1.0d0, 1.0d0, a, b1+b2, d1)
      call wfat_int_expsin(-1.0d0, 1.0d0, a, b1-b2, d2)
      intg_exact = (d1 - d2) / 2
      write(luout, '("2nd integral: ")')
      write(luout, '(2x, a12, f19.10)') 'quadrature:', intg
      write(luout, '(2x, a12, f19.10)') 'exact:', intg_exact


c     ==== 3rd integral check ====
c     I = \int_2^10 dx cos(b1 (x-6)) sin(b2 (x-6)) exp(a(x-6))
      a = -1.0d0
      b1 = 2*pi_const/16
      b2 = 5.0d0
      write(luout, '(/, a)') "== The third integral check =="
      call wfat_assocleg_quad(2.0d0, 10.0d0, m, n, dbl_mb(k_pt),
     &     dbl_mb(k_wt))
      write(luout, '(a)') 'Associated Legendre quadrature nodes and '//
     &     'weights'
      do i = 1, n
         write(luout, '(i5, 2es18.6)')
     &        i, dbl_mb(k_pt+i-1), dbl_mb(k_wt+i-1)
      enddo

      intg = 0.0d0
      do i = 1, n
         xi = dbl_mb(k_pt+i-1)
         yi = 2.0d0/(10-2) * (xi - (2*1 - (-1)*10)/2.0d0)
         intg = intg + dbl_mb(k_wt+i-1) *
     &        cos(b1*(xi-6)) * sin(b2*(xi-6)) * exp(a*(xi-6)) /
     &        (1-yi**2)**m
      enddo

c     Change of integral limit due to the use of u=x-6.
c     2 to 10 becomes -4 to 4.
      call wfat_int_expsin(-4.0d0, 4.0d0, a, b1+b2, d1)
      call wfat_int_expsin(-4.0d0, 4.0d0, a, b1-b2, d2)
      intg_exact = (d1 - d2) / 2
      write(luout, '("3rd integral: ")')
      write(luout, '(2x, a12, f19.10)') 'quadrature:', intg
      write(luout, '(2x, a12, f19.10)') 'exact:', intg_exact      


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate node', 32, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate weight', 32, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_int_expsin(xlow, xup, a, b, d)

      implicit none
      double precision xlow, xup
      double precision a, b
      double precision d

      double precision tlow, tup

      call wfat_int_expsin0(xlow, a, b, tlow)
      call wfat_int_expsin0(xup, a, b, tup)

      d = tup - tlow

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_int_expsin0(x, a, b, t)

      implicit none
      double precision x
      double precision a, b
      double precision t

      t = exp(a*x) * (a*sin(b*x) - b*cos(b*x)) / (a**2 + b**2)

      end
c=======================================================
      
      
      
