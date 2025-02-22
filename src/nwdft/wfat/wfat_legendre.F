c=======================================================
c=======================================================
      double precision function assocleg_norm(l, m)

      implicit none
      integer l
      integer m

      double precision factorial
      external factorial

      assocleg_norm = 2.0d0 * factorial(l+m) / (2*l+1) / factorial(l-m)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_assocleg_fixm(m, lmax, x, pl)

      implicit none
#include "stdio.fh"
      integer m
      integer lmax
      double precision x
      double precision pl(abs(m):lmax)

      integer l
      double precision l_, pl_1, pl_2, double_fact
      external double_fact

      
c     l = |m|, m = |m|
      l = abs(m)
      pl(l) = (-1)**l * double_fact(2*l-1) * (1-x**2)**(dble(l)/2)
c      write(luout, '("P_m = ", f16.8)') pl(l)
      
c     l = |m| + 1, m = |m|
      l = abs(m) + 1
      pl(l) = x * (2*(l-1)+1) * pl(l-1)
c      write(luout, '("P_m+1 = ", f16.8)') pl(l)
      
c     If m < 0, convert the above first two functions with positive m.
      if (m < 0) then
         call wfat_assocleg_convert_negm(m, l, pl(abs(m)))
         call wfat_assocleg_convert_negm(m, l, pl(abs(m)+1))
      endif

c     Obtain the functions with l > |m|+1
      pl_2 = pl(abs(m))
      pl_1 = pl(abs(m)+1)
      do l = abs(m)+2, lmax
         l_ = dble(l)
         pl(l) = ((2*l_-1)*x*pl_1 - (l_-1+m)*pl_2) / (l_-m)
         pl_2 = pl_1
         pl_1 = pl(l)
      enddo

c      write(luout, '("P_l:")')
c      call wfat_dprintmat(lmax-abs(m)+1, 1, 'f', pl)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_assocleg_convert_negm(m, l, pl)
c     DESCRIPTION:
c     Convert from P_{l,m}(x) to P_{l,-m}(x).
      implicit none
      integer m
      integer l
      double precision pl   ![in/output]

      double precision factorial
      external factorial

      pl = (-1)**m * factorial(l-m)/factorial(l+m) * pl

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_assocleg_quad(a, b, m, n, pt, wt)
c     DESCRIPTION:
c     This subroutine computes the quadrature nodes and weights of the
c     Gauss-associated-Legendre quadrature. The nodes are calculated as
c     the zeros of P_{m+n,m}(x).
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      double precision a         ![input], left interval boundary.
      double precision b         ![input], right interval boundary.
      integer m                  ![input], the m parameter in P_{l,m}(x) and m >= 0.
      integer n                  ![input], the no. of quadrature nodes.
      double precision pt(n)     ![output], quadrature nodes.
      double precision wt(n)     ![output], quadrature weights.

      integer i, j, j_, l_dm, k_dm, l_nd, k_nd, i_dm, i_dm_
      double precision wt0, factorial
      external factorial
      
      
      if (m < 0)
     &     call errquit('WFAT: The input m of wfat_assocleg_quad '//
     &     'cannot be negative. The quadrature nodes and weights '//
     &     'for a negative value of m are identical to those of '//
     &     'the corresponding positive m. So, pass the absolute '//
     &     'value instead.', 42, UERR)


c     ==== Construct the dipole moment matrix ====
      if (.not. MA_alloc_get(MT_DBL, n**2, 'DM', l_dm, k_dm))
     &     call errquit('WFAT: DM allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, n, 'node', l_nd, k_nd))
     &     call errquit('WFAT: node allocation failed', 59, MA_ERR)
      dbl_mb(k_dm : k_dm+n**2-1) = 0.0d0
      do j = 1, n
         do i = 1, n
            if (j == i + 1) then
               j_ = m + i
               i_dm = k_dm + (j-1)*n + (i-1)
               dbl_mb(i_dm) = sqrt((j_**2 - m**2)/(4*j_**2 - 1.0d0))   ! See eq. (3.24) to (3.26) in Ref. [1].
               i_dm_ = k_dm + (i-1)*n + (j-1)
               dbl_mb(i_dm_) = dbl_mb(i_dm)
            endif
         enddo
      enddo
      call wfat_dsyev_wrap(n, dbl_mb(k_dm), dbl_mb(k_nd))

      
c     ==== Calculate the shifted nodes and weights ====
      do i = 1, n
         i_dm = k_dm + (i-1)*n          ! The first row of the i-th column.
         wt0 = 2**(2*m+1) * (factorial(m))**2/factorial(2*m+1) *
     &        (dbl_mb(i_dm))**2                                     ! eq. (3.45) in Ref. [1].
         pt(i) = (b-a)/2 * (dbl_mb(k_nd+i-1) + (a+b)/(b-a))         ! The shifted node.
         wt(i) = (b-a)/2 * wt0                                      ! The shifted weight.
      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dm))
     &     call errquit('WFAT: Cannot deallocate DM', 32, MA_ERR)
      if (.not. ma_free_heap(l_nd))
     &     call errquit('WFAT: Cannot deallocate node', 32, MA_ERR)
      
      end
c     REFERENCES:
c     [1] Spectral Methods - Algorithm, Analysis and Applications by Shen et al.
c=======================================================
