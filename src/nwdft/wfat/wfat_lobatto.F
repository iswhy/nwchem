c=======================================================
c=======================================================
      subroutine wfat_lgl(a, b, n, node, weig)
c     This routine returns Legendre-Gauss-Lobatto nodes and weights.
c     a = left end node
c     b = right end node
c     n+1 = number of nodes including the end nodes
c     node = quadrature nodes inlcuding end nodes
c     weig = quadrature weights including those associated with the end nodes.
c     For further reference, see eq. (3.180) and theorem 3.8 in Ref. [1].
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer n              ! The number of nodes and weights are N+1
      double precision a, b
      double precision node(0:n), weig(0:n)
      
      integer m, l_dm, k_dm, l_nd, k_nd, i, j, k, k_
      double precision lgd

      
      m = n - 2
      if (n < 2) then
         write(luout, '("n = ", i3)') n
         call errquit('Argument n of wfat_lgl must be bigger than '//
     &        'or equal to 2', 73, UERR)
      endif

      
c     ==== Construct the dipole moment matrix ====
      if (.not. MA_alloc_get(MT_DBL, (m+1)**2, 'DM', l_dm, k_dm))
     &     call errquit('WFAT: DM allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, m+1, 'node', l_nd, k_nd))
     &     call errquit('WFAT: node allocation failed', 59, MA_ERR)
      dbl_mb(k_dm : k_dm+(m+1)**2-1) = 0.0d0
      do i = 0,m
         do j = 0,m
            if (j == i+1) then
               k = i*(m+1) + j + 1
               dbl_mb(k_dm+k-1) = sqrt(1.0d0 * j*(j+2)/(2*j+1)/(2*j+3))   ! See eq. (3.24) to (3.26) in Ref. [1].
               k_ = j*(m+1) + i + 1
               dbl_mb(k_dm+k_-1) = dbl_mb(k_dm+k-1)
            end if
         end do
      end do
      call wfat_dsyev_wrap(m+1, dbl_mb(k_dm), dbl_mb(k_nd))


c     ==== Calculate shifted nodes and weights ====
      node(1:n-1) = dbl_mb(k_nd:k_nd+(m+1)-1)         !node0(0:M)
      node(0) = -1.0d0
      node(n) = 1.0d0
      do i = 0,n
         call wfat_legendre('n', n, node(i), lgd)
         weig(i) = 2.0d0/n/(n+1)/lgd**2              ! Calculate weights
         node(i) = 0.5d0*((b-a)*node(i) + b + a)     ! Calculate shifted nodes   
         weig(i) = (b-a)/2.0d0*weig(i)               ! Calculate shifted weights
      end do

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dm))
     &     call errquit('WFAT: Cannot deallocate DM', 49, MA_ERR)
      if (.not. ma_free_heap(l_nd))
     &     call errquit('WFAT: Cannot deallocate node', 49, MA_ERR)
      end
c     REFERENCES:
c     [1] Spectral Methods - Algorithm, Analysis and Applications by Shen et al.
c=======================================================
