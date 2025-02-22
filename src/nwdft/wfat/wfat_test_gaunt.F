!==================================================
!==================================================
      subroutine wfat_test_gaunt(l1, m1, l2, m2)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer l1
      integer m1
      integer l2
      integer m2

      integer n, i, l3, l_min, l_gn, k_gn, l_gy, k_gy
      double precision pi_const
      parameter (pi_const=3.14159265359d0)


      write(luout, '("")')
      write(luout, '("Testing Gaunt''s integrals:")')


      l_min = max(abs(l1-l2), abs(m1+m2))
      n = (l1+l2) - l_min + 1
      if (.not. MA_alloc_get(MT_DBL, n, 'gaunt_p', l_gn, k_gn))
     &     call errquit('WFAT: gaunt_p allocation failed', 59, MA_ERR)
      call wfat_gauntp_lscan(l1, m1, l2, m2, dbl_mb(k_gn))
      if (.not. MA_alloc_get(MT_DBL, n, 'gaunt_y', l_gy, k_gy))
     &     call errquit('WFAT: gaunt_y allocation failed', 59, MA_ERR)
      call wfat_gaunty_lscan(l1, m1, l2, m2, dbl_mb(k_gy))

      
      write(luout, '("l1 = ", i4)') l1
      write(luout, '("m1 = ", i4)') m1
      write(luout, '("l2 = ", i4)') l2
      write(luout, '("m2 = ", i4)') m2
      write(luout, '("m3 = ", i4)') m1+m2
      write(luout, '(a5, 2a20)') 'l3', 'Gaunt P_lm', 'Gaunt Y_lm'
      do i = 1, n
         l3 = l_min + i - 1
         write(luout, '(i5, 2f20.8)')
     &        l3, dbl_mb(k_gn+i-1), dbl_mb(k_gy+i-1)
      enddo


      if (.not. MA_free_heap(l_gn))
     &     call errquit('WFAT: gaunt_p deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_gy))
     &     call errquit('WFAT: gaunt_y deallocation failed', 61, MA_ERR)
      
      end
!==================================================
