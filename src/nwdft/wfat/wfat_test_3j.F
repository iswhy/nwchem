!==================================================
!==================================================
      subroutine wfat_test_3j(tj1, tm1, tj2, tm2)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer tj1
      integer tm1
      integer tj2
      integer tm2

      integer i, n, l_cf, k_cf
      double precision j1, m1, j2, m2, j_min, j3, norm

      
      write(luout, '("Testing 3-j symbols:")')
      j1 = tj1 / 2.0d0
      j2 = tj2 / 2.0d0
      m1 = tm1 / 2.0d0
      m2 = tm2 / 2.0d0

      j_min = max(abs(j1-j2), abs(m1+m2))
      if ( idnint(2*j_min) == abs(tj1-tj2) ) then
         write(luout, '("Init. type : |j1-j2|")')
      else if ( idnint(2*j_min) == abs(tm1+tm2) ) then
         write(luout, '("Init. type : |m1+m2|")')
      endif
      n = idnint((j1+j2) - j_min + 1)
      if (.not. MA_alloc_get(MT_DBL, n, '3j', l_cf, k_cf))
     &     call errquit('WFAT: 3j allocation failed', 59, MA_ERR)
      call wfat_3j_jscan(tj1, tm1, tj2, tm2, dbl_mb(k_cf))


      write(luout, '("j1 = ", f7.1)') j1
      write(luout, '("m1 = ", f7.1)') m1
      write(luout, '("j2 = ", f7.1)') j2
      write(luout, '("m2 = ", f7.1)') m2
      write(luout, '("m3 = ", f7.1)') -m1-m2
      norm = 0.0d0
      do i = 1, n
         j3 = j_min + i - 1
         write(luout, '(f7.1, f15.8)') j3, dbl_mb(k_cf+i-1)
         norm = norm + (2*j3+1)*dbl_mb(k_cf+i-1)**2
      enddo
      write(luout, '("norm^2 = ", f12.8)') norm
      write(luout, '("")')
      
      
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: 3j deallocation failed', 61, MA_ERR)
      
      end
!==================================================


!==================================================
!==================================================
      subroutine wfat_test_3j000(tj1, tj2)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer tj1
      integer tj2

      integer i, n, l_cf, k_cf, l_cf0, k_cf0
      double precision j1, j2, j3, norm


      if (mod(tj1,2)==1 .or. mod(tj2,2)==1)
     &     call errquit('WFAT: The argument tj1 and tj2 of '//
     &     'wfat_test_3j000 must be even integers.', 39, UERR)
      
      write(luout, '("Testing 3-j 000 symbols:")')
      j1 = tj1 / 2.0d0
      j2 = tj2 / 2.0d0

      
      n = min(tj1,tj2) + 1
      if (.not. MA_alloc_get(MT_DBL, n, '3j', l_cf, k_cf))
     &     call errquit('WFAT: 3j allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, n, '3j_000', l_cf0, k_cf0))
     &     call errquit('WFAT: 3j_000 allocation failed', 59, MA_ERR)
      call wfat_3j_jscan(tj1, 0, tj2, 0, dbl_mb(k_cf))
      call wfat_3j_jscan000(idnint(j1), idnint(j2), dbl_mb(k_cf0))
      write(luout, '("j1 = ", f7.1)') j1
      write(luout, '("m1 = ", f7.1)') 0
      write(luout, '("j2 = ", f7.1)') j2
      write(luout, '("m2 = ", f7.1)') 0
      write(luout, '("m3 = ", f7.1)') 0
      norm = 0.0d0
      do i = 1, n
         j3 = abs(j1-j2) + i - 1
         write(luout, '(f7.1, 2f15.8)')
     &        j3, dbl_mb(k_cf+i-1), dbl_mb(k_cf0+i-1)
         norm = norm + (2*j3+1)*dbl_mb(k_cf0+i-1)**2
      enddo
      write(luout, '("norm^2 = ", f12.8)') norm
      write(luout, '("")')
      
      
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: 3j deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_cf0))
     &     call errquit('WFAT: 3j_000 deallocation failed', 61, MA_ERR)

      end
!==================================================
      
      
