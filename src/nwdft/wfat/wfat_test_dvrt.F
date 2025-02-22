c==========================================================
c==========================================================
      subroutine wfat_test_dvrt(m, nt, ntpz)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nt
      integer ntpz

      integer i
      integer l_s, l_u, l_pl
      integer k_s, k_u, k_pl


      write(luout, '("Testing polar DVR ...")')
      
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'ovl', l_s, k_s))
     &     call errquit('WFAT: Cannot allocate ovl', 1, MA_ERR)

      
c     ==== Check the product transpose(u_t) * u_t ====
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'u_t', l_u, k_u))
     &     call errquit('WFAT: Cannot allocate u_t', 1, MA_ERR)
      call wfat_dvr_fbr_t(m, nt, dbl_mb(k_u))
      write(luout, '("u_t:")')
      call wfat_dprintmat(nt, nt, 'f', dbl_mb(k_u))
      
      call wfat_dgemm_wrap('t', 'n', nt, nt, nt, nt, 1.0d0, dbl_mb(k_u),
     &     dbl_mb(k_u), 0.0d0, dbl_mb(k_s))
      write(luout, '("transpose(u_t) * u_t:")')
      call wfat_dprintmat(nt, nt, 'f', dbl_mb(k_s))


      if (.not. ma_alloc_get(MT_DBL, nt, 'P_lm', l_pl, k_pl))
     &     call errquit('WFAT: Cannot allocate P_lm', 1, MA_ERR)
      call wfat_assocleg_fixm(m, m+nt-1, 0.4d0, dbl_mb(k_pl))
      write(luout, '("Values of P_lm with m = ", i3, " :")') m
      do i = 1, nt
         write(luout, '("l = ", i4, f15.8)')
     &        abs(m)+i-1, dbl_mb(k_pl+i-1)
      enddo

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_s))
     &     call errquit('WFAT: Cannot deallocate ovl', 1, MA_ERR)
      if (.not. ma_free_heap(l_u))
     &     call errquit('WFAT: Cannot deallocate u_t', 1, MA_ERR)
      if (.not. ma_free_heap(l_pl))
     &     call errquit('WFAT: Cannot deallocate P_lm', 1, MA_ERR)
      
      end
c==========================================================
