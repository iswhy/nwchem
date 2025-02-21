c=======================================================
c=======================================================
      subroutine wfat_calc_q(nr, r, n, m, lmax, nlgl, erg, zq, q_1, q_2)

      implicit none
#include "global.fh"            
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "wfat_chgu.fh"
      integer nr
      double precision r(nr)    ! [input], increasingly ordered radial points
      integer n
      integer m
      integer lmax
      integer nlgl
      double precision erg
      double precision zq
      double precision q_1(abs(m):lmax, nr)       ! [output], orientation-independent
      double precision q_2(abs(m):lmax, nr)       ! [output], orientation-dependent

      integer i, l, i_qp1, i_qp2, i_d1, i_d2, nl
      integer l_qp1, l_qp2, l_d1, l_d2, k_qp1, k_qp2, k_d1, k_d2
      double precision kappa, lmb, a0, b0, aa, g, d1, d2
      double precision chgm(abs(m):lmax)
      

      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a0 = 1.0d0 - zq/kappa
      b0 = 2.0d0
      nl = lmax - abs(m) + 1
      call wfat_reset_lfqlib()
      do l = abs(m), MAXQLIB
         aa = l - zq/kappa
         if (aa >= 0.0d0) call wfat_update_lfqlib(aa)
      enddo
      
      
c     ==== Calculate the particular part of Q_l ====
      if (.not. MA_alloc_get(MT_DBL, nl*nr, 'q_p1', l_qp1, k_qp1))
     &     call errquit('WFAT: q_p1 allocation failed', 29, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl*nr, 'q_p2', l_qp2, k_qp2))
     &     call errquit('WFAT: q_p2 allocation failed', 29, MA_ERR)
      call wfat_q_particular(nr, r, n, m, lmax, erg, zq, nlgl,
     &     dbl_mb(k_qp1), dbl_mb(k_qp2))


c     ==== Calculate the coefficient of the homogeneous part of Q_l ====
      if (.not. MA_alloc_get(MT_DBL, nl, 'd1', l_d1, k_d1))
     &     call errquit('WFAT: d1 allocation failed', 29, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl, 'd2', l_d2, k_d2))
     &     call errquit('WFAT: d2 allocation failed', 29, MA_ERR)
      do l = abs(m), lmax
         i_d1 = k_d1 + l - abs(m)
         i_d2 = k_d2 + l - abs(m)
         call wfat_calc_d(n, m, l, erg, zq, dbl_mb(i_d1), dbl_mb(i_d2))
      enddo
      

c     ==== Calculate Q_l ====
      do i = 1, nr
         call wfat_chgm_1pt(r(i), abs(m), lmax, a0, b0, lmb, chgm)
         do l = abs(m), lmax
            g = (kappa*r(i))**l * exp(-kappa*r(i)) * chgm(l)

            d1 = dbl_mb(k_d1 + l - abs(m))
            d2 = dbl_mb(k_d2 + l - abs(m))
            i_qp1 = k_qp1 + (i-1)*nl + l - abs(m)
            q_1(l, i) = d1 * g + dbl_mb(i_qp1)
            i_qp2 = k_qp2 + (i-1)*nl + l - abs(m)
            q_2(l, i) = d2 * g + dbl_mb(i_qp2)

cOLD            if ( isnan(q_1(l, i)) ) q_1(l, i) = 0.0d0
cOLD            if ( isnan(q_2(l, i)) ) q_2(l, i) = 0.0d0
         enddo
      enddo

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_qp1))
     &     call errquit('WFAT: Cannot deallocate q_p1', 67, MA_ERR)
      if (.not. ma_free_heap(l_qp2))
     &     call errquit('WFAT: Cannot deallocate q_p2', 67, MA_ERR)
      if (.not. ma_free_heap(l_d1))
     &     call errquit('WFAT: Cannot deallocate d1', 67, MA_ERR)
      if (.not. ma_free_heap(l_d2))
     &     call errquit('WFAT: Cannot deallocate d2', 67, MA_ERR)

      end
c=======================================================
