c==========================================================
c==========================================================
      subroutine wfat_test_dvrr(m, nr, ntpz)

      implicit none
#include "errquit.fh"
      integer m
      integer nr
      integer ntpz

      double precision ke


      call errquit('WFAT: The wfat_test_dvrr subroutine has not '//
     &     'been revised following the debugging related to ke. First'//
     &     ' fix this before using this routine.', 32, UERR)
      
      call wfat_test_ovl_dvr(m, nr, ke, ntpz)
      call wfat_test_laplr_dvr(m, nr, ke, ntpz)

      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_test_ovl_dvr(m, nr, ke, ntpz)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      double precision ke
      integer ntpz

      integer i, j, k, l_tm, k_tm, l_s, k_s, i_s, i_tm0, i_tmn, i_tm
      integer j_tm, l_u, k_u, l_dv, k_dv, i_dv, l_x, k_x, l_rw, k_rw
      double precision x0, xn, dx, x, s, s1, sn, w


      call errquit('WFAT: The wfat_test_ovl_dvr subroutine has not '//
     &     'been revised following the debugging related to ke. First'//
     &     ' fix this before using this routine.', 32, UERR)

      ke = 1.0d0
      write(luout, '("Testing radial DVR overlap ...")')
      
      if (.not. ma_alloc_get(MT_DBL, ntpz*nr, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'ovl', l_s, k_s))
     &     call errquit('WFAT: Cannot allocate ovl', 1, MA_ERR)


c     ==== Check the overlap of the DVR bases ====
c     == Construct the trapezoid points ==
      x0 = 0.0d0
      xn = 50.0d0
      dx = (xn-x0) / (ntpz-1)
      do i = 1, ntpz
         x = x0 + (i-1)*dx
         i_tm0 = k_tm + (1-1)*ntpz + (i-1)      ! The first DVR at x.
         i_tmn = k_tm + (nr-1)*ntpz + (i-1)     ! The last DVR at x.
         call wfat_eval_dvrr(m, nr, ke, x/ke, dbl_mb(i_tm0:i_tmn:ntpz))   ! 1)
      enddo
c     Notes:
c     1) The row of tmp are the trapezoid points, the columns are the DVR functions.

      
      do j = 1, nr
         do i = 1, nr
c           ==== Trapezoid rule ====
            s = 0.0d0
            do k = 2, ntpz-1
               i_tm = k_tm+(i-1)*ntpz+(k-1)
               j_tm = k_tm+(j-1)*ntpz+(k-1)
               s = s + dbl_mb(i_tm) * dbl_mb(j_tm)
            enddo
            s1 = dbl_mb(k_tm+(i-1)*ntpz) * dbl_mb(k_tm+(j-1)*ntpz)
            sn = dbl_mb(k_tm+i*ntpz-1) * dbl_mb(k_tm+j*ntpz-1)

            i_s = k_s + (j-1)*nr + (i-1)
            dbl_mb(i_s) = dx * (0.5d0*s1 + s + 0.5d0*sn)
         enddo
      enddo

      write(luout, '("DVR overlap - trapezoid:")')
      call wfat_dprintmat(nr, nr, 'f', dbl_mb(k_s))


c     ==== Check the product transpose(u_r) * u_r ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'u_r', l_u, k_u))
     &     call errquit('WFAT: Cannot allocate u_r', 1, MA_ERR)
      call wfat_dvr_fbr_r(m, nr, dbl_mb(k_u))
      call wfat_dgemm_wrap('t', 'n', nr, nr, nr, nr, 1.0d0, dbl_mb(k_u),
     &     dbl_mb(k_u), 0.0d0, dbl_mb(k_s))
      write(luout, '("transpose(u_r) * u_r:")')
      call wfat_dprintmat(nr, nr, 'f', dbl_mb(k_s))


c     ==== Check the values of the DVR functions at the quadrature points ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'dvr_at_pt', l_dv, k_dv))
     &     call errquit('WFAT: Cannot allocate dvr_at_pt', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'xweight', l_rw, k_rw))
     &     call errquit('WFAT: Cannot allocate xweight', 1, MA_ERR)
      call wfat_lf_quadrature_rule(nr, abs(dble(m))+2, dbl_mb(k_x),
     &     dbl_mb(k_rw))
      write(luout, '("The values of the DVR functions (row) at the '//
     &     'quadrature points (column):")')
      write(luout, '(a5)', advance='no') 'Ref:'
      do i = 1, nr              ! Loop over radial quadrature points
         x = dbl_mb(k_x + i - 1)
         i_dv = k_dv + (i-1)*nr
         call wfat_eval_dvrr(abs(m), nr, ke, x/ke, dbl_mb(i_dv))
         w = dbl_mb(k_rw + i - 1)
         write(luout, '(f10.5)', advance='no')
     &        exp(-x/2) * x**(abs(dble(m))/2) / sqrt(w)
      enddo
      write(luout, '("")')
      call wfat_dprintmat(nr, nr, 'f', dbl_mb(k_dv))
         

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)
      if (.not. ma_free_heap(l_s))
     &     call errquit('WFAT: Cannot deallocate ovl', 1, MA_ERR)
      if (.not. ma_free_heap(l_u))
     &     call errquit('WFAT: Cannot deallocate u_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_dv))
     &     call errquit('WFAT: Cannot deallocate dvr_at_pt', 1, MA_ERR)
      if (.not. ma_free_heap(l_x))
     &     call errquit('WFAT: Cannot deallocate xquad', 1, MA_ERR)
      if (.not. ma_free_heap(l_rw))
     &     call errquit('WFAT: Cannot deallocate xweight', 1, MA_ERR)
      
      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_test_fbrr(m, nr, ntpz)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      integer ntpz

      integer i, j, k, l_tm, k_tm, l_s, k_s, i_s, i_tm0, i_tmn, i_tm
      integer j_tm
      double precision x0, xn, dx, x, s, s1, sn, ke


      call errquit('WFAT: The wfat_test_fbrr subroutine has not been '//
     &     'revised following the debugging related to ke. First fix '//
     &     'this before using this routine.', 32, UERR)
      
      write(luout, '("Testing radial FBR overlap ...")')
      ke = 1.0d0
      
      if (.not. ma_alloc_get(MT_DBL, ntpz*nr, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'ovl', l_s, k_s))
     &     call errquit('WFAT: Cannot allocate ovl', 1, MA_ERR)


c     ==== Check the overlap of the FBR bases ====
c     == Construct the trapezoid points ==
      x0 = 0.0d0
      xn = 50.0d0
      dx = (xn-x0) / (ntpz-1)
      do i = 1, ntpz
         x = x0 + (i-1)*dx
         i_tm0 = k_tm + (1-1)*ntpz + (i-1)      ! The first FBR at x.
         i_tmn = k_tm + (nr-1)*ntpz + (i-1)     ! The last FBR at x.
         call wfat_eval_fbrr(m, nr, ke, x/ke, dbl_mb(i_tm0:i_tmn:ntpz))
      enddo
      

      do j = 1, nr
         do i = 1, nr
c           ==== Trapezoid rule ====
            s = 0.0d0
            do k = 2, ntpz-1
               i_tm = k_tm+(i-1)*ntpz+(k-1)
               j_tm = k_tm+(j-1)*ntpz+(k-1)
               s = s + dbl_mb(i_tm) * dbl_mb(j_tm)
            enddo
            s1 = dbl_mb(k_tm+(i-1)*ntpz) * dbl_mb(k_tm+(j-1)*ntpz)
            sn = dbl_mb(k_tm+i*ntpz-1) * dbl_mb(k_tm+j*ntpz-1)

            i_s = k_s + (j-1)*nr + (i-1)
            dbl_mb(i_s) = dx * (0.5d0*s1 + s + 0.5d0*sn)
         enddo
      enddo

      write(luout, '("FBR overlap - trapezoid:")')
      call wfat_dprintmat(nr, nr, 'f', dbl_mb(k_s))


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)
      if (.not. ma_free_heap(l_s))
     &     call errquit('WFAT: Cannot deallocate ovl', 1, MA_ERR)
      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_test_laplr_dvr(m, nr, ke, ntpz)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      integer ntpz
      double precision ke

      integer l_dv, l_dv2, l_lr, l_fb, l_ur
      integer k_dv, k_dv2, k_lr, k_fb, k_ur
      integer i, j, k, i_dv0, i_dvn, i_dv, i_dv2, i_lr
      double precision x0, xn, dx, x, s, s1, sn


      call errquit('WFAT: The wfat_test_laplr_dvr subroutine has not '//
     &     'been revised following the debugging related to ke. First'//
     &     ' fix this before using this routine.', 32, UERR)
      
      write(luout, '("Testing r-Laplacian in DVR basis ...")')
      
      
c     ==== Calculate the DVR at the trapezoid points ==
      if (.not. ma_alloc_get(MT_DBL, ntpz*nr, 'dvr', l_dv, k_dv))
     &     call errquit('WFAT: Cannot allocate dvr', 1, MA_ERR)
      x0 = 0.0d0
      xn = 50.0d0
      dx = (xn-x0) / (ntpz-1)
      do i = 1, ntpz
         x = x0 + (i-1)*dx
         i_dv0 = k_dv + (1-1)*ntpz + (i-1)      ! The first DVR at x.
         i_dvn = k_dv + (nr-1)*ntpz + (i-1)     ! The last DVR at x.
         call wfat_eval_dvrr(m, nr, ke, x/ke, dbl_mb(i_dv0:i_dvn:ntpz))   ! 1)
         dbl_mb(i_dv0:i_dvn:ntpz) = x * dbl_mb(i_dv0:i_dvn:ntpz)   ! x \rho_i^{|m|}(x), for i = 1,2,...,nr.
      enddo
c     Notes:
c     1) The row of dvr are the trapezoid points, the columns are the DVR functions.


c     ==== Calculate the 2nd derivative of the DVR at the trapezoid points ====
      if (.not. ma_alloc_get(MT_DBL, ntpz*nr, 'd2dvr', l_dv2, k_dv2))
     &     call errquit('WFAT: Cannot allocate d2dvr', 1, MA_ERR)
      do i = 1, nr              ! Loop over DVR functions
         i_dv = k_dv + (i-1)*ntpz
         i_dv2 = k_dv2 + (i-1)*ntpz
         call wfat_num_derivative2(ntpz, dx, dbl_mb(i_dv),
     &        dbl_mb(i_dv2))
      enddo


c     ==== Construct the T_r^{DVR} matrix numerically ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'lapr_dv', l_lr, k_lr))
     &     call errquit('WFAT: Cannot allocate lapr_dv', 1, MA_ERR)
      do j = 1, nr         ! Loop over the 2nd derivative of x rho_j(x)
         do i = 1, nr         ! Loop over x rho_i(x)
c           ==== Trapezoid rule ====
            s = 0.0d0
            do k = 2, ntpz-1
               i_dv = k_dv + (i-1)*ntpz + (k-1)
               i_dv2 = k_dv2 + (j-1)*ntpz + (k-1)
               s = s + dbl_mb(i_dv) * dbl_mb(i_dv2)
            enddo
            s1 = dbl_mb(k_dv+(i-1)*ntpz) * dbl_mb(k_dv2+(j-1)*ntpz)
            sn = dbl_mb(k_dv+i*ntpz-1) * dbl_mb(k_dv2+j*ntpz-1)

            i_lr = k_lr + (j-1)*nr + (i-1)
            dbl_mb(i_lr) = dx * (0.5d0*s1 + s + 0.5d0*sn)
         enddo
      enddo

      write(luout, '("The r-Laplacian matrix in the DVR basis '//
     &     '(trapezoid rule): ")')
      call wfat_dprintmat(nr, nr, 'f', dbl_mb(k_lr))


c     ==== Construct the T_r^{DVR} matrix analytically ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'lapr_fb', l_fb, k_fb))
     &     call errquit('WFAT: Cannot allocate lapr_fb', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'u_r', l_ur, k_ur))
     &     call errquit('WFAT: Cannot allocate u_r', 1, MA_ERR)
      call wfat_laplr_fbrr(nr, abs(m), dbl_mb(k_fb))
      call wfat_dvr_fbr_r(abs(m), nr, dbl_mb(k_ur))
c     transpose(U_r) . T_r . U_r
      call wfat_represent_sym(nr, nr, dbl_mb(k_ur), dbl_mb(k_fb),
     &     dbl_mb(k_lr))
      write(luout, '("The r-Laplacian matrix in the DVR basis '//
     &     '(exact): ")')
      call wfat_dprintmat(nr, nr, 'f', dbl_mb(k_lr))
         

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dv))
     &     call errquit('WFAT: Cannot deallocate dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dv2))
     &     call errquit('WFAT: Cannot deallocate d2dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_lr))
     &     call errquit('WFAT: Cannot deallocate lapr_dv', 1, MA_ERR)
      if (.not. ma_free_heap(l_fb))
     &     call errquit('WFAT: Cannot deallocate lapr_fb', 1, MA_ERR)
      if (.not. ma_free_heap(l_ur))
     &     call errquit('WFAT: Cannot deallocate u_r', 1, MA_ERR)
      
      end
c==========================================================

      
c==========================================================
c==========================================================


c==========================================================
      
      
