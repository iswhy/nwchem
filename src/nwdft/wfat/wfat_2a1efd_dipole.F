!========================================================
!========================================================
      subroutine wfat_2a1efd_dpo(comp, m, p, nr, nt1, nt2, k1, k2,
     &     cf1, cf2, rq, dt_in, dt, numt, dpo)
cc     Calculates the dipole matrix element between two wave functions
cc     using their coefficients in the FEDVR basis, that is
cc          D_{ij} = < u_{m+p,k'} | comp | v_{m,k} >,
cc     where
cc          v_{m,k} = \sum_i c_i \chi_{i,m,k}
cc          \chi_{i01}(r,\theta,\phi) = f_i(r,\theta)/\sqrt(2*pi)
cc          \chi_{i02}(r,\theta,\phi) = 0
cc          \chi_{im1}(r,\theta,\phi) = f_i(r,\theta) cos(m\phi)/\sqrt(pi)
cc          \chi_{im2}(r,\theta,\phi) = f_i(r,\theta) sin(m\phi)/\sqrt(pi)
cc     where comp = x,y,z and the size of D is n1 by n2.
      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      character*1 comp
      integer m
      integer p
      integer nr
      integer nt1
      integer nt2
      integer k1
      integer k2
      double precision cf1(nr*nt1)    ! [input], m-1 or m+1
      double precision cf2(nr*nt2)    ! [input], m
      double precision rq(nr)
      logical dt_in
      double precision dt(1:*)         ! [input], required only if dt_in=.true., regardless of dt_in, its minimum size should be nt1*nt2.
      logical numt                     ! [input], meaningful only if dt_in=.false.
      double precision dpo             ! [output],

      integer l_tm, k_tm
      integer i, ic1, ic2, delkron
      double precision s0, cc
      external delkron


      if (comp == 'x' .and. (abs(p) /= 1 .or. k1 /= k2)) then
         dpo = 0.0d0
         return
      endif
      if (comp == 'y' .and. (abs(p) /= 1 .or. k1 == k2)) then
         dpo = 0.0d0
         return
      endif
      if (comp == 'z' .and. (p /= 0 .or. k1 /= k2)) then
         dpo = 0.0d0
         return
      endif


      if (m < 0)
     &     call errquit('WFAT: The argument m of wfat_2a1e_dpo '//
     &     'must be zero or a positive integer.', 85, UERR)
      if (m == 0 .and. p == -1)
     &     call errquit('WFAT: If the argument m of wfat_2a1e_dpo '//
     &     'has the value of zero, then the allowed values of p is 0 '//
     &     'or 1.', 62, UERR)
      if (k1 /= 1 .and. k1 /= 2)
     &     call errquit('WFAT: The argument k1 of wfat_2a1e_dpo '//
     &     'must be equal to either 1 or 2.', 59, UERR)
      if (k2 /= 1 .and. k2 /= 2)
     &     call errquit('WFAT: The argument k2 of wfat_2a1e_dpo '//
     &     'must be equal to either 1 or 2.', 59, UERR)
      if (comp /= 'x' .and. comp /= 'y' .and. comp /= 'z')
     &     call errquit('WFAT: The argument comp of wfat_2a1e_dpo '//
     &     'must be equal to either x, y, or z.', 59, UERR)
      
      if (.not. dt_in) then
         if (numt) then
            call wfat_dipole_dvrt_num(m, p, nt1, nt2, dt)
         else
            call wfat_dipole_dvrt(m, p, nt1, nt2, dt)
         endif
      endif
      

      if (.not. ma_alloc_get(MT_DBL, nt1, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
      dpo = 0.0d0
      do i = 1, nr
         ic1 = 1 + (i-1)*nt1
         ic2 = 1 + (i-1)*nt2

         call wfat_dgemm_wrap('n', 'n', nt1, nt2, nt2, 1, 1.0d0, 
     &        dt, cf2(ic2), 0.0d0, dbl_mb(k_tm))
         s0 = sum(cf1(ic1:ic1+nt1-1) * dbl_mb(k_tm:k_tm+nt1-1))

         dpo = dpo + s0 * rq(i)
      enddo

      call wfat_2a1e_dpo_phi(comp, k1, m+p, k2, m, cc)
      dpo = cc * dpo
      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)      
      
      end
!========================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1efd_dpo_col(c, m1, m2, s1, s2, nr, nt1, 
     &     nt2, rq, dt, g_cf1, cf2, g_dp)

      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      character*1 c                      ! [input]
      integer m1                         ! [input]
      integer m2                         ! [input]
      integer s1                         ! [input]
      integer s2                         ! [input]
      integer nr                         ! [input]
      integer nt1                        ! [input]
      integer nt2                        ! [input]
      double precision rq(nr)            ! [input]
      double precision dt(nt1, nt2)      ! [input]
      integer g_cf1                      ! [input], GA of size nr*nt1 x nr*nt1. 1)
      double precision cf2(nr*nt2)       ! [input]
      integer g_dp                       ! [output], 1D GA of size nr*nt1.
c     1) It is recommended that the blocks of g_cf1 consist of columns (i.e. the
c        number of rows of each block is equal to the number of rows of g_dp), and
c        that the distribution of columns of g_cf1 is identical to that of g_dp.

      logical iscol
      integer p, i, l_cf1, k_cf1, itype, dim1, dim2, ilo, ihi
      integer jlo, jhi, klo, khi
      double precision dp

      call ga_sync()
      call ga_inquire(g_dp, itype, dim1, dim2)
      if (dim1 /= 1 .and. dim2 /= 1)
     &     call errquit('WFAT: Input argument g_dp of '//
     &     'wfat_2a1efd_dpo_col must one be dimensional.', 29, GA_ERR)
      
      if (c /= 'x' .and. c /= 'y' .and. c /= 'z')
     &     call errquit('WFAT: The argument c of wfat_2a1efd_dpo_col '//
     &     'must be equal to either x, y, or z.', 59, UERR)
      
      p = m1 - m2
      if (p < -1 .or. p > +1) then
         dp = 0.0d0
         return
      endif

      if (.not. ma_push_get(mt_dbl, nr*nt1, 'cf1', l_cf1, k_cf1))
     &     call errquit('WFAT: Could not allocate cf1.', 56, MA_ERR)

      call ga_distribution(g_dp, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (jlo == jhi) then
         klo = ilo
         khi = ihi
         iscol = .true.
      else
         klo = jlo
         khi = jhi
         iscol = .false.
      endif
cOLD      do i = 1, nr*nt1
      do i = klo, khi
         call ga_get(g_cf1, 1, nr*nt1, i, i, dbl_mb(k_cf1), nr*nt1)
         call wfat_2a1efd_dpo(c, m2, p, nr, nt1, nt2, s1, s2,
     &        dbl_mb(k_cf1), cf2, rq, .true., dt, .false., dp)
         if (iscol) then
            call ga_put(g_dp, i, i, 1, 1, dp, 1)
         else
            call ga_put(g_dp, 1, 1, i, i, dp, 1)
         endif
      enddo

      if (.not. ma_pop_stack(l_cf1))
     &     call errquit('WFAT: Could not deallocate cf1', 11, MA_ERR)
      call ga_sync()
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1efd_polb(mch, sphi, nr, nt, nt_m, nt_p, rq, cf,
     &     cfx, cfy, cfz, polb)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer mch
      integer sphi
      integer nr
      integer nt
      integer nt_m
      integer nt_p
      double precision rq(nr)
      double precision cf(nr*nt)
      double precision cfx(1:*)
      double precision cfy(1:*)
      double precision cfz(nr*nt)
      double precision polb(3)

      integer l_dtm, k_dtm, l_dtp, k_dtp, l_dt, k_dt
      integer n_m, n_p, sphi_
      double precision pb0
      

c     ==== m - 1 ====
      if (mch >= 1) then
         n_m = nr * nt_m
         if (.not. ma_alloc_get(MT_DBL, nt_m*nt, 'dt_dvr-', l_dtm,
     &        k_dtm))
     &        call errquit('WFAT: Cannot allocate dt_dvr-', 1, MA_ERR)
         call wfat_dipole_dvrt(mch-1, 1, nt, nt_m, dbl_mb(k_dtm))
      else
         n_m = 0
      endif

      
c     ==== m + 1 ====
      n_p = nr * nt_p
      if (.not. ma_alloc_get(MT_DBL, nt_p*nt, 'dt_dvr+', l_dtp, k_dtp))
     &     call errquit('WFAT: Cannot allocate dt_dvr+', 1, MA_ERR)
      call wfat_dipole_dvrt(mch+1, -1, nt, nt_p, dbl_mb(k_dtp))


c     ==== m,m ====
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'dt_dvr', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt_dvr', 1, MA_ERR)
      call wfat_dipole_dvrt(mch, 0, nt, nt, dbl_mb(k_dt))
      

c     ==== xx ====
      if (mch >= 1) then
         call wfat_2a1efd_dpo('x', mch-1, 1, nr, nt, nt_m, sphi, sphi,
     &        cf, cfx(1:n_m), rq, .true., dbl_mb(k_dtm), .false., pb0)
      else
         pb0 = 0.0d0
      endif
      call wfat_2a1efd_dpo('x', mch+1, -1, nr, nt, nt_p, sphi, sphi,
     &     cf, cfx(n_m+1:n_m+n_p), rq, .true., dbl_mb(k_dtp), .false.,
     &     polb(1))
      polb(1) = pb0 + polb(1)


c     ==== yy ====
      if (sphi==1) sphi_ = 2
      if (sphi==2) sphi_ = 1
      if (mch >= 1) then
         call wfat_2a1efd_dpo('y', mch-1, 1, nr, nt, nt_m, sphi, sphi_,
     &        cf, cfy(1:n_m), rq, .true., dbl_mb(k_dtm), .false., pb0)
      else
         pb0 = 0.0d0
      endif
      call wfat_2a1efd_dpo('y', mch+1, -1, nr, nt, nt_p, sphi, sphi_,
     &     cf, cfy(n_m+1:n_m+n_p), rq, .true., dbl_mb(k_dtp), .false.,
     &     polb(2))
      polb(2) = pb0 + polb(2)

      
c     ==== zz ====
      call wfat_2a1efd_dpo('z', mch, 0, nr, nt, nt, sphi, sphi,
     &     cf, cfz, rq, .true., dbl_mb(k_dt), .false., polb(3))
      polb(3) = polb(3)  


      polb = -2 * polb


c     ==== Deallocations ====
      if (mch >= 1) then
         if (.not. ma_free_heap(l_dtm))
     &        call errquit('WFAT: Cannot deallocate dt_dvr-',
     &        12, MA_ERR)         
      endif
      if (.not. ma_free_heap(l_dtp))
     &     call errquit('WFAT: Cannot deallocate dt_dvr+', 12, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt_dvr', 12, MA_ERR)
      
      end
c===============================================================
      

      
