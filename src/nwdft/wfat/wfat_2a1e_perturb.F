c===============================================================
c===============================================================
      subroutine wfat_2a1e_perturb(geom, m, sphi, orbid, pottype, nel, 
     &     potparam, ke, mult, nr, nt, nr_m, nt_m, nr_p, nt_p, g_cf, 
     &     erg, cfx, cfy, cfz, polb)

      implicit none
#include "global.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                         ! [input],
      integer m                            ! [input],
      integer sphi                         ! [input],
      integer orbid                        ! [input],
      character*3 pottype
      integer nel                    ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)
      double precision ke                  ! [input],
      double precision mult                ! [input],
      integer nr                           ! [input],
      integer nt                           ! [input],
      integer nr_m                         ! [input],
      integer nt_m                         ! [input],
      integer nr_p                         ! [input],
      integer nt_p                         ! [input],
      integer g_cf                         ! [input], GA of size nr*nt x nr*nt
      double precision erg(nr*nt)          ! [input],
      double precision cfx(1:*)            ! [output], The minimum size is nr_p*nt_p if m=0, if m>0, the minimum size is nr_m*nt_m + nr_p*nt_p.
      double precision cfy(1:*)            ! [output], The minimum size is nr_p*nt_p if m=0, if m>0, the minimum size is nr_m*nt_m + nr_p*nt_p.
      double precision cfz(nr*nt)          ! [output],
      double precision polb(3)             ! [output],

      character*5 sm
      integer g_cfm, g_cfp, l_egm, k_egm, l_egp, k_egp, l_cf, k_cf
      integer n_m, n_p, ndid

      
      ndid = ga_nodeid()
      
c     ==== Solve the m-1 eigenvectors ====
      if (m >= 1) then
         n_m = nr_m * nt_m
         if (.not. ga_create(mt_dbl, n_m, n_m, 'cf-', n_m, 0, g_cfm))
     &        call errquit('WFAT: Cannot create cf- global array.',
     &        7, GA_ERR)
         if (.not. ma_alloc_get(MT_DBL, n_m, 'erg-', l_egm, k_egm))
     &        call errquit('WFAT: Cannot allocate erg-', 1, MA_ERR)
         call wfat_2a1e_solve(m-1, nr_m, nt_m, geom, ke, mult, .false., 
     &        1, 1, pottype, nel, potparam, g_cfm, dbl_mb(k_egm))

         if (ndid == 0) then
            write(sm, '(i5)') m-1
            write(luout, '(2/)')
            call util_print_centered(luout,
     &           'Eigenvalue spectrum of m = '//trim(adjustl(sm))//
     &           ' manifold (needed for WFAT(1))', 60, .true.)
            call wfat_2a1e_print_eigvals(n_m, 20, dbl_mb(k_egm))
         endif
      endif

c     ==== Solve the m+1 eigenvectors ====
      n_p = nr_p * nt_p
      if (.not. ga_create(mt_dbl, n_p, n_p, 'cf+', n_p, 0, g_cfp))
     &     call errquit('WFAT: Cannot create cf+ global array.',
     &     7, GA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n_p, 'erg+', l_egp, k_egp))
     &     call errquit('WFAT: Cannot allocate erg+', 1, MA_ERR)
      call wfat_2a1e_solve(m+1, nr_p, nt_p, geom, ke, mult, .false., 1,
     &     1, pottype, nel, potparam, g_cfp, dbl_mb(k_egp))

      if (ndid == 0) then
         write(sm, '(i5)') m+1
         write(luout, '(2/)')
         call util_print_centered(luout,
     &        'Eigenvalue spectrum of m = '//trim(adjustl(sm))//
     &        ' manifold (needed for WFAT(1))', 60, .true.)
         call wfat_2a1e_print_eigvals(n_p, 20, dbl_mb(k_egp))
      endif
      
c     ==== Obtain the coefficients of the x, y, z ====
c     ====    corrections to the wavefunction     ====
      call wfat_2a1e_perturbcf(geom, m, sphi, orbid, ke, nr,
     &     nt, nr_m, nt_m, nr_p, nt_p, g_cf, g_cfm, g_cfp,
     &     erg, dbl_mb(k_egm), dbl_mb(k_egp), cfx, cfy, cfz)
      
c     ==== Obtain the polarizability ====
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate coef.', 67, MA_ERR)
      call ga_get(g_cf, 1, nr*nt, orbid, orbid, dbl_mb(k_cf), nr*nt)
      call wfat_2a1e_polb(m, sphi, nr, nt, nr_m, nt_m, nr_p, nt_p, ke,
     &     dbl_mb(k_cf), cfx, cfy, cfz, .false., polb)
      write(luout, '("Polarizability : ", 3f15.6)') polb
      call errquit('The subroutine wfat_2a1e_polb above gives a '//
     &     'diverging result, this bug has not yet been resolved.',
     &     7, UERR)

c     ==== Deallocations ====
      if (m >= 1) then
         if (.not. ga_destroy(g_cfm))
     &        call errquit('WFAT: Cannot destroy cf-.', 22, GA_ERR)
         if (.not. ma_free_heap(l_egm))
     &        call errquit('WFAT: Cannot deallocate erg-', 1, MA_ERR)
      endif
      if (.not. ga_destroy(g_cfp))
     &     call errquit('WFAT: Cannot destroy cf+.', 22, GA_ERR)
      if (.not. ma_free_heap(l_egp))
     &     call errquit('WFAT: Cannot deallocate erg+', 1, MA_ERR)
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Cannot deallocate coef', 41, MA_ERR)

      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_perturbcf(geom, m, sphi, orbid, ke, nr,
     &     nt, nr_m, nt_m, nr_p, nt_p, g_cf, g_cfm, g_cfp, erg, erg_m,
     &     erg_p, cfx, cfy, cfz)

      implicit none
#include "global.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer m
      integer sphi
      integer orbid
      double precision ke
      integer nr
      integer nt
      integer nr_m
      integer nt_m
      integer nr_p
      integer nt_p
      integer g_cf             ! [input], GA of size nr*nt x nr*nt
      integer g_cfm            ! [input], GA of size nr_m*nt_m x nr_m*nt_m
      integer g_cfp            ! [input], GA of size nr_p*nt_p x nr_p*nt_p
      double precision erg(nr*nt)
      double precision erg_m(nr_m*nt_m)
      double precision erg_p(nr_p*nt_p)
      double precision cfx(1:*)
      double precision cfy(1:*)
      double precision cfz(nr*nt)

      integer j0, n, n_m, n_p, i, j, i_, sphi_, nblock(2), g_dp0, g_eg0
      integer g_t0, g_dpm, g_egm, g_tm, g_dpp, g_egp, g_tp, ndid, nnode
      integer l_drm, l_dtm, l_drp, l_dtp, l_dr, l_dt, l_ilo, l_jlo, l_cf
      integer k_drm, k_dtm, k_drp, k_dtp, k_dr, k_dt, k_ilo, k_jlo, k_cf


      ndid = ga_nodeid()

      j0 = orbid
      n = nr * nt
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate coef.', 14, MA_ERR)
      call ga_get(g_cf, 1, nr*nt, j0, j0, dbl_mb(k_cf), nr*nt)
      nnode = ga_nnodes()
      if (.not. ma_alloc_get(MT_INT, nnode, 'ilo', l_ilo, k_ilo))
     &     call errquit('WFAT: Cannot allocate ilo', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_INT, nnode, 'jlo', l_jlo, k_jlo))
     &     call errquit('WFAT: Cannot allocate jlo', 1, MA_ERR)
      
      
c     ==== m - 1 ====
      if (m >= 1) then
         n_m = nr_m * nt_m
         if (.not. ma_alloc_get(MT_DBL, nr_m*nr, 'dr_dvr-', l_drm, 
     &        k_drm))
     &        call errquit('WFAT: Cannot allocate dr_dvr-', 1, MA_ERR)
         call wfat_dipole_dvrr(m-1, nr_m, m, nr, ke, dbl_mb(k_drm))
   
         if (.not. ma_alloc_get(MT_DBL, nt_m*nt, 'dt_dvr-', l_dtm,
     &        k_dtm))
     &        call errquit('WFAT: Cannot allocate dt_dvr-', 1, MA_ERR)
         call wfat_dipole_dvrt(m, -1, nt_m, nt, dbl_mb(k_dtm))

c        ==== m-1 dipole vectors, 1/delta erg, and tmp ====
         call ga_nblock(g_cfm, nblock)
         call wfat_ga_get_lo(g_cfm, int_mb(k_ilo), int_mb(k_jlo))
         if (.not. ga_create_irreg(mt_dbl, n_m, 1, 'dp-', int_mb(k_jlo),
     &        nblock(2), 1, 1, g_dpm))
     &        call errquit('WFAT: Cannot create dp-', 1, GA_ERR)
         if (.not. ga_duplicate(g_dpm, g_egm, 'eg-'))
     &        call errquit('WFAT: Cannot create eg-', 88, GA_ERR)
         if (ndid == 0)
     &        call ga_put(g_egm, 1, n_m, 1, 1, 1/(erg(j0)-erg_m), n_m)
         if (.not. ga_duplicate(g_dpm, g_tm, 'tmp-'))
     &        call errquit('WFAT: Cannot create tmp-', 88, GA_ERR)
      else
         n_m = 0
      endif

      
c     ==== m + 1 ====
      n_p = nr_p * nt_p
      if (.not. ma_alloc_get(MT_DBL, nr_p*nr, 'dr_dvr+', l_drp, k_drp))
     &     call errquit('WFAT: Cannot allocate dr_dvr+', 1, MA_ERR)
      call wfat_dipole_dvrr(m+1, nr_p, m, nr, ke, dbl_mb(k_drp))

      if (.not. ma_alloc_get(MT_DBL, nt_p*nt, 'dt_dvr+', l_dtp, k_dtp))
     &     call errquit('WFAT: Cannot allocate dt_dvr+', 1, MA_ERR)
      call wfat_dipole_dvrt(m, 1, nt_p, nt, dbl_mb(k_dtp))

c     ==== m+1 dipole vectors, 1/delta erg, and tmp ====
      call ga_nblock(g_cfp, nblock)
      call wfat_ga_get_lo(g_cfp, int_mb(k_ilo), int_mb(k_jlo))
      if (.not. ga_create_irreg(mt_dbl, n_p, 1, 'dp+', int_mb(k_jlo), 
     &     nblock(2), 1, 1, g_dpp))
     &     call errquit('WFAT: Cannot create dp+', 1, GA_ERR)
      if (.not. ga_duplicate(g_dpp, g_egp, 'eg+'))
     &     call errquit('WFAT: Cannot create eg+', 88, GA_ERR)
      if (ndid == 0)
     &     call ga_put(g_egp, 1, n_p, 1, 1, 1/(erg(j0)-erg_p), n_p)
      if (.not. ga_duplicate(g_dpp, g_tp, 'tmp+'))
     &     call errquit('WFAT: Cannot create tmp+', 88, GA_ERR)
      

c     ==== x component of the corrected wave function ====
c     == m-1 sub-component ==
      if (m >= 1) then
         call wfat_2a1e_dpo_col('x', ke, m-1, m, sphi, sphi, nr_m, nt_m,
     &        nr, nt, dbl_mb(k_drm), dbl_mb(k_dtm), g_cfm, dbl_mb(k_cf),
     &        g_dpm)
         call ga_elem_multiply(g_dpm, g_egm, g_dpm)
         call ga_dgemm('n', 'n', n_m, 1, n_m, 1.0d0, g_cfm, g_dpm,
     &        0.0d0, g_tm)
         call ga_get(g_tm, 1, n_m, 1, 1, cfx(1), n_m)
         call ga_sync()
      endif
c     == m+1 sub-component ==
      call wfat_2a1e_dpo_col('x', ke, m+1, m, sphi, sphi, nr_p, nt_p,
     &     nr, nt, dbl_mb(k_drp), dbl_mb(k_dtp), g_cfp, dbl_mb(k_cf),
     &     g_dpp)
      call ga_elem_multiply(g_dpp, g_egp, g_dpp)
      call ga_dgemm('n', 'n', n_p, 1, n_p, 1.0d0, g_cfp, g_dpp,
     &     0.0d0, g_tp)
      call ga_get(g_tp, 1, n_p, 1, 1, cfx(n_m+1), n_p)
      call ga_sync()


c     ==== y component of the corrected wave function ====
      if (sphi==1) sphi_ = 2
      if (sphi==2) sphi_ = 1
c     == m-1 sub-component ==
      if (m >= 1) then
         call wfat_2a1e_dpo_col('y', ke, m-1, m, sphi_, sphi, nr_m, 
     &        nt_m, nr, nt, dbl_mb(k_drm), dbl_mb(k_dtm), g_cfm, 
     &        dbl_mb(k_cf), g_dpm)
         call ga_elem_multiply(g_dpm, g_egm, g_dpm)
         call ga_dgemm('n', 'n', n_m, 1, n_m, 1.0d0, g_cfm, g_dpm,
     &        0.0d0, g_tm)
         call ga_get(g_tm, 1, n_m, 1, 1, cfy(1), n_m)
         call ga_sync()
      endif
c     == m+1 sub-component ==
      call wfat_2a1e_dpo_col('y', ke, m+1, m, sphi_, sphi, nr_p, nt_p, 
     &     nr, nt, dbl_mb(k_drp), dbl_mb(k_dtp), g_cfp, dbl_mb(k_cf),
     &     g_dpp)
      call ga_elem_multiply(g_dpp, g_egp, g_dpp)
      call ga_dgemm('n', 'n', n_p, 1, n_p, 1.0d0, g_cfp, g_dpp,
     &     0.0d0, g_tp)
      call ga_get(g_tp, 1, n_p, 1, 1, cfy(n_m+1), n_p)
      call ga_sync()


c     ==== m,m ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'dr_dvr', l_dr, k_dr))
     &     call errquit('WFAT: Cannot allocate dr_dvr', 1, MA_ERR)
      call wfat_dipole_dvrr(m, nr, m, nr, ke, dbl_mb(k_dr))

      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'dt_dvr', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt_dvr', 1, MA_ERR)
      call wfat_dipole_dvrt(m, 0, nt, nt, dbl_mb(k_dt))
      
c     ==== m dipole vectors, 1/delta erg, and tmp ====
      call ga_nblock(g_cf, nblock)
      call wfat_ga_get_lo(g_cf, int_mb(k_ilo), int_mb(k_jlo))
      if (.not. ga_create_irreg(mt_dbl, n, 1, 'dp0', int_mb(k_jlo),
     &     nblock(2), 1, 1, g_dp0))
     &     call errquit('WFAT: Cannot create dp0', 1, GA_ERR)
      if (.not. ga_duplicate(g_dp0, g_eg0, 'eg0'))
     &     call errquit('WFAT: Cannot create eg0', 88, GA_ERR)
      if (ndid == 0) then
         call ga_put(g_eg0, 1, n, 1, 1, erg(j0)-erg, n)
         call ga_put(g_eg0, j0, j0, 1, 1, 1.0d0, 1)
      endif
      call ga_recip(g_eg0)
      if (ndid == 0) call ga_put(g_eg0, j0, j0, 1, 1, 0.0d0, 1)
      if (.not. ga_duplicate(g_dp0, g_t0, 'tmp0'))
     &     call errquit('WFAT: Cannot create tmp0', 88, GA_ERR)

      
c     ==== z component of the corrected wave function ====
      call wfat_2a1e_dpo_col('z', ke, m, m, sphi, sphi, nr, nt, nr, nt,
     &     dbl_mb(k_dr), dbl_mb(k_dt), g_cf, dbl_mb(k_cf), g_dp0)
      call ga_elem_multiply(g_dp0, g_eg0, g_dp0)
      call ga_dgemm('n', 'n', n, 1, n, 1.0d0, g_cf, g_dp0,
     &     0.0d0, g_t0)
      call ga_get(g_t0, 1, n, 1, 1, cfz(1), n)
      call ga_sync()
      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_ilo))
     &     call errquit('WFAT: Cannot deallocate ilo', 1, MA_ERR)
      if (.not. ma_free_heap(l_jlo))
     &     call errquit('WFAT: Cannot deallocate jlo', 1, MA_ERR)

      if (m >= 1) then
         if (.not. ma_free_heap(l_drm))
     &        call errquit('WFAT: Cannot deallocate dr_dvr-', 1, MA_ERR)
         if (.not. ma_free_heap(l_dtm))
     &        call errquit('WFAT: Cannot deallocate dt_dvr-', 1, MA_ERR)
         if (.not. ga_destroy(g_dpm))
     &        call errquit('WFAT: Cannot destroy dp-.', 45, GA_ERR)
         if (.not. ga_destroy(g_egm))
     &        call errquit('WFAT: Cannot destroy eg-.', 45, GA_ERR)
         if (.not. ga_destroy(g_tm))
     &        call errquit('WFAT: Cannot destroy tmp-.', 45, GA_ERR)
      endif
      if (.not. ga_destroy(g_dpp))
     &     call errquit('WFAT: Cannot destroy dp+.', 45, GA_ERR)
      if (.not. ga_destroy(g_egp))
     &     call errquit('WFAT: Cannot destroy eg+.', 45, GA_ERR)
      if (.not. ga_destroy(g_tp))
     &     call errquit('WFAT: Cannot destroy tmp.', 45, GA_ERR)
      if (.not. ga_destroy(g_dp0))
     &     call errquit('WFAT: Cannot destroy dp0.', 45, GA_ERR)
      if (.not. ga_destroy(g_eg0))
     &     call errquit('WFAT: Cannot destroy eg0.', 45, GA_ERR)
      if (.not. ga_destroy(g_t0))
     &     call errquit('WFAT: Cannot destroy tmp0.', 45, GA_ERR)
      
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Cannot deallocate coef', 77, MA_ERR)
      if (.not. ma_free_heap(l_drp))
     &     call errquit('WFAT: Cannot deallocate dr_dvr+', 1, MA_ERR)
      if (.not. ma_free_heap(l_dtp))
     &     call errquit('WFAT: Cannot deallocate dt_dvr+', 1, MA_ERR)
      if (.not. ma_free_heap(l_dr))
     &     call errquit('WFAT: Cannot deallocate dr_dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)
      
      end
c===============================================================


      
      
