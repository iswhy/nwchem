c===============================================================
c===============================================================
      subroutine wfat_2a1e_prep_system(rtdb, mch, sphi, spinid, orbid, 
     &     pottype, nel, potparam)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb                  ! [input],
      integer mch                   ! [output],
      integer sphi                  ! [output],
      integer spinid                ! [output],
      integer orbid                 ! [output],
      character*3 pottype           ! [output],
      integer nel                   ! [output], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)  ! [output],

      
c     ==== The m (magnetic) quantum number ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_m', mt_int, 1, mch))
     &     call errquit('WFAT: Failed reading 2a1e_m from rtdb',
     &     1, RTDB_ERR)

c     ==== The k quantum number (the xz-plane reflection symmetry) ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_xz', mt_int, 1, sphi))
     &     call errquit('WFAT: Failed reading 2a1e_xz from rtdb',
     &     1, RTDB_ERR)
      
c     ==== The eigenstate ID within the above (m,k) ====
      if (.not. rtdb_get(rtdb, 'wfat:orbid', mt_int, 1, orbid))
     &     call errquit('WFAT: Failed reading orbid from rtdb',
     &     18, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:spinid', mt_int, 1, spinid))
     &     call errquit('WFAT: Failed reading spinid from rtdb',
     &     18, INPUT_ERR)

c     ==== The Potential parameters ====
      if (.not. rtdb_cget(rtdb, 'wfat:2a1e_pot', 1, pottype))
     &     call errquit('WFAT: Failed reading 2a1e_pot from RTDB', 1,
     &     RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'wfat:2a1e_potparam', mt_dbl, 2,
     &     potparam))
     &     call errquit('WFAT: Failed reading 2a1e_potparam from RTDB',
     &     1, RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'wfat:2a1e_nel', mt_int, 1, nel))
     &     call errquit('WFAT: Failed reading 2a1e_nel from RTDB', 1, 
     &     RTDB_ERR)
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_prep_bas(rtdb, fedvr, init_ke, mult_ke, 
     &     maxiter, nr, nt, nfe, maxn, l_ndvr, k_ndvr, l_fen, k_fen, 
     &     l_xquad, k_xquad, l_wquad, k_wquad)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb              ! [input],
      logical fedvr             ! [input],
      double precision init_ke   ! [output],
      double precision mult_ke   ! [output],
      integer maxiter            ! [output],
      integer nr                 ! [output],
      integer nt                ! [output],
      integer nfe               ! [output],
      integer maxn               ! [output],
      integer l_ndvr, k_ndvr     ! [output],
      integer l_fen, k_fen       ! [output],
      integer l_xquad, k_xquad   ! [output],
      integer l_wquad, k_wquad   ! [output],

      integer nr_temp(3), nt_temp(3)
      

      if (fedvr) then
c        ==== No. of FEs ====
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_nfe', mt_int, 1, nfe))
     &        call errquit('WFAT: Failed reading 2a1e_nfe from rtdb',
     &        1, RTDB_ERR)
         
         init_ke = 0.0d0
         maxiter = -1
      else
         nfe = 1         
         
c        ==== Initial basis exponent ===
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_inke', mt_dbl, 1, init_ke))
     &        call errquit('WFAT: Failed reading 2a1e_inke from rtdb',
     &        1, RTDB_ERR)

c        ==== Exponential multipler ===
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_multke', mt_dbl, 1,
     &        mult_ke))
     &        call errquit('WFAT: Failed reading 2a1e_multke from rtdb',
     &        1, RTDB_ERR)
         
c        ==== Maximum iteration ====
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_maxit', mt_int, 1,
     &        maxiter))
     &        call errquit('WFAT: Failed reading 2a1e_maxit from rtdb',
     &        1, RTDB_ERR)
      endif


c     ==== No. of DVRs in each FE ====
      if (.not. ma_alloc_get(MT_INT, nfe, 'ndvr', l_ndvr, k_ndvr))
     &     call errquit('WFAT: Cannot allocate ndvr', 3, MA_ERR)
      if (fedvr) then
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_ndvr', mt_int, nfe,
     &        int_mb(k_ndvr))) call errquit('WFAT: Failed reading '//
     &        '2a1e_ndvr from rtdb', 1, RTDB_ERR)
      endif

         
c     ==== The boundary nodes of the FEs ====
      if (.not. ma_alloc_get(MT_DBL, nfe, 'fenodes', l_fen, k_fen))
     &     call errquit('WFAT: Cannot allocate fenodes', 3, MA_ERR)
      if (fedvr) then
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_fenodes', mt_dbl, nfe, 
     &        dbl_mb(k_fen))) call errquit('WFAT: Failed reading '//
     &        '2a1e_fenodes from rtdb', 1, RTDB_ERR)
      endif


c     ==== FEDVR points and weights in all FEs ====
      if (fedvr) then
         maxn = maxval(int_mb(k_ndvr:k_ndvr+nfe-1))
      else
         maxn = 1
      endif
      if (.not. ma_alloc_get(MT_DBL, maxn*nfe, 'xquad',
     &     l_xquad, k_xquad))
     &     call errquit('WFAT: Cannot allocate xquad', 3, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, maxn*nfe, 'wquad',
     &     l_wquad, k_wquad))
     &     call errquit('WFAT: Cannot allocate wquad', 3, MA_ERR)
      if (fedvr)
     &     call wfat_fedvr_rule(nfe, int_mb(k_ndvr), maxn, 
     &     dbl_mb(k_fen), dbl_mb(k_xquad), dbl_mb(k_wquad))
      
         
c     ==== The number of radial bases ====
      if (fedvr) then
cOLD         nr = (sum(int_mb(k_ndvr:k_ndvr+nfe-1)) - nfe*2) + (nfe - 1)
         call wfat_get_fedvrnr(nfe, int_mb(k_ndvr), nr)
      else
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_nr', mt_int, 3, nr_temp))
     &        call errquit('WFAT: Failed reading 2a1e_nr from rtdb',
     &        1, RTDB_ERR)
         nr = nr_temp(1)
      endif
         

c     ==== The number of polar DVR basis ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_nt', mt_int, 3, nt_temp))
     &     call errquit('WFAT: Failed reading 2a1e_nt from rtdb',
     &     1, RTDB_ERR)
      nt = nt_temp(1)
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_erase_bas(l_ndvr, l_fen, l_xquad, l_wquad)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer l_ndvr, l_fen, l_xquad, l_wquad

      if (.not. MA_free_heap(l_ndvr))
     &     call errquit('WFAT: ndvr deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_fen))
     &     call errquit('WFAT: fenodes deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_xquad))
     &     call errquit('WFAT: xquad deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_wquad))
     &     call errquit('WFAT: wquad deallocation failed', 3, MA_ERR)

      end
c===============================================================
      

c===============================================================
c===============================================================
      subroutine wfat_2a1e_prep_cf(nr, nt, g_cf, l_cf, k_cf, l_eg, k_eg)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer nr           ! [input]
      integer nt           ! [input]
      integer g_cf         ! [output]
      integer l_cf, k_cf   ! [output]
      integer l_eg, k_eg   ! [output]

      if (.not. ga_create(mt_dbl, nr*nt, nr*nt, '2a1e_coef', nr*nt, 0,
     &     g_cf))
     &     call errquit('WFAT: Cannot create coef global array.',
     &     1, GA_ERR)

      if (.not. ma_alloc_get(MT_DBL, nr*nt, '2a1e_ion_coef', l_cf,
     &     k_cf))
     &     call errquit('WFAT: Cannot allocate 2a1e_ion_coef',
     &     1, MA_ERR)

      if (.not. ma_alloc_get(MT_DBL, nr*nt, '2a1e_erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate 2a1e_erg', 1, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_prep_cf1(wfat1, fedvr, rtdb, m, nr, nt, nr_m,
     &     nt_m, nr_p, nt_p, l_cfx, k_cfx, l_cfy, k_cfy, l_cfz, k_cfz)

      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      logical wfat1          ! [input].
      logical fedvr          ! [input].
      integer rtdb           ! [input].
      integer m              ! [input].
      integer nr, nt         ! [input].
      integer nr_m, nt_m     ! [output].
      integer nr_p, nt_p     ! [output].
      integer l_cfx, k_cfx   ! [output].
      integer l_cfy, k_cfy   ! [output].
      integer l_cfz, k_cfz   ! [output].

      integer n_m, n_p, nz, nr_temp(3), nt_temp(3)


      if (.not. rtdb_get(rtdb, 'wfat:2a1e_nt', mt_int, 3, nt_temp))
     &     call errquit('WFAT: Failed reading 2a1e_nt from rtdb.',
     &     1, RTDB_ERR)
      if (.not. fedvr) then
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_nr', mt_int, 3, nr_temp))
     &        call errquit('WFAT: Failed reading 2a1e_nr from rtdb.',
     &        1, RTDB_ERR)
      endif
      
      if (wfat1) then
c        ==== The number of radial DVR basis for m-1 spectrum ====
         if (m >= 1) then
            if (fedvr) then
               nr_m = nr
            else
               nr_m = nr_temp(2)
            endif
            
c           ==== The number of polar DVR basis for m-1 spectrum ====
            nt_m = nt_temp(2)
         else
            nr_m = 0
            nt_m = 0
         endif
         
c        ==== The number of radial DVR basis for m+1 spectrum ====
         if (fedvr) then
            nr_p = nr
         else
            nr_p = nr_temp(3)
         endif
         
c        ==== The number of polar DVR basis for m+1 spectrum ====
         nt_p = nt_temp(3)
      else
         nr_m = 1
         nt_m = 1
         nr_p = 1
         nt_p = 1
      endif
      
      
c     ==== Allocations ====
      if (wfat1) then
         if (m >= 1) then
            n_m = nr_m * nt_m
         else
            n_m = 0
         endif
         n_p = nr_p * nt_p
         nz = nr * nt
      else
         n_m = 1
         n_p = 1
         nz = 1
      endif

      if (.not. ma_alloc_get(MT_DBL, n_m+n_p, '2a1e_coefx',
     &     l_cfx, k_cfx))
     &     call errquit('WFAT: Cannot allocate 2a1e_coefx', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n_m+n_p, '2a1e_coefy',
     &     l_cfy, k_cfy))
     &     call errquit('WFAT: Cannot allocate 2a1e_coefy', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nz, '2a1e_coefz',
     &     l_cfz, k_cfz))
     &     call errquit('WFAT: Cannot allocate 2a1e_coefz', 1, MA_ERR)

      end
c===============================================================
      
      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_erase_cf(g_cf, l_cf, l_eg)

      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      integer g_cf, l_cf, l_eg

      if (.not. ga_destroy(g_cf))
     &     call errquit('WFAT: Cannot destroy 2a1e_coef.', 45, GA_ERR)
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: 2a1e_ion_coef deallocation failed', 3,
     &     MA_ERR)
      if (.not. MA_free_heap(l_eg))
     &     call errquit('WFAT: 2a1e_erg deallocation failed', 3, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_erase_cf1(l_cfx, l_cfy, l_cfz)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer l_cfx, l_cfy, l_cfz

      if (.not. MA_free_heap(l_cfx))
     &     call errquit('WFAT: 2a1e_coefx deallocation failed', 3,
     &     MA_ERR)
      if (.not. MA_free_heap(l_cfy))
     &     call errquit('WFAT: 2a1e_coefy deallocation failed', 3,
     &     MA_ERR)
      if (.not. MA_free_heap(l_cfz))
     &     call errquit('WFAT: 2a1e_coefz deallocation failed', 3,
     &     MA_ERR)

      end
c===============================================================
      
      
