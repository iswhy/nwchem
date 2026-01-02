c=================================================
c=================================================
      subroutine wfat_prep_tasktitle(ti)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      character*7 ti   ! [input]

      integer ndid

      
      ndid = ga_nodeid()
      if (trim(adjustl(ti)) /= 'OE' .and. trim(adjustl(ti)) /= 'ME'
     &     .and. ti /= 'OE-2a1e')
     &     call errquit('wfat_prep_tasktitle: Invalid value for the '//
     &     'argument ti. A valid value is ''OE'', ''ME'', or 
     &     ''OE-2a1e''.', 1, UERR)


      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(
     &        luout, '--------------------------------', 40,
     &        .false.)
         call util_print_centered(
     &        luout, '<<<<<< Entering '//trim(adjustl(ti))//
     &        '-WFAT task >>>>>>', 40, .false.)
         call util_print_centered(
     &        luout, '--------------------------------', 40,
     &        .false.)
      endif
      
      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_close_task(ti)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      character*2 ti

      integer ndid

      
      ndid = ga_nodeid()
      if (ti /= 'OE' .and. ti /= 'ME' .and. ti /= 'OE-2a1e')
     &     call errquit('wfat_close_task: Invalid value for the '//
     &     'argument ti. A valid value is ''OE'', ''ME'', or 
     &     ''OE-2a1e''.', 1, UERR)
      

      if (ndid == 0) then
         write(luout, '(/, 3x, "*** INFO ***")')
         write(luout, '(5x, a2, "WFAT task successfully carried out")')
     &        ti
         write(luout, '(/)')
         call util_print_centered(
     &        luout, '-------------------------------', 40,
     &        .false.)
         call util_print_centered(
     &        luout, '<<<<<< Exiting '//trim(adjustl(ti))//
     &        'WFAT task >>>>>>', 40, .false.)
         call util_print_centered(
     &        luout, '-------------------------------', 40,
     &        .false.)
      endif

      end
c=================================================
      
      
c=================================================
c=================================================
      subroutine wfat_prep_general(rtdb, gq, lmax_)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb   ! [input]
      integer gq     ! [output]
      integer lmax_  ! [output]


      if (.not. rtdb_get(rtdb, 'wfat:gqual', mt_int, 1, gq))
     &     call errquit('WFAT: Failed reading gqual from rtdb',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:lmax', mt_int, 1, lmax_))
     &     call errquit('WFAT: Failed reading lmax from rtdb',
     &     1, RTDB_ERR)
      
      end
c=================================================

      
c=================================================
c=================================================
      subroutine wfat_prep_channels(rtdb, n0, nch, l_or0, k_or0,
     &     l_sp, k_sp, l_or, k_or)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "global.fh"
      integer rtdb              ! [input],
      integer n0(2)          ! [input],
      integer nch            ! [output],
      integer l_or0, k_or0   ! [output],
      integer l_sp, k_sp     ! [output],
      integer l_or, k_or     ! [output],

      integer i, is


c     ==== Obtaining nch, the number of MO's ====
      if (.not. rtdb_get(rtdb, 'wfat:n_orbid', mt_int, 1, nch))
     &     call errquit('WFAT: Failed reading n_orbid from rtdb',
     &     1, INPUT_ERR)

      
c     ==== Obtaining an array of MO ID's ====
      if (.not. ma_alloc_get(MT_INT, nch, 'orbid0', l_or0, k_or0))
     &     call errquit('WFAT: Failed allocating orbid0', 1, MA_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:orbid', mt_int, nch,
     &     int_mb(k_or0)))
     &     call errquit('WFAT: Failed reading orbid from rtdb',
     &     1, INPUT_ERR)


c     ==== Obtaining an array of spin ID's ====
      if (.not. ma_alloc_get(MT_INT, nch, 'spinid', l_sp, k_sp))
     &     call errquit('WFAT: Failed allocating spinid', 1, MA_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:spinid', mt_int, nch,
     &     int_mb(k_sp)))
     &     call errquit('WFAT: Failed reading spinid from rtdb',
     &     1, INPUT_ERR)
      

c     ==== Constructing the MO ID array starting from 1 ====
      if (.not. ma_alloc_get(MT_INT, nch, 'orbid', l_or, k_or))
     &     call errquit('WFAT: Failed allocating orbid', 1, MA_ERR)
      do i = 1, nch
         is = int_mb(k_sp+i-1)
         int_mb(k_or+i-1) = int_mb(k_or0+i-1) + n0(is)
      enddo      

      
c     ==== Print channels ====
      call wfat_info_banner(luout, .true.)
      if (ga_nodeid() == 0)
     &     call wfat_print_channels(nch, int_mb(k_sp), int_mb(k_or0))
      call ga_sync()

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_channels(l_or0, l_sp, l_or)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_or0, l_sp, l_or

      if (.not. MA_free_heap(l_or0))
     &     call errquit('WFAT: orbid0 deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_sp))
     &     call errquit('WFAT: spinid deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_or))
     &     call errquit('WFAT: orbid deallocation failed', 3, MA_ERR)

      end
c=================================================

      
c=================================================
c=================================================
      subroutine wfat_prep_paraqn(rtdb, nqn, l_npa, k_npa, l_mpa, k_mpa)

      implicit none
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb           ! [input],
      integer nqn            ! [output],
      integer l_npa, k_npa   ! [output],
      integer l_mpa, k_mpa   ! [output],

      
c     ==== Obtaining nqn, the number of parabolic quantum numbers ====
      if (.not. rtdb_get(rtdb, 'wfat:n_npara', mt_int, 1, nqn))
     &     call errquit('WFAT: Failed reading n_npara from rtdb',
     &     1, INPUT_ERR)

      
c     ==== Obtaining n parabolic quantum numbers from RTDB ====
      if (.not. ma_alloc_get(MT_INT, nqn, 'npara', l_npa, k_npa))
     &     call errquit('WFAT: Failed allocating npara', 1, MA_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:npara', mt_int, nqn,
     &     int_mb(k_npa)))
     &     call errquit('WFAT: Failed reading npara from rtdb',
     &     1, INPUT_ERR)


c     ==== Obtaining m parabolic quantum numbers from RTDB ====
      if (.not. ma_alloc_get(MT_INT, nqn, 'mpara', l_mpa, k_mpa))
     &     call errquit('WFAT: Failed allocating mpara', 1, MA_ERR)    
      if (.not. rtdb_get(rtdb, 'wfat:mpara', mt_int, nqn,
     &     int_mb(k_mpa)))
     &     call errquit('WFAT: Failed reading mpara from rtdb',
     &     1, INPUT_ERR)      


c     ==== Printing ====
      call wfat_info_banner(luout, .true.)
      if (ga_nodeid() == 0)
     &     call wfat_print_paraqn(nqn, int_mb(k_npa), int_mb(k_mpa))
      call ga_sync()
      
      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_paraqn(l_npa, l_mpa)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_npa, l_mpa

      if (.not. MA_free_heap(l_npa))
     &     call errquit('WFAT: npara deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_mpa))
     &     call errquit('WFAT: mpara deallocation failed', 3, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_wfat1(rtdb, nch, nqn, wfat1, nlgl, nw1x,
     &     l_w1x, k_w1x)

      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb                  ! [input],
      integer nch                   ! [input],
      integer nqn                   ! [input],
      logical wfat1                 ! [output],
      integer nlgl                  ! [output],
      integer nw1x, l_w1x, k_w1x    ! [output],

      if (.not. rtdb_get(rtdb, 'wfat:wfat1', mt_log, 1, wfat1))
     &     call errquit('WFAT: Cannot read wfat1 from RTDB', 27,
     &     RTDB_ERR)
      if (wfat1) then
         if (.not. rtdb_get(rtdb, 'wfat:nlgl', mt_int, 1, nlgl))
     &        call errquit('WFAT: Cannot read nlgl from RTDB', 27,
     &        RTDB_ERR)
         nw1x = nqn * nch
      else
         nlgl = 0
         nw1x = 1
      endif

      if (.not. MA_alloc_get(MT_LOG, nw1x, 'wfat1_exc', l_w1x, k_w1x))
     &     call errquit('WFAT: wfat1_exc allocation failed', 59, MA_ERR)
      if (wfat1) then
         if (.not. rtdb_get(rtdb, 'wfat:wfat1_exclude', mt_log, nw1x, 
     &        log_mb(k_w1x))) call errquit('WFAT: Failed reading '//
     &        'dthr from rtdb', 22, RTDB_ERR)
      endif
      
      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_wfat1(l_w1x)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_w1x

      if (.not. MA_free_heap(l_w1x))
     &     call errquit('WFAT: wfat_exc deallocation failed', 3, MA_ERR)

      end
c=================================================

      
c=================================================
c=================================================
      subroutine wfat_prep_cf(nch, nb, wfat1, ncf, nch1, l_cf, k_cf)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nch         ! [input]
      integer nb          ! [input]
      logical wfat1       ! [input]
      integer ncf         ! [output]
      integer nch1        ! [output]
      integer l_cf, k_cf  ! [output]


      if (wfat1) then
         nch1 = 3*nch
      else
         nch1 = 0
      endif
      ncf = nb * (nch + nch1)
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vch', l_cf, k_cf))
     &     call errquit('WFAT: vch allocation failed', 59, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_cf(l_cf)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_cf

      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: vch deallocation failed', 3, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_acoef(nat, nqn, nch, wfat1, n_acg, l_acg,
     &     k_acg, n_acj, l_acj, k_acj, n_ack, l_ack, k_ack, n_dum,
     &     l_dum, k_dum)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nat
      integer nqn
      integer nch                   ! [input]. The number of ionized channels. It cannot be extended to include their 1st order corrections.
      logical wfat1
      integer n_acg, l_acg, k_acg
      integer n_acj, l_acj, k_acj
      integer n_ack, l_ack, k_ack
      integer n_dum, l_dum, k_dum

      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      
      n_acg = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, n_acg, 'acg', l_acg, k_acg))
     &     call errquit('WFAT: acg allocation failed', 81, MA_ERR)
      if (wfat1) then
         n_acj = nat * nqn * 3 * nch
         if (.not. MA_alloc_get(MT_DCPL, n_acj, 'acj', l_acj, k_acj))
     &        call errquit('WFAT: acj allocation failed', 81, MA_ERR)
         n_ack = nat * nqn * nch
         if (.not. MA_alloc_get(MT_DCPL, n_ack, 'ack', l_ack, k_ack))
     &        call errquit('WFAT: ack allocation failed', 81, MA_ERR)
         n_dum = nat * nqn * nch
         if (.not. MA_alloc_get(MT_DCPL, n_dum, 'dummy', l_dum, k_dum))
     &        call errquit('WFAT: dummy allocation failed', 59, MA_ERR)
         dcpl_mb(k_dum:k_dum+n_dum-1) = czero
      else
         n_acj = 1
         n_ack = 1
         n_dum = 1
      endif

      end
c=================================================

      
c=================================================
c=================================================
      subroutine wfat_erase_acoef(wfat1, l_acg, l_acj, l_ack, l_dum)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      logical wfat1
      integer l_acg, l_acj, l_ack, l_dum

      if (.not. MA_free_heap(l_acg))
     &     call errquit('WFAT: acg deallocation failed', 3, MA_ERR)
      if (wfat1) then
         if (.not. MA_free_heap(l_acj))
     &        call errquit('WFAT: acj deallocation failed', 3, MA_ERR)
         if (.not. MA_free_heap(l_ack))
     &        call errquit('WFAT: ack deallocation failed', 3, MA_ERR)
         if (.not. MA_free_heap(l_dum))
     &        call errquit('WFAT: dum deallocation failed', 3, MA_ERR)
      endif
      
      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_rate(nat, nqn, nch, n_gn, l_gn, k_gn, nfc,
     &     l_fc, k_fc, nrt, l_rt, k_rt)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nat
      integer nqn
      integer nch
      integer n_gn, l_gn, k_gn
      integer nfc, l_fc, k_fc
      integer nrt, l_rt, k_rt

      n_gn = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, n_gn, 'srfc', l_gn, k_gn))
     &     call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
      nfc = nqn * nch
      if (.not. MA_alloc_get(MT_DBL, nfc, 'ffc', l_fc, k_fc))
     &     call errquit('WFAT: ffc allocation failed', 2, MA_ERR)
      nrt = nqn * nch
      if (.not. MA_alloc_get(MT_DBL, nrt, 'p_rate', l_rt, k_rt))
     &     call errquit('WFAT: p_rate allocation failed', 2, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_rate(l_gn, l_fc, l_rt)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_gn, l_fc, l_rt
      
      if (.not. MA_free_heap(l_gn))
     &     call errquit('WFAT: srfc deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_fc))
     &     call errquit('WFAT: ffc deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_rt))
     &     call errquit('WFAT: p_rate deallocation failed', 3, MA_ERR)
      
      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_allgrid(nptp, n_pt, l_pt, k_pt, n_wt,
     &     l_wt, k_wt)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nptp
      integer n_pt, l_pt, k_pt
      integer n_wt, l_wt, k_wt

      n_pt = 3 * nptp
      if (.not. MA_alloc_get(MT_DBL, n_pt, 'grid_pt', l_pt, k_pt))
     &     call errquit('WFAT: grid_pt allocation failed', 2, MA_ERR)
      n_wt = nptp
      if (.not. MA_alloc_get(MT_DBL, n_wt, 'grid_wt', l_wt, k_wt))
     &     call errquit('WFAT: grid_wt allocation failed', 2, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_allgrid(l_pt, l_wt)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_pt, l_wt
      
      if (.not. MA_free_heap(l_pt))
     &     call errquit('WFAT: grid_pt deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_wt))
     &     call errquit('WFAT: grid_wt deallocation failed', 3, MA_ERR)
      
      end
c=================================================
      

c=================================================
c=================================================
      subroutine wfat_prep_wfat1_intgl(wfat1, lmax, natm, nqn, nch, 
     &           nintj, g_j, nintk, g_k, npb, l_pb, k_pb)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
      logical wfat1        ! [input],
      integer lmax         ! [input],
      integer natm         ! [input],
      integer nqn          ! [input],
      integer nch          ! [input],
      integer nintj        ! [output],
      integer g_j          ! [output],
      integer nintk        ! [output],
      integer g_k(2)       ! [output],
      integer npb          ! [output], 
      integer l_pb, k_pb   ! [output], 

      integer nsh, i

      nsh = (1 + lmax)**2
      if (wfat1) then
         nintj = nsh * natm * nqn * 3 * nch
         nintk = nsh * natm * nqn * nch
         npb = 3 * 3 * nch
      else
         nintj = 1     !! 1)
         nintk = 1     !! 2)
         npb = 1       !! 3)
      endif


      if (.not. ga_create(mt_dcpl, nintj, 1, 'intl_j', 0, 0, g_j))
     &        call errquit('WFAT: Cannot create intl_j global array.',
     &        1, GA_ERR)
      call ga_zero(g_j)

      do i = 1, 2
         if (.not. ga_create(mt_dcpl, nintk, 1, 'intl_k', 0, 0, g_k(i)))
     &        call errquit('WFAT: Cannot create intl_k global '//
     &        'array.', 1, GA_ERR)
         call ga_zero(g_k(i))
      enddo
         
      if (.not. MA_alloc_get(MT_DBL, npb, 'polb', l_pb, k_pb))
     &     call errquit('WFAT: polb allocation failed', 59, MA_ERR)

      end
c     NOTES:
c     1) It's actually okay not to allocate the global array g_j when wfat1=.false.. nintj is still computed anyway because it will be input to wfat_prep_wfat1_ti later.
c     2) It's actually okay not to allocate the global array g_k when wfat1=.false.. nintk is still computed anyway because it will be input to wfat_prep_wfat1_ti later.
c     3) polb (static polarizabilities) needs to be allocated anyway because it will be fed into the wfat_calc_intgl subroutine later.
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_wfat1_intgl(g_j, g_k, l_pb)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_j
      integer g_k(2)
      integer l_pb

      integer i

      if (.not. ga_destroy(g_j))
     &     call errquit('WFAT: Cannot destroy intl_j.', 1, GA_ERR)
      do i = 1, 2
         if (.not. ga_destroy(g_k(i)))
     &        call errquit('WFAT: Cannot destroy intl_k.', 1, GA_ERR)
      enddo
      if (.not. MA_free_heap(l_pb))
     &     call errquit('WFAT: polb deallocation failed', 3, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_wfat1_ti(nintj, nintk, l_inj, k_inj,
     &     l_ink1, k_ink1, l_ink2, k_ink2)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nintj             ! [input],
      integer nintk             ! [input],
      integer l_inj, k_inj      ! [output],
      integer l_ink1, k_ink1    ! [output],
      integer l_ink2, k_ink2    ! [output],

      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      
      
      if (.not. MA_alloc_get(MT_DCPL, nintj, 'intgl_j', l_inj, k_inj))          !! 1)
     &     call errquit('WFAT: intgl_j allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nintk, 'intgl_k1',
     &     l_ink1, k_ink1))
     &     call errquit('WFAT: intgl_k1 allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nintk, 'intgl_k2',
     &     l_ink2, k_ink2))
     &     call errquit('WFAT: intgl_k2 allocation failed', 59, MA_ERR)

      dcpl_mb(k_inj:k_inj+nintj-1) = czero
      dcpl_mb(k_ink1:k_ink1+nintk-1) = czero
      dcpl_mb(k_ink2:k_ink2+nintk-1) = czero
      
      end
c     NOTES:
c     1) intgl_j, intgl_k1, and intgl_k2 are allocated regardless of the value
c        of wfat1 because these arrays will be input to the wfat_calc_ti subroutine
c        later. However, when wfat1 is not active, the values of nintj and nintk
c        returned by wfat_prep_wfat1_intgl (that must be called before calling this
c        subroutine) are both 1.
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_wfat1_dgn(wfat1, nsub, maxn, l_dgs, k_dgs,
     &     l_dgss, k_dgss, l_ndgs, k_ndgs, l_dps, k_dps, l_pbs, k_pbs)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      logical wfat1            ! [input]
      integer nsub             ! [input]
      integer maxn             ! [input]
      integer l_dgs, k_dgs     ! [output]
      integer l_dgss, k_dgss   ! [output]
      integer l_ndgs, k_ndgs   ! [output]
      integer l_dps, k_dps     ! [output]
      integer l_pbs, k_pbs     ! [output]

      integer n_dgs, n_dgss, n_ndgs, n_dps, n_pbs

      
      if (wfat1) then
         n_dgs = maxn * nsub
         n_dgss = nsub
         n_ndgs = nsub
         n_dps = maxn**2 * 3 * nsub
         n_pbs = maxn**2 * 6 * nsub
      else
         n_dgs = 1
         n_dgss = 1
         n_ndgs = 1
         n_dps = 1
         n_pbs = 1
      endif

      if (.not. MA_alloc_get(MT_INT, n_dgs, 'degsub_id', l_dgs, k_dgs))
     &     call errquit('WFAT: degsub_id allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_INT, n_dgss, 'degsub_spin',
     &     l_dgss, k_dgss)) call errquit('WFAT: degsub_spin '//
     &     'allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_INT, n_ndgs, 'degsub_n',
     &     l_ndgs, k_ndgs)) call errquit('WFAT: degsub_n allocation '//
     &     'failed', 59, MA_ERR)

      if (.not. MA_alloc_get(MT_DBL, n_dps, 'degsub_dip', l_dps, k_dps))
     &     call errquit('WFAT: degsub_dip allocation failed', 
     &     59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, n_pbs, 'degsub_polb', 
     &     l_pbs, k_pbs)) call errquit('WFAT: degsub_polb allocation '//
     &     'failed', 59, MA_ERR)
      
      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_wfat1_ti(l_inj, l_ink1, l_ink2)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_inj
      integer l_ink1
      integer l_ink2

      if (.not. MA_free_heap(l_inj))
     &     call errquit('WFAT: intgl_j deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_ink1))
     &     call errquit('WFAT: intgl_k1 deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_ink2))
     &     call errquit('WFAT: intgl_k2 deallocation failed', 3, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_erase_wfat1_dgn(l_dgs, l_dgss, l_ndgs, l_dps,
     &     l_pbs)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_dgs, l_dgss, l_ndgs, l_dps, l_pbs

      if (.not. MA_free_heap(l_dgs))
     &     call errquit('WFAT: degsub_id deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_dgss))
     &     call errquit('WFAT: degsub_spin deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_ndgs))
     &     call errquit('WFAT: degsub_n deallocation failed',
     &     3, MA_ERR)

      if (.not. MA_free_heap(l_dps))
     &     call errquit('WFAT: degsub_dip deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_pbs))
     &     call errquit('WFAT: degsub_polb deallocation failed',
     &     3, MA_ERR)

      end
c=================================================

      
c=================================================
c=================================================
      subroutine wfat_prep_resmode(rtdb, rmode)

      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb     ! [input]
      logical rmode    ! [output]

      logical rmode0
      integer ndid


      ndid = ga_nodeid()

      if (rtdb_get(rtdb, 'wfat:res_mode', mt_log, 1, rmode0)) then
         rmode = rmode0
         if (ndid==0) write(luout, 301)
      endif
 301  format(/,
     &      3x, '*** INFO ***', /,
     &      5x, 'WFAT module is run in ''research mode''.')

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_printsw(rtdb, s, p_srfc, p_rate, p_yild,
     &     p_avbeta, p_avgamma, p_tmd, p_tmd_atom)

      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb         ! [input],
      character*(*) s      ! [input],
      logical p_srfc       ! [output],
      logical p_rate       ! [output],
      logical p_yild       ! [output],
      logical p_avbeta     ! [output],
      logical p_avgamma    ! [output],
      logical p_tmd        ! [output],
      logical p_tmd_atom   ! [output],


      if (s/='' .and. s/='inp:') call errquit('WFAT: Bad value for '//
     &     'argument s of wfat_prep_printsw.', 30, UERR)

c     ==== Obtaining print switches ====
      if (.not. rtdb_get(rtdb, 'wfat:'//s//'p_srfc', mt_log, 1, p_srfc))
     &     p_srfc = .false.
      if (.not. rtdb_get(rtdb, 'wfat:'//s//'p_rate', mt_log, 1, p_rate))
     &     p_rate = .false.
      if (.not. rtdb_get(rtdb, 'wfat:'//s//'p_yild', mt_log, 1, p_yild))
     &     p_yild = .false.
      if (.not. rtdb_get(rtdb, 'wfat:'//s//'p_avbeta', mt_log,
     &     1, p_avbeta)) p_avbeta = .false.
      if (.not. rtdb_get(rtdb, 'wfat:'//s//'p_avgamma', mt_log,
     &     1, p_avgamma)) p_avgamma = .false.
      if (.not. rtdb_get(rtdb, 'wfat:'//s//'p_tmd', mt_log, 1, p_tmd))
     &     p_tmd = .false.
      if (.not. rtdb_get(rtdb, 'wfat:'//s//'p_tmd_atom', mt_log, 1,
     &     p_tmd_atom)) p_tmd_atom = .false.

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_dcfield(rtdb, efield)

      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
      integer rtdb              ! [input],
      double precision efield  ! [output],

      integer tg

      if (ga_nodeid() == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Ionizing field', 20, .true.)
         write(luout, '("")')
      endif
      
      if (.not. rtdb_get(rtdb, 'wfat:field:typeg', mt_int, 1, tg))
     &     call errquit('WFAT: Cannot read field:typeg', 1, RTDB_ERR)
      if (tg /= 1) call errquit('WFAT: typeg must be 1 (DC) when '//
     &     'wfat_prep_dcfield is called.', 27, UERR)

      if (.not. rtdb_get(rtdb, 'wfat:field:max', mt_dbl, 1, efield))
     &     call errquit('WFAT: Cannot read field:max', 1, RTDB_ERR)

      if (ga_nodeid() == 0) then
         write(luout, '(5x, "* ", a10, " = ", a)') 'type      ', 'DC'
         write(luout, '(5x, "* ", a10, " = ", f10.3, " a.u.")')
     &        'strength  ', efield
      endif

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_prep_tdfield(rtdb, tmax, ntime, efield)

      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "wfat.fh"
#include "wfat_field_param.fh"
      integer rtdb                 ! [input]
      double precision tmax        ! [output]
      integer ntime                ! [output]
      double precision efield(ntime)  ! [output]

      character*2 tp(max_nfield)
      character*1000 ext_fn, ext_fnb, fname
      integer nf, tg, i, ndid
      double precision dt, t, mx(max_nfield), fr(max_nfield)
      double precision cp(max_nfield), ce(max_nfield), wi(max_nfield)

      double precision wfat_convert_time
      external wfat_convert_time
      

      ndid = ga_nodeid()

      if (.not. rtdb_get(rtdb, 'wfat:field:typeg', mt_int, 1, tg))
     &     call errquit('WFAT: Cannot read field:typeg from RTDB', 70,
     &     RTDB_ERR)
      
c     ==== Get no. field blocks ====
      if (.not. rtdb_get(rtdb, 'wfat:field:n', mt_int, 1, nf))
     &     call errquit('WFAT: Cannot read field:n', 1, RTDB_ERR)

c     ==== Get no. of time points ====
      if (.not. rtdb_get(rtdb, 'wfat:field:ntime', mt_int, 1, ntime))
     &     call errquit('WFAT: Cannot read field:ntime', 1, RTDB_ERR)

c     ==== Get time span ====
      if (.not. rtdb_get(rtdb, 'wfat:field:tmax', mt_dbl, 1, tmax))
     &     call errquit('WFAT: Cannot read field:tmax', 1, RTDB_ERR)
      dt = tmax / (ntime-1)
      
      
      if (tg == 3) then
c        ==== Get laser amplitudes ====
         if (.not. rtdb_get(rtdb, 'wfat:field:max', mt_dbl, nf, mx))
     &        call errquit('WFAT: Cannot read field:max', 1, RTDB_ERR)
c        ==== Get laser envelope type ====
         if (.not. rtdb_cget(rtdb, 'wfat:field:type', nf, tp))
     &        call errquit('WFAT: Cannot read field:type', 1, RTDB_ERR)
c        ==== Get laser frequencies ====
         if (.not. rtdb_get(rtdb, 'wfat:field:freq', mt_dbl, nf, fr))
     &        call errquit('WFAT: Cannot read field:freq', 1, RTDB_ERR)
c        ==== Get laser CEPs ====
         if (.not. rtdb_get(rtdb, 'wfat:field:cep', mt_dbl, nf, cp))
     &        call errquit('WFAT: Cannot read field:cep', 1, RTDB_ERR)
c        ==== Get laser delays ====
         if (.not. rtdb_get(rtdb, 'wfat:field:center', mt_dbl, nf, ce))
     &        call errquit('WFAT: Cannot read field:center', 1,
     &        RTDB_ERR)
c        ==== Get laser durations ====
         if (.not. rtdb_get(rtdb, 'wfat:field:width', mt_dbl, nf, wi))
     &        call errquit('WFAT: Cannot read field:width', 1, RTDB_ERR)

         do i = 1, ntime
            t = (i-1) * dt
            call wfat_eval_laser(nf, tp, mx, fr, cp, ce, wi, t,
     &           efield(i))
         enddo
      else if (tg == 2) then
         if (.not. rtdb_cget(rtdb, 'wfat:ext_fnb', 1, ext_fnb)) 
     &        call errquit('WFAT: Failed reading ext_fnb from rtdb',
     &        1, RTDB_ERR)

         open(1000, file=ext_fnb, action='read', status='old',
     &        form='unformatted')
         do i = 1, ntime
            read(1000) efield(i)
         enddo
         close(1000)
      endif
            

c     ==== Print laser parameters ====
      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Ionizing field', 20, .true.)
         write(luout, '("")')

         if (tg == 3) then
            call wfat_print_laser_info(nf, mx, tp, fr, cp, ce, wi, 
     &           ntime, tmax)
         else if (tg == 2) then
            if (.not. rtdb_cget(rtdb, 'wfat:ext_fn', 1, ext_fn)) 
     &           call errquit('WFAT: Failed reading ext_fn from rtdb',
     &           1, RTDB_ERR)
            write(luout, '(5x, "* ", a16, " = ", a)')
     &           'Laser file name ', trim(adjustl(ext_fn))
         endif

c        ==== Print laser field as a function of time ====
         call util_file_name('wlas', .false., .false., fname)
         open(WF_LAS_UNIT, file=fname, action='write')
         do i = 1, ntime
            t = (i-1) * dt   
            write(WF_LAS_UNIT, '(i4, 2x, es14.6, 2x es14.6)') i,
     &           wfat_convert_time('au',t,'fs'), efield(i)
         enddo
         close(WF_LAS_UNIT)
      endif      

      end
c=================================================


cOLDc=================================================
cOLDc=================================================
cOLD      subroutine wfat_erase_tdfield(l_al, l_wl, l_pl, l_tl, l_dl)
cOLD      consider whether to remove this subroutine
cOLD
cOLD      implicit none
cOLD#include "errquit.fh"
cOLD#include "mafdecls.fh"
cOLD      integer l_al, l_wl, l_pl, l_tl, l_dl
cOLD
cOLD
cOLD      if (.not. ma_free_heap(l_al))
cOLD     &     call errquit('WFAT: alas deallocation failed', 3, MA_ERR)
cOLD      if (.not. ma_free_heap(l_wl))
cOLD     &     call errquit('WFAT: wlas deallocation failed', 3, MA_ERR) 
cOLD      if (.not. ma_free_heap(l_pl))
cOLD     &     call errquit('WFAT: plas deallocation failed', 3, MA_ERR)
cOLD      if (.not. ma_free_heap(l_tl))
cOLD     &     call errquit('WFAT: tlas deallocation failed', 3, MA_ERR)
cOLD      if (.not. ma_free_heap(l_dl))
cOLD     &     call errquit('WFAT: dlas deallocation failed', 3, MA_ERR)
cOLD         
cOLD      end
cOLDc=================================================


c=================================================
c=================================================
      subroutine wfat_prep_channelprop(nch, nerg, ncpos, ncmm, niond,
     &     l_eg, k_eg, l_ps, k_ps, l_cmm, k_cmm, l_idp, k_idp)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nch            ! [input],
      integer nerg           ! [output],
      integer ncpos          ! [output],
      integer ncmm           ! [output],
      integer niond          ! [output],
      integer l_eg, k_eg     ! [output],
      integer l_ps, k_ps     ! [output],
      integer l_cmm, k_cmm   ! [output],
      integer l_idp, k_idp   ! [output],


      nerg = nch
      if (.not. MA_alloc_get(MT_DBL, nerg, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: erg allocation failed', 2, MA_ERR)
      ncpos = 3 * nch
      if (.not. MA_alloc_get(MT_DBL, ncpos, 'chpos', l_ps, k_ps))
     &     call errquit('WFAT: chpos allocation failed', 59, MA_ERR)
      ncmm = 3 * nch
      if (.not. MA_alloc_get(MT_DBL, ncmm, 'chmmt', l_cmm, k_cmm))
     &     call errquit('WFAT: chmmt allocation failed', 2, MA_ERR)
      niond = 3 * nch
      if (.not. ma_alloc_get(MT_DBL, niond, 'iondip', l_idp, k_idp))
     &     call errquit('WFAT: iondip allocation failed', 1, MA_ERR)

      end
c=================================================
      

c=================================================
c=================================================
      subroutine wfat_erase_channelprop(l_eg, l_ps, l_cmm, l_idp)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_eg, l_ps, l_cmm, l_idp
      
      if (.not. MA_free_heap(l_eg))
     &     call errquit('WFAT: erg deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_ps))
     &     call errquit('WFAT: chpos deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_cmm))
     &     call errquit('WFAT: chmmt deallocation failed', 3, MA_ERR)
      if (.not. ma_free_heap(l_idp))
     &     call errquit('WFAT: Deallocating iondip failed', 3, MA_ERR)

      end
c=================================================
