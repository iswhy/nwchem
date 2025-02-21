c=================================================
c=================================================
      logical function task_mewfat(rtdb)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "wfat.fh"
#include "mewfat_qm_param.fh"
#include "wfat_field_param.fh"
      integer rtdb

      character(len=FLD_TP_LEN) las_env(max_nfield)
      character*4096 fname
      character*2 td
      character*7 jobtype
      character*26 date
      logical res_mode, rtdb_stat, os, oc, or, extf, isthere, wfat1
      logical pr_srfc, pr_yild, pr_rate, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom, ismask_i, ismask_e, int_nosave
      integer gqual, lmax
      integer ndid, nqn, nlas, ntime, nsh, nintg
      integer nerg, ncpos, ncmm, niond, cx, cy, nocmo0, nocmo1, npcf
      integer nrcf, nqcf, nscf, nn, ma
      integer l_npa, l_mpa, l_al, l_ftl, l_wl, l_pl, l_tl, l_dl
      integer l_ieg, l_ps, l_cmm, l_idp, l_in, l_mo0, l_mo1, l_pcf
      integer l_rcf, l_qcf, l_scf, l_eg0, l_eg1, l_ed0, l_ed1
      integer k_npa, k_mpa, k_al, k_ftl, k_wl, k_pl, k_tl, k_dl
      integer k_ieg, k_ps, k_cmm, k_idp, k_in, k_mo0, k_mo1, k_pcf
      integer k_rcf, k_qcf, k_scf, k_eg0, k_eg1, k_ed0, k_ed1
      integer g_il
      double precision tmax, fieldmag, smh, mask_param(2), ergmask
      double precision amsk, bmsk

      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )

      
      
      ndid = ga_nodeid()
      task_mewfat = .false.

      
c     ================================
c     ==== Initiate MEWFAT module ====
c     ================================
      call mewfat_init(rtdb)


c     ==== Task title ====
      call ga_sync()
      jobtype = 'ME'
      call wfat_prep_tasktitle(jobtype)
      

c     ==== Reading general parameters ====
      call wfat_prep_general(rtdb, gqual, lmax)
      

c     ==== Checking the numerical density ====



c     ==== Obtaining parabolic quantum numbers ====
      call wfat_prep_paraqn(rtdb, nqn, l_npa, k_npa, l_mpa, k_mpa)


c     ==== Obtaining WFAT(1) information ====
cforlater      call wfat_prep_wfat1(rtdb, wfat1, nlgl)
      wfat1 = .false.
      

c     ==== Check if research mode is enabled ====
      call wfat_prep_resmode(rtdb, res_mode)
      

c     ==== Obtaining print switches ====
      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)
      if (pr_tmd .or. pr_tmd_atom)
     &     call errquit('MEWFAT: At the moment, TMD calculation is '//
     &     'not yet implemented within MEWFAT', 1, UERR)


c     ==== When TD field is true ====
      if (.not. rtdb_cget(rtdb, 'wfat:tdfield', 1, td))
     &     td = 'no'
      if (td=='td') then
         call wfat_prep_tdfield(rtdb, nlas, l_al, k_al, las_env,
     &        l_wl, k_wl, l_pl, k_pl, l_tl, k_tl, l_dl, k_dl, ntime,
     &        tmax, extf)
      else if (td=='ti') then
         call wfat_prep_field(rtdb, fieldmag)
      else if (td=='no') then
         fieldmag = 0.0d0
      endif


c     ==== Get channel properties ====
      call wfat_prep_channelprop(1, nerg, ncpos, ncmm, niond,
     &     l_ieg, k_ieg, l_ps, k_ps, l_cmm, k_cmm, l_idp, k_idp)
      call mewfat_prep_channelprop(1, l_eg0, k_eg0, l_eg1, k_eg1,
     &     l_ed0, k_ed0, l_ed1, k_ed1)
      call mewfat_channel_property(rtdb, geom_hd, bas_hd, nbas, gqual,
     &     wf_msm, nel0, nel1, g_mov0, g_mov1, charge0, charge1, 
     &     multip0, multip1, dbl_mb(k_eg0), dbl_mb(k_eg1), 
     &     dbl_mb(k_ieg), dbl_mb(k_ed0), dbl_mb(k_ed1), dbl_mb(k_ps), 
     &     dbl_mb(k_cmm), dbl_mb(k_idp))
      
      

c      call wfat_test_dtrm(rtdb, geom_hd, bas_hd, nbas0, nspin0, nel, 
c     &                   g_mov, charge0, multip0, dbl_mb(k_moerg), tol2e, 
c     &                   oskel, ixcfit, iolgc)


c     ==== Print Dyson orbital ====
      call mewfat_print_dyson(rtdb, geom_hd, bas_hd, nel0, spinx,
     &     g_mov0, g_mov1, oskel)
      

c     ==== Diagonalize total dipole matrix ====
      if (rtdb_get_info(rtdb, 'mewfat:dipdiag', ma, nn, date))
     &     call mewfat_dipdiag(rtdb, geom_hd, bas_hd, nn)      
      
      
c     ==== Calculate determinant coefficients and get MO's ====
      call mewfat_alloc_ocmo(nbas, nel0, nel1, nocmo0, nocmo1,
     &     l_mo0, k_mo0, l_mo1, k_mo1)
      call mewfat_alloc_detcoef(spinx, spiny, nel0, cx, cy, npcf, nrcf, 
     &     nqcf, nscf, l_pcf, k_pcf, l_rcf, k_rcf, l_qcf, k_qcf, l_scf,
     &     k_scf)
      call mewfat_calc_detcoef(geom_hd, bas_hd, nbas, spinx, spiny,
     &     nel0, nel1, cx, cy, g_mov0, g_mov1, oskel, dbl_mb(k_mo0),
     &     dbl_mb(k_mo1), dbl_mb(k_pcf), dbl_mb(k_rcf), dbl_mb(k_qcf),
     &     dbl_mb(k_scf))
      if (ndid == 0)
     &     call mewfat_print_detcoef(spinx, nel0, cx, cy, dbl_mb(k_pcf),
     &     dbl_mb(k_rcf), dbl_mb(k_qcf), dbl_mb(k_scf))

      
c     ==== Create MEWFAT integrals array ====
      nsh = (1 + lmax)**2
      nintg = nsh * natm * nqn
      if (.not. ga_create(mt_dcpl, nintg, 1, 'intgl', 0, 0, g_il))
     &     call errquit('WFAT: Cannot create intgl global array.',
     &     1, GA_ERR)
      call ga_zero(g_il)

      
c     ==== Get mask parameters ====
      call wfat_get_masks(rtdb, dbl_mb(k_ieg), ismask_i, mask_param, 
     &     ismask_e, ergmask, amsk, bmsk)


      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('MEWFAT: Cannot read wfat:intgl_name from '//
     &     'RTDB.', 1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      

c     ==== Calculate or store MEWFAT integrals ====
      smh = 0.0d0
      call util_get_rtdb_state(rtdb, os, oc, or)
      if (lmax >= 0) then
         if (os .and. .not. or) then
cNEWRESTART      if (.not. isthere .or. os) then
c           ==== Calculate ME-WFAT integrals on 'start' startup mode ====
            call mewfat_calc_intgl(rtdb, geom_hd, bas_hd, nbas, gqual, 
     &           spinx, spiny, nel0, nel1, dbl_mb(k_mo0), dbl_mb(k_mo1),
     &           dbl_mb(k_pcf), dbl_mb(k_rcf), dbl_mb(k_qcf), 
     &           dbl_mb(k_scf), zion, dbl_mb(k_ps), dbl_mb(k_ieg), 
     &           wf_msm, nqn, int_mb(k_npa), int_mb(k_mpa), lmax, smh, 
     &           ismask_i, mask_param, wfphase, res_mode, g_il)
         
c           ==== Store ME-WFAT integrals freshly calculated ====
            if (.not. rtdb_get(rtdb, 'wfat:intgl_nosave', mt_log, 1,
     &           int_nosave)) call errquit('MEWFAT: Cannot read '//
     &           'intgl_nosave from RTDB', 1, RTDB_ERR)
            if (ndid == 0 .and. .not. int_nosave) then
               call mewfat_store_intgl(rtdb, nsh, natm, nqn, 1, 
     &              dbl_mb(k_eg0), dbl_mb(k_eg1), dbl_mb(k_ieg),
     &              dbl_mb(k_ed0), dbl_mb(k_ed1), dbl_mb(k_ps),
     &              dbl_mb(k_cmm), dbl_mb(k_idp), g_il)
            else if (ndid == 0 .and. int_nosave) then
               write(luout, '("ATTENTION: MEWFAT integrals are not '//
     &              'saved.")')
            endif
            call ga_sync()
            
         else if (.not. os .and. or) then
cNEWRESTART      else
c           ==== Binary-load MEWFAT integrals on 'restart' startup mode ====
            if (ndid == 0) then
               call mewfat_load_intgl(rtdb, nsh, natm, nqn, 1,
     &              spinx, dbl_mb(k_ieg), dbl_mb(k_ps), dbl_mb(k_cmm),
     &              dbl_mb(k_idp), g_il)
            endif
            call ga_sync()
            
            call ga_brdcst(msg_wfat_eg, dbl_mb(k_ieg),
     &                     ma_sizeof(MT_DBL, 1, MT_BYTE), 0)
            call ga_brdcst(msg_wfat_cmm, dbl_mb(k_cmm),
     &                     ma_sizeof(MT_DBL, 3*1, MT_BYTE), 0)
         else
            call errquit('MEWFAT: MEWFAT module requires the startup '//
     &           'directive to be either ''start'' or ''restart''.',
     &           1, UERR)
         endif
      else
         if (.not. os .and. or) 
     &        call errquit('MEWFAT: When performing WFAT calculation '//
     &        'without partial wave expansion, the job must be a '//
     &        'start type. Restart is not allowed', 38, UERR)
      endif      


c     ==== Printings ====
      if (lmax >= 0 .and. ndid == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &        call errquit('WFAT: Cannot allocate intgl.', 2, MA_ERR)
         call ga_get(g_il, 1, nintg, 1, 1, dcpl_mb(k_in), nintg)
         if (wfat1) then
         else
            call wfat_print_wfat0_ints(1, nqn, natm, lmax, 
     &           int_mb(k_npa), int_mb(k_mpa), dcpl_mb(k_in))
         endif
cOLD         call wfat_print_wfat1_ints(1, nqn, natm, lmax, int_mb(k_npa),
cOLD     &     int_mb(k_mpa), wfat1, dcpl_mb(k_in), czero, czero, czero)


c        ============================================
c        ==== Calculate and print key quantities ====
c        ============================================
         if (td=='td') then
            if (pr_srfc .or. pr_yild .or. pr_avbeta .or. pr_avgamma)
     &           call mewfat_calc_td(rtdb, geom_hd, nel0, natm, pr_srfc,
     &           pr_yild, pr_avbeta, pr_avgamma, nlas, las_env, 
     &           dbl_mb(k_al), dbl_mb(k_wl), dbl_mb(k_pl), dbl_mb(k_tl),
     &           dbl_mb(k_dl), ntime, tmax, extf, 1,
     &           spinx, dbl_mb(k_ieg), dbl_mb(k_cmm), nqn,
     &           int_mb(k_npa), int_mb(k_mpa), lmax, zion,
     &           dcpl_mb(k_in))
         else
            if (pr_srfc .or. pr_rate)
     &           call mewfat_calc_ti(rtdb, geom_hd, nel0, natm, pr_srfc,
     &           pr_rate, fieldmag, 1, spinx, dbl_mb(k_ieg), 
     &           dbl_mb(k_cmm), nqn, int_mb(k_npa), int_mb(k_mpa),
     &           lmax, zion, dcpl_mb(k_in))
         endif
         
         
         if (.not. MA_free_heap(l_in))
     &        call errquit('WFAT: intgl deallocation failed', 3, MA_ERR)
         rtdb_stat = rtdb_parallel(.true.)
      else if (lmax < 0) then
         call mewfat_calc_ti_nopw(rtdb, geom_hd, bas_hd, pr_srfc, 
     &        pr_rate, fieldmag, nbas, spinx, spiny, nel0, nel1, 
     &        dbl_mb(k_mo0), dbl_mb(k_mo1), dbl_mb(k_pcf), 
     &        dbl_mb(k_rcf), dbl_mb(k_qcf), dbl_mb(k_scf), 1, 0, 
     &        dbl_mb(k_ieg), dbl_mb(k_cmm), dbl_mb(k_ps), nqn, 
     &        int_mb(k_npa), int_mb(k_mpa), zion, gqual, wf_msm, smh, 
     &        ismask_e, ergmask, amsk, bmsk, ismask_i, mask_param)
      endif
      

c     ==== Deallocations ====
      call wfat_erase_paraqn(l_npa, l_mpa)
      call mewfat_erase_channelprop(l_eg0, l_eg1, l_ed0, l_ed1)
      call mewfat_dealloc_ocmo(l_mo0, l_mo1)
      call mewfat_dealloc_detcoef(l_pcf, l_rcf, l_qcf, l_scf)
      if (td=='td')
     &     call wfat_erase_tdfield(l_al, l_wl, l_pl, l_tl, l_dl)
      call wfat_erase_channelprop(l_ieg, l_ps, l_cmm, l_idp)
      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 1, GA_ERR)


c     ==== Clean up and terminate task ====
      call wfat_close_task('ME')
      call mewfat_terminate(rtdb)
      

      task_mewfat = .true.
      
      
      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_calc_td(rtdb, geom, n0, nat, pr_srfc, pr_yild,
     &           pr_avbeta, pr_avgamma, nlas, ftlas, alas, wlas,
     &           plas, tlas, dlas, nt, tmax, extf, nch, spinid, erg,
     &           chmmt, nqn, nid, mid, lmax, zi, intgl)

c     DESCRIPTION:
      
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "wfat_field_param.fh"
      integer rtdb              ! [input]
      integer geom
      integer n0(2)             ! [input], no. of electrons in the neutral.
      integer nat
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_yild                 ! [input], if .true. then print ion yield.
      logical pr_avbeta               ! [input], if .true. then print beta-averaged ion yield.
      logical pr_avgamma              ! [input], if .true. then print gamma-averaged ion yield.
      integer nlas
      character(len=FLD_TP_LEN) ftlas(nlas)
      double precision alas(nlas)   ! Peak intensity, TW/cm^2                              
      double precision wlas(nlas)   ! Wavelength, nm                                       
      double precision plas(nlas)   ! Carrier-envelope phase, rad                          
      double precision tlas(nlas)   ! Delay (temporal position of the envelope peak), fs   
      double precision dlas(nlas)   ! Field FWHM, fs                                       
      double precision nt
      double precision tmax
      logical extf
      integer nch                     ! [input], number of ionized channels.
      integer spinid(nch)             ! [input], spin value (1 or 2) of ionized channels.
      double precision erg(nch)       ! [input]
      double precision chmmt(3, nch)  ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax                    ! [input]
      double precision zi             ! [input]
      double complex intgl(1:*)       ! [input]

      integer nsh, nintg, dummy_int, l_in, k_in
      double precision nall

      logical wfat1, wfat1_exc
      integer ndgsub, max_ndgsub, dgs, ndgs
      double precision dps, pbs
      double complex intgl_j, intgl_k1, intgl_k2, czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      

      nsh = (1 + lmax)**2
      nintg = nsh * nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &     call errquit('MEWFAT: Cannot allocate intgl.', 2, MA_ERR)
      nall = dble(sum(n0))
      dcpl_mb(k_in : k_in+nintg-1) = sqrt(nall) * intgl(1:nintg)     ! 1)
c     1) The prefactor \sqrt{N} actually appears in the expression for ionization
c     rates, not in the asymptotic coefficients or structure factors. However,
c     this line will include this prefactor in the calculated values of asymptotic
c     coefficients and structure factors. This is done because the same subroutine to
c     calculate ionization rates, structure factors, and asymptotic coefficients
c     are also used in the OE-WFAT codes, in other words, for facilitating code
c     reusability.


      dummy_int = 0    ! Channel ID.
c     WARNING!!!!: Update the arguments of wfat_calc_td below!!!
cOLD      call wfat_calc_td(rtdb, geom, pr_srfc, pr_yild, 
cOLD     &     pr_avbeta, pr_avgamma, nlas, ftlas, alas, wlas, 
cOLD     &     plas, tlas, dlas, nt, tmax, extf, nch, dummy_int, spinid, 
cOLD     &     erg, chmmt, nqn, nid, mid, lmax, zi, dcpl_mb(k_in))


      
c     ==== These dummy parameters will need to be properly computed ====
c     ====      later when WFAT(1) is implemented for ME-WFAT       ====
      wfat1 = .false.
      ndgsub = 1
      max_ndgsub = 1
      dgs = 1
      ndgs = 1
      dps = 0.0d0
      pbs = 0.0d0
      wfat1_exc = .false.
      intgl_j = czero
      intgl_k1 = czero
      intgl_k2 = czero
c     ==================================================================
      
cOLD      call wfat_calc_td(rtdb, geom, pr_srfc, pr_yild, 
cOLD     &     pr_avbeta, pr_avgamma, nlas, ftlas, alas, wlas, 
cOLD     &     plas, tlas, dlas, nt, tmax, extf, nch, dummy_int, spinid, 
cOLD     &     erg, chmmt, nqn, nid, mid, lmax, zi, dcpl_mb(k_in), 
cOLD     &     wfat1, nb, ns, orbid, dpmo, pb, chpos, intgl_j, intgl_k1,
cOLD     &     intgl_k2)
      call wfat_calc_td(rtdb, geom, pr_srfc, pr_yild, 
     &     pr_avbeta, pr_avgamma, nlas, ftlas, alas, wlas, 
     &     plas, tlas, dlas, nt, tmax, extf, nch, dummy_int, spinid, 
     &     erg, chmmt, nqn, nid, mid, lmax, zi, dcpl_mb(k_in), 
     &     wfat1, ndgsub, max_ndgsub, dgs, ndgs, dps, pbs, wfat1_exc, 
     &     intgl_j, intgl_k1, intgl_k2)
      

      if (.not. ma_free_heap(l_in))
     &     call errquit('MEWFAT: Cannot deallocate intgl', 31, MA_ERR)
      
      end
c=================================================


c===============================================================
c===============================================================
      subroutine mewfat_calc_ti(rtdb, geom, n0, nat, pr_srfc, pr_rate, 
     &           fieldm, nch, spinid, erg, chmmt, nqn, nid, mid,
     &           lmax, zi, intgl)

c     DESCRIPTION:
      
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input]
      integer geom
      integer n0(2)             ! [input], no. of electrons in the neutral.
      integer nat
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_rate                 ! [input], if .true. then print ion yield.
      double precision fieldm
      integer nch                     ! [input], number of ionized channels.
      integer spinid(nch)             ! [input], spin value (1 or 2) of ionized channels.
      double precision erg(nch)       ! [input]
      double precision chmmt(3, nch)  ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax                    ! [input]
      double precision zi             ! [input]
      double complex intgl(1:*) ! [input]

      integer nsh, nintg, dummy_int, l_in, k_in
      double precision nall

      logical wfat1, wfat1_exc
      integer ndgsub, max_ndgsub, dgs, ndgs
      double precision dps, pbs
      double complex intgl_j, intgl_k1, intgl_k2, czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      

      nsh = (1 + lmax)**2
      nintg = nsh * nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &     call errquit('MEWFAT: Cannot allocate intgl.', 2, MA_ERR)
      nall = dble(sum(n0))
      dcpl_mb(k_in : k_in+nintg-1) = sqrt(nall) * intgl(1:nintg)     ! 1)
c     1) The prefactor \sqrt{N} actually appears in the expression for ionization
c     rates, not in the asymptotic coefficients or structure factors. However,
c     this line will include this prefactor in the calculated values of asymptotic
c     coefficients and structure factors. This is done because the same subroutine to
c     calculate ionization rates, structure factors, and asymptotic coefficients
c     are also used in the OE-WFAT codes, in other words, for facilitating code
c     reusability.

      dummy_int = 0         ! Channel ID.
cOLD      call wfat_calc_ti(rtdb, geom, pr_srfc, pr_rate, fieldm, 
cOLD     &     1, dummy_int, spinid, erg, chmmt, nqn, nid, mid, lmax,
cOLD     &     zi, dcpl_mb(k_in))


c     ==== These dummy parameters will need to be properly computed ====
c     ====      later when WFAT(1) is implemented for ME-WFAT       ====
      wfat1 = .false.
      ndgsub = 1
      max_ndgsub = 1
      dgs = 1
      ndgs = 1
      dps = 0.0d0
      pbs = 0.0d0
      wfat1_exc = .false.
      intgl_j = czero
      intgl_k1 = czero
      intgl_k2 = czero
c     ==================================================================
      
cOLD      call wfat_calc_ti(rtdb, geom, pr_srfc, pr_rate, fieldm, 
cOLD     &     1, dummy_int, spinid, erg, chmmt, nqn, nid, mid, lmax,
cOLD     &     zi, dcpl_mb(k_in), wfat1, nb, ns, orbid, dpmo, pb, chpos, 
cOLD     &     intgl_j, intgl_k1, intgl_k2)
      call wfat_calc_ti(rtdb, geom, pr_srfc, pr_rate, fieldm, 
     &     1, dummy_int, spinid, erg, chmmt, nqn, nid, mid, lmax,
     &     zi, dcpl_mb(k_in), wfat1, ndgsub, max_ndgsub, dgs, ndgs, dps,
     &     pbs, wfat1_exc, intgl_j, intgl_k1, intgl_k2)

      
      if (.not. ma_free_heap(l_in))
     &     call errquit('MEWFAT: Cannot deallocate intgl', 83, MA_ERR)
      
      end
c===============================================================
