c===============================================================
c===============================================================
      logical function task_wfat(rtdb)

      implicit none
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb

      character*32 theory
      logical task_oewfat, task_wfat_2a1e
      external task_oewfat, task_wfat_2a1e

      if (.not. rtdb_cget(rtdb, 'wfat:theory', 1, theory))
     &     call errquit('WFAT: Error reading theory from rtdb',
     &     1, RTDB_ERR)

      if (trim(adjustl(theory)) == '2a1e') then
         task_wfat = task_wfat_2a1e(rtdb)
      else
         task_wfat = task_oewfat(rtdb)
      endif
      
      end
c===============================================================

      
c===============================================================
c===============================================================
      logical function task_oewfat(rtdb)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "wfat.fh"
#include "wfat_qm_param.fh"
#include "wfat_field_param.fh"
      integer rtdb

      character(len=FLD_TP_LEN) las_env(max_nfield)
      character*2 td
      character*7 jobtype
      logical pr_srfc, pr_yild, pr_rate, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom
      logical rtdb_stat, wfat1, res_mode, os, oc, or, int_nosave
      logical dd, extf, corb_print, perturb, ismask_i, ismask_e
      character*32 theory
      character*2 prtb_type
      integer ndid, gqual, i, lmax, nlgl, is, npb, ndpo, ncf, nch1, nrd
      integer l_nd, l_or0, l_sp, l_or, l_npa, l_mpa, l_eg, l_cf, l_rd
      integer l_ps, l_cmm, l_idp, l_in, l_al, l_wl, l_pl, l_tl
      integer l_dl, l_pb, l_dpo, l_dps, l_pbs, l_inj, l_ink1, l_ink2
      integer l_dgs, l_dgss, l_ndgs, l_w1x
      integer k_nd, k_or0, k_sp, k_or, k_npa, k_mpa, k_eg, k_cf, k_rd
      integer k_ps, k_cmm, k_idp, k_in, k_al, k_wl, k_pl, k_tl
      integer k_dl, k_pb, k_dpo, k_dps, k_pbs, k_inj, k_ink1, k_ink2
      integer k_dgs, k_dgss, k_ndgs, k_w1x
      integer g_il, g_j, g_k(2), nlas, ntime
      integer nch, nqn, nerg, ncpos, ncmm, niond, nsh, nintg, nintj
      integer nintk, prtb_ord, ndgsub, max_ndgsub, nw1x
      double precision fieldmag, smh, prtb_bt, prtb_gm, ctr(3)
      double precision mask_param(2), ergmask, amsk, bmsk, tmax

      
      ndid = ga_nodeid()
      task_oewfat = .false.

      
c     ==============================
c     ==== Initiate WFAT module ====
c     ==============================
      call wfat_init(rtdb)

      res_mode = .false.

      
c     ==== Check whether the current theory is DFT ====
      if(.not. rtdb_cget(rtdb,'wfat:theory',1,theory))
     &     call errquit('WFAT: Error reading previous task theory '//
     &     'from RTDB', 1, INPUT_ERR)
      if (trim(adjustl(theory)) /= 'RHF' .and.
     &    trim(adjustl(theory)) /= 'UHF') then
         if (ndid==0) write(luout, '("WFAT: Error because previous '//
     &        'SCF type is ", a, " instead of RHF.")')
     &        trim(theory)
         call errquit('      WFAT task operation is specified while '//
     &        'the theory is not DFT. At the moment WFAT can only be '//
     &        'used in conjunction with single determinant methods '//
     &        'to obtain quantities such as molecular orbitals and '//
     &        'effective potential. Thus terminating NWChem ...',
     &        1, INPUT_ERR)
      endif
  

      
c     ==== Task title ====
      call ga_sync()
      jobtype = 'OE'
      call wfat_prep_tasktitle(jobtype)


c     ==== Reading general parameters ====
      call wfat_prep_general(rtdb, gqual, lmax)

      
c     ==== Checking the numerical density ====
      if (.not. MA_alloc_get(MT_DBL, 2, 'numden', l_nd, k_nd))
     &     call errquit('WFAT: numden allocation failed', 59, MA_ERR)
      call dfill(2, 0.0d0, dbl_mb(k_nd), 1)
      call wfat_numdens(bas_hd, geom_hd, nspin0, gqual,
     &        dbl_mb(k_nd : k_nd+1))

      
c     ==== Obtaining ionization channels ====
      call wfat_prep_channels(rtdb, nel, nch, l_or0, k_or0, l_sp, k_sp,
     &     l_or, k_or)


c     ==== Obtaining parabolic quantum numbers ====
      call wfat_prep_paraqn(rtdb, nqn, l_npa, k_npa, l_mpa, k_mpa)
      

c     ==== Obtaining WFAT(1) information ====
      call wfat_prep_wfat1(rtdb, nch, nqn, wfat1, nlgl, nw1x, l_w1x,
     &     k_w1x)

      
c     ==== Check if research mode is enabled ====
      call wfat_prep_resmode(rtdb, res_mode)


c     ==== Obtaining print switches ====
      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)

      
c     ==== Whether using TD field or no ====
      if (.not. rtdb_cget(rtdb, 'wfat:tdfield', 1, td)) td = 'no'
      if (td=='td') then
         call wfat_prep_tdfield(rtdb, nlas, l_al, k_al, las_env,
     &        l_wl, k_wl, l_pl, k_pl, l_tl, k_tl, l_dl, k_dl, ntime,
     &        tmax, extf)
      else if (td=='ti') then
         call wfat_prep_field(rtdb, fieldmag)
      else if (td=='no') then
         fieldmag = 0.0d0
      endif


c     ==== Perturbation theory analysis ====
      if (.not. rtdb_get(rtdb, 'wfat:perturb', mt_log, 1, perturb))
     &     perturb = .false.
      if (td=='ti' .and. perturb)
     &     call wfat_calc_perturb(rtdb, fieldmag, pr_rate,
     &     pr_avbeta, pr_avgamma, pr_tmd, pr_tmd_atom, res_mode,
     &     nch, int_mb(k_or), int_mb(k_sp))


cfordebug     ==== Some tests ====
cfordebug      call wfat_test_numerfexch(rtdb)
      

c     ==== Allocate channel properties ====
      call wfat_prep_channelprop(nch, nerg, ncpos, ncmm, niond,
     &     l_eg, k_eg, l_ps, k_ps, l_cmm, k_cmm, l_idp, k_idp)
      

c     ==== Calculate (or load) WFAT integrals ====
      nsh = (1 + lmax)**2
      nintg = nsh * natm * nqn * nch
      if (.not. ga_create(mt_dcpl, nintg, 1, 'intgl', 0, 0, g_il))
     &     call errquit('WFAT: Cannot create intgl global array.',
     &     1, GA_ERR)
      call ga_zero(g_il)

      
c     ==== Preparing WFAT(1) integrals ====
      call wfat_prep_wfat1_intgl(wfat1, lmax, natm, nqn, nch, nintj, 
     &     g_j, nintk, g_k, npb, l_pb, k_pb)


c     ==== Printing the corrected ionizing orbital ====
      if (.not. rtdb_get(rtdb, 'wfat:corb_print', mt_log, 
     &     1, corb_print)) corb_print = .false.
      if (corb_print)
     &     call wfat_print_orb1(rtdb, nch, int_mb(k_or), int_mb(k_sp))

      
c     ==== GTO coefficients of ionized orbitals and their properties ====
      call wfat_prep_cf(nch, nbas0, wfat1, ncf, nch1, l_cf, k_cf)
      call wfat_orb_property(rtdb, nch, int_mb(k_or), int_mb(k_sp), 
     &     dbl_mb(k_cf), dbl_mb(k_eg), dbl_mb(k_ps), dbl_mb(k_cmm),
     &     dbl_mb(k_idp))
      if (wfat1) then
         call wfat_orb1_property(nch, int_mb(k_sp), int_mb(k_or),
     &        dbl_mb(k_cf+nbas0*nch), dbl_mb(k_pb))
         if (ndid == 0)
     &        call wfat_print_w1x(nch, int_mb(k_or0), int_mb(k_sp), nqn,
     &        int_mb(k_npa), int_mb(k_mpa), log_mb(k_w1x))
         call ga_sync()
      endif

      
c     ==== Get the reduced density matrix ====
      nrd = nbas0 * nbas0 * nspin0
      if (.not. MA_alloc_get(MT_DBL, nrd, 'rdm', l_rd, k_rd))
     &     call errquit('WFAT: rdm allocation failed', 81, MA_ERR)
      call wfat_get_rdm(dbl_mb(k_rd))   ! Consider creating a more generic RDM calculator.


c     ==== Get mask parameters ====
      call wfat_get_masks(rtdb, dbl_mb(k_eg), ismask_i, mask_param, 
     &     ismask_e, ergmask, amsk, bmsk)
      if (ismask_i .and. ndid==0) call wfat_print_mask(mask_param(1),
     &     mask_param(2))

      
c     ==== Print grid info ====#
c      call wfat_print_grid()
      
      
c     ==== Calculate or store OEWFAT integrals ====
      smh = 0.0d0
      call util_get_rtdb_state(rtdb, os, oc, or)
      if (lmax >= 0) then
         if (os .and. .not. or) then
c           ==== Calculate OE-WFAT integrals on 'start' startup mode ====
            call wfat_calc_intgl(rtdb, geom_hd, bas_hd, nch, 
     &           int_mb(k_sp), nqn, int_mb(k_npa), int_mb(k_mpa), lmax, 
     &           gqual, zion, wf_msm, smh, res_mode, wfat1, nlgl, 
     &           log_mb(k_w1x), nbas0, nspin0, dbl_mb(k_rd), 
     &           dbl_mb(k_cf), dbl_mb(k_eg), dbl_mb(k_ps), 
     &           dbl_mb(k_cmm), dbl_mb(k_idp), ismask_i, mask_param,
     &           g_il, g_j, g_k)
            
c           ==== Store the OE-WFAT integrals freshly calculated ====
            if (.not. rtdb_get(rtdb, 'wfat:intgl_nosave', mt_log, 1,
     &           int_nosave)) call errquit('WFAT: Cannot read '//
     &           'intgl_nosave from RTDB', 1, RTDB_ERR)
            if (ndid == 0 .and. .not. int_nosave) then
               call wfat_store_intgl(rtdb, nsh, natm, nqn, nch,
     &              dbl_mb(k_eg), dbl_mb(k_ps), dbl_mb(k_cmm),
     &              dbl_mb(k_idp), g_il)
            else if (ndid == 0 .and. int_nosave) then
               write(luout, '("ATTENTION: WFAT integrals are not '//
     &              'saved.")')
            endif
            call ga_sync()
            
         else if (.not. os .and. or) then
c           ==== Binary-load WFAT integrals on 'restart' startup mode ====
            if (ndid == 0) then
               call wfat_load_intgl(rtdb, nsh, natm, nqn, nch,
     &              dbl_mb(k_eg), dbl_mb(k_ps), dbl_mb(k_cmm),
     &              dbl_mb(k_idp), g_il)
               call wfat_print_orb_property(nch, int_mb(k_or), 
     &              int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_ps),
     &              dbl_mb(k_cmm), dbl_mb(k_idp))
            endif
            call ga_sync()
            
            call ga_brdcst(msg_wfat_eg, dbl_mb(k_eg),
     &                     ma_sizeof(MT_DBL, nch, MT_BYTE), 0)
            call ga_brdcst(msg_wfat_cmm, dbl_mb(k_cmm),
     &                     ma_sizeof(MT_DBL, 3*nch, MT_BYTE), 0)
         else
            call errquit('WFAT: WFAT module requires the startup '//
     &           'directive to be either ''start'' or ''restart''.',
     &           1, UERR)
         endif
      else
         if (.not. os .and. or) 
     &        call errquit('WFAT: When performing WFAT calculation '//
     &        'without partial wave expansion, the job must be a '//
     &        'start type. Restart is not allowed', 38, UERR)
      endif


c     ==== Take care of the dipole moments and polarizabilities ====
c     ====      in case of WFAT(1) and degenerate channels      ====
      if (wfat1) then
c        ==== Rearrange ionization channel IDs ====
c        ====    into degenerate subspaces     ====
         call wfat_dgn_nsub_ch(nspin0, nch, int_mb(k_or),
     &        int_mb(k_sp), ndgsub, max_ndgsub)
         call wfat_prep_wfat1_dgn(wfat1, ndgsub, max_ndgsub, l_dgs,
     &        k_dgs, l_dgss, k_dgss, l_ndgs, k_ndgs, l_dps, k_dps,
     &        l_pbs, k_pbs)
         call wfat_dgn_sub(nch, int_mb(k_or), int_mb(k_sp), ndgsub, 
     &        max_ndgsub, int_mb(k_dgs), int_mb(k_dgss),
     &        int_mb(k_ndgs))

         ndpo = nbas0 * nbas0 * nspin0 * 3
         if (.not. MA_alloc_get(MT_DBL, ndpo, 'dpmo', l_dpo, k_dpo))
     &        call errquit('WFAT: dpmo allocation failed', 59, MA_ERR)
         ctr = 0.0d0
         call wfat_calc_dpmo(ctr, dbl_mb(k_dpo))

c        ==== Obtain degenerate subspaces dipole moment ====
         call wfat_dpmo_dgn_sub(nch, int_mb(k_or), ndgsub, 
     &        max_ndgsub, int_mb(k_dgs), int_mb(k_dgss), 
     &        int_mb(k_ndgs), nbas0, nspin0, dbl_mb(k_dpo),
     &        dbl_mb(k_dps))

c        ==== Obtain degenerate subspaces polarizability ====
         call wfat_polb_dgn_sub(nch, int_mb(k_or), ndgsub, 
     &        max_ndgsub, int_mb(k_dgs), int_mb(k_dgss), 
     &        int_mb(k_ndgs), nbas0, nspin0, dbl_mb(k_dpo),
     &        dbl_mb(k_moerg), dbl_mb(k_pbs))
      endif
         

c     ==== Ionization calculations and printings ====
      if (lmax >= 0 .and. ndid == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &        call errquit('WFAT: Cannot allocate intgl.', 2, MA_ERR)
         call ga_get(g_il, 1, nintg, 1, 1, dcpl_mb(k_in), nintg)

c        ==== When WFAT(1) is in use ====
         call wfat_prep_wfat1_ti(nintj, nintk, l_inj, k_inj, l_ink1, 
     &        k_ink1, l_ink2, k_ink2)
         if (wfat1) then
            call ga_get(g_j, 1, nintj, 1, 1, dcpl_mb(k_inj), nintj)
            call ga_get(g_k(1), 1, nintk, 1, 1, dcpl_mb(k_ink1), nintk)
            call ga_get(g_k(2), 1, nintk, 1, 1, dcpl_mb(k_ink2), nintk)
         endif
         
c        ==== Printing the WFAT integrals ====
         if (wfat1) then
            call wfat_print_wfat1_ints(nch, nqn, natm, lmax, 
     &           int_mb(k_npa), int_mb(k_mpa), wfat1, dcpl_mb(k_in), 
     &           dcpl_mb(k_inj), dcpl_mb(k_ink1), dcpl_mb(k_ink2))
         else
            call wfat_print_wfat0_ints(nch, nqn, natm, lmax, 
     &           int_mb(k_npa), int_mb(k_mpa), dcpl_mb(k_in))
         endif


c        =============================================
c        ==== Calculate and print some quantities ====
c        =============================================
         if (td=='td') then
            if (pr_srfc .or. pr_yild .or. pr_avbeta .or. pr_avgamma)
     &           call wfat_calc_td(rtdb, geom_hd, pr_srfc, pr_yild, 
     &           pr_avbeta, pr_avgamma, nlas, las_env, dbl_mb(k_al), 
     &           dbl_mb(k_wl), dbl_mb(k_pl), dbl_mb(k_tl), dbl_mb(k_dl),
     &           ntime, tmax, extf, nch, int_mb(k_or0), int_mb(k_sp),
     &           dbl_mb(k_eg), dbl_mb(k_cmm), nqn, int_mb(k_npa), 
     &           int_mb(k_mpa), lmax, zion, dcpl_mb(k_in), wfat1, 
     &           ndgsub, max_ndgsub, int_mb(k_dgs), int_mb(k_ndgs),
     &           dbl_mb(k_dps), dbl_mb(k_pbs), log_mb(k_w1x), 
     &           dcpl_mb(k_inj), dcpl_mb(k_ink1), dcpl_mb(k_ink2))
         else
            if (pr_srfc .or. pr_rate)
     &           call wfat_calc_ti(rtdb, geom_hd, pr_srfc, pr_rate, 
     &           fieldmag, nch, int_mb(k_or0), int_mb(k_sp), 
     &           dbl_mb(k_eg), dbl_mb(k_cmm), nqn, int_mb(k_npa),
     &           int_mb(k_mpa), lmax, zion, dcpl_mb(k_in), wfat1, 
     &           ndgsub, max_ndgsub, int_mb(k_dgs), int_mb(k_ndgs),
     &           dbl_mb(k_dps), dbl_mb(k_pbs), log_mb(k_w1x), 
     &           dcpl_mb(k_inj), dcpl_mb(k_ink1), dcpl_mb(k_ink2))

            if (pr_tmd .or. pr_tmd_atom) 
     &           call wfat_calcprint_tmd(rtdb, geom_hd, pr_tmd, 
     &           pr_tmd_atom, fieldmag, nch, int_mb(k_or0), 
     &           int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_cmm), nqn, 
     &           int_mb(k_npa), int_mb(k_mpa), lmax, dcpl_mb(k_in),
     &           res_mode)
         endif
         
         
         if (.not. MA_free_heap(l_in))
     &        call errquit('WFAT: intgl deallocation failed', 3, MA_ERR)
         call wfat_erase_wfat1_ti(l_inj, l_ink1, l_ink2)
         rtdb_stat = rtdb_parallel(.true.)

      else if (lmax < 0) then
         call wfat_calc_ti_nopw(rtdb, geom_hd, bas_hd, pr_srfc, pr_rate,
     &        fieldmag, nbas0, nspin0, dbl_mb(k_rd), nch, int_mb(k_or0),
     &        dbl_mb(k_cf), int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_cmm),
     &        dbl_mb(k_ps), nqn, int_mb(k_npa), int_mb(k_mpa), zion,
     &        gqual, wf_msm, smh, ismask_e, ergmask, amsk, bmsk,
     &        ismask_i, mask_param, wfat1, ndgsub, max_ndgsub,
     &        int_mb(k_dgs), int_mb(k_ndgs), dbl_mb(k_dps),
     &        dbl_mb(k_pbs), log_mb(k_w1x))
      endif
      if (wfat1)
     &     call wfat_erase_wfat1_dgn(l_dgs, l_dgss, l_ndgs, l_dps,
     &     l_pbs)
      call ga_sync()
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_nd))
     &     call errquit('WFAT: numden deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_rd))
     &     call errquit('WFAT: rdm deallocation failed', 61, MA_ERR)
      call wfat_erase_cf(l_cf)
      call wfat_erase_channels(l_or0, l_sp, l_or)
      call wfat_erase_paraqn(l_npa, l_mpa)
      if (td=='td') 
     &     call wfat_erase_tdfield(l_al, l_wl, l_pl, l_tl, l_dl)
      call wfat_erase_channelprop(l_eg, l_ps, l_cmm, l_idp)
      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 1, GA_ERR)
      call wfat_erase_wfat1(l_w1x)
      call wfat_erase_wfat1_intgl(g_j, g_k, l_pb)
      if (wfat1) then
         if (.not. MA_free_heap(l_dpo))
     &        call errquit('WFAT: dpmo deallocation failed', 3, MA_ERR)
      endif

      call wfat_terminate(rtdb)


cforCHECK      if (.not. rtdb_print(rtdb, .false.))
cforCHECK     &     call errquit('WFAT: Cannot print RTDB.', 1, RTDB_ERR)

      
      task_oewfat = .true.
      call wfat_close_task('OE')
      
      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_calcprint_tmd(rtdb, geom, pr_tmd, pr_tmd_atom, ef,
     &           nch, orbid0, spinid, erg, chmmt, nqn, nid, mid, lmax,
     &           intgl, rmode)
      
      implicit none
#include "global.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
ccOLD #include "wfat_qm_param.fh"            
      integer rtdb
      integer geom
      logical pr_tmd
      logical pr_tmd_atom
      double precision ef
      integer nch
      integer orbid0(nch)
      integer spinid(nch)
      double precision erg(nch)
      double precision chmmt(3, nch)
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      integer lmax
      double complex intgl(1:*)
      logical rmode

      logical gap
      integer ix, iy, ndid, nat
      integer nsf, ntmd, nktmd, ntms
      integer l_sf, l_tm, l_tms
      integer k_sf, k_tm, k_tms
      double precision kx, ky, kmax, dk, bt_tmd, gm_tmd


      ndid = ga_nodeid()
      if (ndid /= 0) return

      
c     ==== Obtaining beta and gamma ====
      if (.not. rtdb_get(rtdb, 'wfat:beta_tmd', mt_dbl, 1, bt_tmd))
     &     call errquit('WFAT: Failed reading bt_tmd from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:gamma_tmd', mt_dbl, 1, gm_tmd))
     &     call errquit('WFAT: Failed reading gm_tmd from rtdb',
     &     1, INPUT_ERR)

      
c     ==== Allocation of structure factors and TMD ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 77, GEOM_ERR)
      nsf = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfc', l_sf, k_sf))
     &     call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
      call wfat_get_srfc0(geom, bt_tmd, gm_tmd, nch, erg,
     &                    nqn, mid, lmax, chmmt, intgl, .false.,
     &                    dcpl_mb(k_sf))
      ntmd = nat * nch
      if (.not. MA_alloc_get(MT_DCPL, ntmd, 'tmd', l_tm, k_tm))
     &     call errquit('WFAT: tmd allocation failed', 2, MA_ERR)
      ntms = nat * nch
      if (.not. MA_alloc_get(MT_DBL, ntmd, 'tmd_ksum', l_tms, k_tms))
     &     call errquit('WFAT: tmd_ksum allocation failed', 2, MA_ERR)

      
c     ==== Calculating max. momentum to calculate the TMD at ====
      call wfat_get_kmax(ef, nch, erg, nqn, nid, mid, kmax)
      if (.not. rtdb_get(rtdb, 'wfat:nktmd', mt_int, 1, nktmd))
     &     call errquit('WFAT: Failed reading nktmd from rtdb',
     &     1, RTDB_ERR)
      dk = 2 * kmax / (nktmd - 1)
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Max. transverse momentum = ", e14.6, '//
     &      '" rad/Bohr")') kmax


c     ==== Initializing the printings ====      
      if (ndid==0 .and. pr_tmd) call wfat_initprint_tmd(nch, orbid0,
     &                                  spinid)
      if (ndid==0 .and. pr_tmd_atom) call wfat_initprint_tmdar(nch, 
     &                                       nat, orbid0, spinid)


c     ==== Calculate the TMDs ====
      dbl_mb(k_tms:k_tms+ntms-1) = 0.0d0
      ky = -kmax
      do iy = 1, nktmd
         kx = -kmax
         do ix = 1, nktmd

            if (ix /= nktmd) then
               gap = .false.
            else
               gap = .true.
            endif
            call wfat_calc_tmd(kx, ky, nch, erg, nqn, nid, mid,
     &                         ef, dcpl_mb(k_sf), dcpl_mb(k_tm))
            
            if (ndid == 0) then
               if (pr_tmd) call wfat_print_tmd(kx, ky, nch,
     &                          dcpl_mb(k_tm), gap)
               if (pr_tmd_atom) call wfat_print_tmdar(kx, ky, nch,
     &                               nat, dcpl_mb(k_tm), gap)
            endif

            if (pr_tmd_atom .and. rmode)
     &           dbl_mb(k_tms:k_tms+ntms-1) = dbl_mb(k_tms:k_tms+ntms-1)
     &                              + abs(dcpl_mb(k_tm:k_tm+ntmd-1))

            kx = kx + dk
         enddo
         ky = ky + dk
      enddo


      if (pr_tmd_atom .and. rmode .and. ndid==0) then
         dbl_mb(k_tms:k_tms+ntms-1) = dbl_mb(k_tms:k_tms+ntms-1)*dk**2
         call wfat_print_tms(geom, nch, nat, dbl_mb(k_tms))
      endif

      
c     ==== Terminating printings ====
      if (ndid==0 .and. pr_tmd) call wfat_endprint_tmd()
      if (ndid==0 .and. pr_tmd_atom) call wfat_endprint_tmdar()
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sf))
     &     call errquit('WFAT: srfc deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_tm))
     &     call errquit('WFAT: tmd deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_tms))
     &     call errquit('WFAT: tmd_ksum deallocation failed',
     &     3, MA_ERR)



      
      end
c===============================================================
