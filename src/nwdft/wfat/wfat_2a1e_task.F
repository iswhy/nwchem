c===============================================================
c===============================================================
      logical function task_wfat_2a1e(rtdb)

      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "wfat_qm_param.fh"
#include "wfat.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "wfat_field_param.fh"
#include "geom.fh"
      integer rtdb

      character(len=FLD_TP_LEN) las_env(max_nfield)
      character*2 td
      character*3 pottype
      character*7 jobtype
      character*16 tag
      logical pr_srfc, pr_yild, pr_rate, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom, os, oc, or, int_nosave, l1, extf
      logical res_mode, wfat1, rtdb_stat, fedvr, ismask_i, ismask_e
      integer l_cf, l_eg, l_npa, l_mpa, l_in, l_inj, l_ink1
      integer l_ink2, l_pb, l_cfx, l_cfy, l_cfz, l_ndvr, l_fen, l_xquad
      integer l_wquad, l_al, l_wl, l_pl, l_tl, l_dl, l_w1x
      integer k_cf, k_eg, k_npa, k_mpa, k_in, k_inj, k_ink1
      integer k_ink2, k_pb, k_cfx, k_cfy, k_cfz, k_ndvr, k_fen, k_xquad
      integer k_wquad, k_al, k_wl, k_pl, k_tl, k_dl, k_w1x
      integer i, mch, sphi, spinid, orbid, nr, nt, maxiter, i_eg, g_cf
      integer lmax, gqual, nqn, nsh, nintg, g_il, ndid, nlgl, nel_2a
      integer nintj, nintk, g_j, g_k(2), npb, nr_m, nt_m, nr_p, nt_p
      integer ndgsub, max_ndgsub, dgs(1,1), ndgs(1), nfe, maxn, idum
      integer nlas, ntime, nw1x
      double precision init_ke, out_ke, mult_ke, edip(3), chpos(3)
      double precision chmmt(3), fieldmag, pb(3,3), polb(3), polb33(3,3)
      double precision dps(1,1,3,1), pbs(1,1,6,1), ddum, potparam(2)
      double precision mask_param(2), ergmask, amsk, bmsk, iondip(3)
      double precision tmax

      integer nch
      parameter(nch=1)


      ndid = ga_nodeid()
      task_wfat_2a1e = .false.
      

c     ==============================
c     ==== Initiate WFAT module ====
c     ==============================
      call wfat_2a1e_init(rtdb)


      res_mode = .false.
      if (.not. rtdb_get(rtdb, 'wfat:fedvr', mt_log, 1, fedvr)) then
         fedvr = .false.
         if (ndid == 0) then
            write(luout, '(/)')
            write(luout, '("WARNING: Proceeding without using FEDVR '//
     &           'radial basis. Laguerre DVR radial basis will be '//
     &           'used. Be warned that Laguerre DVR basis is not '//
     &           'thoroughly tested and is unstable for the latter '//
     &           'WFAT integral computations.")')
         endif
      endif

      
c     ==== Task title ====
      jobtype = 'OE-2a1e'
      call wfat_prep_tasktitle('OE-2a1e')


c     ==== Extract the general parameters ====
      call wfat_prep_general(rtdb, gqual, lmax)

      
c     ==== Extract the parameters of the quantum system ====
      call wfat_2a1e_prep_system(rtdb, mch, sphi, spinid, orbid,
     &     pottype, nel_2a, potparam)
      if (ndid == 0) call wfat_2a1e_stateinfo(mch, sphi, spinid, orbid)
      call wfat_2a1e_prep_bas(rtdb, fedvr, init_ke, mult_ke, maxiter, 
     &     nr, nt, nfe, maxn, l_ndvr, k_ndvr, l_fen, k_fen, l_xquad, 
     &     k_xquad, l_wquad, k_wquad)
            

c     ==== Allocate the coefficients and energies array ====
      call wfat_2a1e_prep_cf(nr, nt, g_cf, l_cf, k_cf, l_eg, k_eg)
      i_eg = k_eg + (orbid-1)
      call wfat_manerg(rtdb, 1, dbl_mb(i_eg))


c     ==== Obtain parabolic quantum numbers ====
      call wfat_prep_paraqn(rtdb, nqn, l_npa, k_npa, l_mpa, k_mpa)
      

c     ==== Obtain WFAT(1) information ====
      call wfat_prep_wfat1(rtdb, nch, nqn, wfat1, nlgl, nw1x, l_w1x,
     &     k_w1x)
      

c     ==== Obtain print switches ====
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


c     ==== Print DVR grids and weights ====
      if (ndid == 0) then
         if (fedvr) then
            call wfat_2a1efd_print_rgrid(nfe, int_mb(k_ndvr), maxn, 
     &           dbl_mb(k_xquad), dbl_mb(k_wquad))
         else
            call wfat_2a1e_print_rgrid(mch, nr)
         endif
         call wfat_2a1e_print_tgrid(mch, nt)
      endif
c      call errquit('stopping here', 4, UERR)
      
      
c     ==== Solve the Hamiltonian for the target ====
c     ====   eigenstate (the ionized channel)   ====
      if (fedvr) then
         call wfat_2a1efd_target(geom_hd, mch, sphi, orbid, pottype, 
     &        nel_2a, potparam, nfe, int_mb(k_ndvr), maxn, 
     &        dbl_mb(k_xquad), dbl_mb(k_wquad), nr, nt, g_cf,
     &        dbl_mb(k_eg), edip)
      else
         call wfat_2a1e_target(geom_hd, mch, sphi, orbid, pottype, 
     &        nel_2a, potparam, nr, nt, init_ke, mult_ke, maxiter, g_cf,
     &        dbl_mb(k_eg), edip, out_ke)
      endif
      call wfat_2a1e_mandip(rtdb, edip)
      call ga_get(g_cf, 1, nr*nt, orbid, orbid, dbl_mb(k_cf), nr*nt)
      
      call wfat_2a1e_prep_cf1(wfat1, fedvr, rtdb, mch, nr, nt, nr_m, 
     &     nt_m, nr_p, nt_p, l_cfx, k_cfx, l_cfy, k_cfy, l_cfz, k_cfz)
      if (wfat1) then
         if (fedvr) then
            call wfat_2a1efd_perturb(geom_hd, mch, sphi, orbid, nfe, 
     &           int_mb(k_ndvr), maxn, dbl_mb(k_xquad), dbl_mb(k_wquad),
     &           nr, nt, nt_m, nt_p, pottype, nel_2a, potparam, g_cf, 
     &           dbl_mb(k_eg), dbl_mb(k_cfx), dbl_mb(k_cfy), 
     &           dbl_mb(k_cfz), polb)
         else
            call wfat_2a1e_perturb(geom_hd, mch, sphi, orbid, pottype, 
     &           nel_2a, potparam, out_ke, mult_ke, nr, nt, nr_m, nt_m, 
     &           nr_p, nt_p, g_cf, dbl_mb(k_eg), dbl_mb(k_cfx), 
     &           dbl_mb(k_cfy), dbl_mb(k_cfz), polb)
         endif
         if (ndid == 0)
     &        call wfat_print_w1x(nch, orbid, spinid, nqn,
     &        int_mb(k_npa), int_mb(k_mpa), log_mb(k_w1x))
      endif
      polb33 = 0.0d0
      do i = 1, 3
         polb33(i,i) = polb(i)
      enddo


c     ==== Obtain properties of the ionized channel ====
      call wfat_2a1e_orbcent(rtdb, geom_hd, edip, chpos, chmmt)
      iondip = (/ 0.0d0, 0.0d0, 0.0d0 /)
      if (ndid == 0) then
         call wfat_print_orb_property(nch, orbid, 1, dbl_mb(i_eg),
     &        chpos, chmmt, iondip)
         if (wfat1) call wfat_print_polb_tensors(nch, polb33)
      endif
      call ga_sync()
c     edip = Unshifted dipole moment.
c     chmmt = Shifted dipole moment.

      
c     ==== Get mask parameters ====
      call wfat_get_masks(rtdb, dbl_mb(k_eg), ismask_i, mask_param, 
     &     ismask_e, ergmask, amsk, bmsk)
      if (ismask_i .and. ndid==0) call wfat_print_mask(mask_param(1),
     &     mask_param(2))
      

c     ==== Calculate (or load) WFAT integrals ====
      nsh = (1 + lmax)**2
      nintg = nsh * natm * nqn * nch
      if (.not. ga_create(mt_dcpl, nintg, 1, 'intgl', 0, 0, g_il))
     &     call errquit('WFAT: Cannot create intgl global array.',
     &     1, GA_ERR)
      call ga_zero(g_il)
      

c     ==== Prepare WFAT(1) integrals ====
      call wfat_prep_wfat1_intgl(wfat1, lmax, natm, nqn, nch, nintj, 
     &     g_j, nintk, g_k, npb, l_pb, k_pb)
      if (wfat1)
     &     call wfat_2a1e_assign_polb(polb, dbl_mb(k_pb))   ! k_pb is not needed apparently, check if this line can be removed.
      

      call util_get_rtdb_state(rtdb, os, oc, or)
      if (lmax >= 0) then
         if (os .and. .not. or) then
c           ==== Compute the OEWFAT integrals ====
            call wfat_2a1e_intgl(rtdb, geom_hd, nqn, int_mb(k_npa), 
     &           int_mb(k_mpa), lmax, gqual, zion, wf_msm, sphi, mch, 
     &           nr, nt, out_ke, dbl_mb(k_cf), dbl_mb(i_eg), chmmt, 
     &           chpos, pottype, nel_2a, potparam, wfat1, nlgl, nr_m, 
     &           nt_m, nr_p, nt_p, dbl_mb(k_cfx), dbl_mb(k_cfy), 
     &           dbl_mb(k_cfz), log_mb(k_w1x), fedvr, nfe, 
     &           int_mb(k_ndvr), maxn, dbl_mb(k_xquad), dbl_mb(k_wquad),
     &           ismask_i, mask_param, g_il, g_j, g_k)

c           ==== Store the OE-WFAT integrals freshly calculated ====
            if (.not. rtdb_get(rtdb, 'wfat:intgl_nosave', mt_log, 1,
     &           int_nosave)) call errquit('WFAT: Cannot read '//
     &           'intgl_nosave from RTDB', 1, RTDB_ERR)
            if (ndid == 0 .and. .not. int_nosave) then
               call wfat_store_intgl(rtdb, nsh, natm, nqn, nch,
     &              dbl_mb(i_eg), chpos, chmmt, iondip, g_il)
            else if (ndid == 0 .and. int_nosave) then
               write(luout, '("ATTENTION: WFAT integrals are not '//
     &              'saved.")')
            endif
            call ga_sync()
         else if (.not. os .and. or) then
c           ==== Binary-load WFAT integrals on 'restart' startup mode ====
            if (ndid == 0) then
               call wfat_load_intgl(rtdb, nsh, natm, nqn, nch,
     &              dbl_mb(i_eg), chpos, chmmt, iondip, g_il)
clater               call wfat_print_orb_property(nch, int_mb(k_or), 
clater     &              int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_ps),
clater     &              dbl_mb(k_cmm), dbl_mb(k_idp))
            endif
            call ga_sync()
            
            call ga_brdcst(msg_wfat_eg, dbl_mb(i_eg),
     &                     ma_sizeof(MT_DBL, nch, MT_BYTE), 0)
            call ga_brdcst(msg_wfat_cmm, chmmt,
     &                     ma_sizeof(MT_DBL, 3, MT_BYTE), 0)
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
            
      
      
cOLD     ==== No partial wave mode ====
cOLD      if (.not. fedvr)
cOLD     &     call wfat_2a1e_task_nopw(rtdb, geom_hd, nqn, int_mb(k_npa),
cOLD     &     int_mb(k_mpa), gqual, zion, wf_msm, sphi, mch, nr, nt, orbid,
cOLD     &     dbl_mb(k_cf), out_ke, dbl_mb(i_eg), edip, pottype, nel_2a, 
cOLD     &     potparam, wfat1, nlgl, nr_m, nt_m, nr_p, nt_p, dbl_mb(k_cfx),
cOLD     &     dbl_mb(k_cfy), dbl_mb(k_cfz))
      

      if (wfat1) then
clater            ctr = 0.0d0
clater            call wfat_calc_dpmo(ctr, dbl_mb(k_dpo))
         ndgsub = 1             ! 1)
         max_ndgsub = 1
         dgs(1,1) = 1
         ndgs(1) = 1
         dps(1,1,:,1) = chmmt
         pbs(1,1,:,1) = 0.0d0
         pbs(1,1,1:3,1) = polb
c        1) The dipole moment matrix in the eigenstate basis of linear
c           molecules is always diagonal. Hence, any degenerate eigenstate
c           maybe treated in the same way as nondegenerate ones in the
c           perturbation analysis.
      endif

      
c     ==== Ionization calculations and printings ====
      if (lmax >= 0 .and. ndid == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &        call errquit('WFAT: Cannot allocate intgl.', 41, MA_ERR)
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

c        ==================================================
c        ==== Calculate and print important quantities ====
c        ==================================================
         if (td == 'td') then
            if (pr_srfc .or. pr_yild .or. pr_avbeta .or. pr_avgamma)
     &           call wfat_calc_td(rtdb, geom_hd, pr_srfc, pr_yild, 
     &           pr_avbeta, pr_avgamma, nlas, las_env, dbl_mb(k_al), 
     &           dbl_mb(k_wl), dbl_mb(k_pl), dbl_mb(k_tl), dbl_mb(k_dl),
     &           ntime, tmax, extf, nch, orbid, spinid, dbl_mb(i_eg), 
     &           chmmt, nqn, int_mb(k_npa), int_mb(k_mpa), lmax, zion, 
     &           dcpl_mb(k_in), wfat1, ndgsub, max_ndgsub, dgs, ndgs,
     &           dps, pbs, log_mb(k_w1x), dcpl_mb(k_inj), 
     &           dcpl_mb(k_ink1), dcpl_mb(k_ink2))
         else
            if (pr_srfc .or. pr_rate)
     &           call wfat_calc_ti(rtdb, geom_hd, pr_srfc, pr_rate,
     &           fieldmag, nch, orbid, spinid, dbl_mb(i_eg), chmmt, 
     &           nqn, int_mb(k_npa), int_mb(k_mpa), lmax, zion, 
     &           dcpl_mb(k_in), wfat1, ndgsub, max_ndgsub, dgs, ndgs,
     &           dps, pbs, log_mb(k_w1x), dcpl_mb(k_inj), 
     &           dcpl_mb(k_ink1), dcpl_mb(k_ink2))

clater            if (pr_tmd .or. pr_tmd_atom) 
clater     &           call wfat_calcprint_tmd(rtdb, geom_hd, pr_tmd, 
clater     &           pr_tmd_atom, fieldmag, nch, int_mb(k_or0), 
clater     &           int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_cmm), nqn, 
clater     &           int_mb(k_npa), int_mb(k_mpa), lmax, dcpl_mb(k_in),
clater     &           res_mode)
         endif
         
         if (.not. MA_free_heap(l_in))
     &        call errquit('WFAT: intgl deallocation failed', 3, MA_ERR)
         call wfat_erase_wfat1_ti(l_inj, l_ink1, l_ink2)
         rtdb_stat = rtdb_parallel(.true.)
      else if (lmax < 0) then
         call wfat_2a1e_ti_nopw(rtdb, geom_hd, nqn, int_mb(k_npa),
     &        int_mb(k_mpa), gqual, zion, wf_msm, sphi, mch, nr, nt, 
     &        out_ke, dbl_mb(k_cf), dbl_mb(i_eg), chmmt, chpos, orbid,
     &        spinid, pottype, nel_2a, potparam, pr_srfc, pr_rate, 
     &        fieldmag, wfat1, nr_m, nt_m, nr_p, nt_p, dbl_mb(k_cfx), 
     &        dbl_mb(k_cfy), dbl_mb(k_cfz), pbs, log_mb(k_w1x), fedvr, 
     &        nfe, int_mb(k_ndvr), maxn, dbl_mb(k_xquad), 
     &        dbl_mb(k_wquad), ismask_i, mask_param)
      endif
      call ga_sync()
      
      
c     ==== Deallocations ====
      call wfat_erase_paraqn(l_npa, l_mpa)
      call wfat_2a1e_erase_bas(l_ndvr, l_fen, l_xquad, l_wquad)
      call wfat_2a1e_erase_cf(g_cf, l_cf, l_eg)
      call wfat_2a1e_erase_cf1(l_cfx, l_cfy, l_cfz)
      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 45, GA_ERR)
      call wfat_erase_wfat1(l_w1x)
      call wfat_erase_wfat1_intgl(g_j, g_k, l_pb)
      if (td=='td') 
     &     call wfat_erase_tdfield(l_al, l_wl, l_pl, l_tl, l_dl)


      task_wfat_2a1e = .true.
      call wfat_close_task('OE-2a1e')
      
      end
c===============================================================
      
