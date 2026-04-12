c=======================================================
c=======================================================
      subroutine wfat_calc_ti_nopw(rtdb, geom, bas, pr_srfc, pr_rate,
     &     pr_tmd, nb, ns, rdm, nch, orbid0, vch, spinid, erg, 
     &     chmmt, chpos, nqn, nid, mid, zi, gq, celm, smh, ismask_e, 
     &     ergmask, amsk, bmsk, ismask, mask_param, wfat1, nsub, maxn, 
     &     dgs, ndgs, dps, pbs, wfat1_exc)

      implicit none
#include "util.fh"
#include "global.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb                    ! [input]
      integer geom              ! [input]
      integer bas
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_rate                 ! [input], if .true. then print ion yield.
      logical pr_tmd                  ! [input], if .true. then print TMD.
      integer nb
      integer ns
      double precision rdm(nb,nb,ns)
      integer nch               ! [input], number of ionized channels.
      integer orbid0(nch)       ! [input], orbital ID of ionized channels.
      double precision vch(nb,*)      ! [input], the minimum number of columns is nch if wfat1 is .false., otherwise 4*nch.
      integer spinid(nch)             ! [input], spin value (1 or 2) of ionized channels.
      double precision erg(nch)       ! [input]
      double precision chmmt(3,nch)   ! [input]
      double precision chpos(3,nch) ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      double precision zi       ! [input]
      integer gq
      integer celm
      double precision smh
      logical ismask_e
      double precision ergmask        ! [input]. Only used when ismask_e=.true.. Must always be supplied regardless.
      double precision amsk           ! [input]. Only used when ismask_e=.true.. Must always be supplied regardless.
      double precision bmsk           ! [input]. Only used when ismask_e=.true.. Must always be supplied regardless.
      logical ismask                  ! [input],
      double precision mask_param(2)  ! [input]. Only used when ismask=.true.. Must always be supplied regardless.
      logical wfat1                   ! [input]
      integer nsub                    ! [input]. Only used when wfat1=.true..
      integer maxn                    ! [input]. Only used when wfat1=.true..
      integer dgs(1:*)                ! [input]. Only used when wfat1=.true..
      integer ndgs(1:*)               ! [input]. Only used when wfat1=.true..
      double precision dps(1:*)       ! [input]. Only used when wfat1=.true..
      double precision pbs(1:*)       ! [input]. Only used when wfat1=.true..
      logical wfat1_exc(1:*)          ! [input]. Only used when wfat1=.true..

      character*12 dt_str
      character*16 tagt
      logical gap, interm, cgrid(nch), rtdb_stat
      integer l_sg, l_fc, l_rt, l_icp, l_ipp, l_pt, l_wt, l_vc, l_cl
      integer l_ma, l_tu0, l_tu1, l_acg, l_acj, l_ack, l_dum, l_sc, l_sx
      integer l_pbd, l_sj, l_sk, l_on, l_si, l_pos, l_pto, l_wto, l_vco
      integer l_clo, l_mao, l_tu0o, l_tu1o
      integer k_sg, k_fc, k_rt, k_icp, k_ipp, k_pt, k_wt, k_vc, k_cl
      integer k_ma, k_tu0, k_tu1, k_acg, k_acj, k_ack, k_dum, k_sc, k_sx
      integer k_pbd, k_sj, k_sk, k_on, k_si, k_pos, k_pto, k_wto, k_vco
      integer k_clo, k_mao, k_tu0o, k_tu1o
      integer natom, n_sg, nfc, nrt, n_pt, n_wt, n_pto, n_wto, n_acg
      integer n_acj, n_ack, n_dum, nsc, nsx, npbd, n_sj, n_sk, n_on
      integer nch1, ndid, nbt, ngm, natp, natop, nptp, nptop, i, k, iqn
      integer i_on, ig, ib, i_tu0, i_tu0o, i_vc, i_vco, i_acj, icpo
      integer ippo(2), ozn, geomo, pr_chk(20), nang, kk
      double precision bt_start, bt_last, dbt, dbtr, gm_start, gm_last
      double precision dgm, dgmr, dthr, bt, gm, rate, vdum(3), fieldm
      double precision deg2rad, t_now, t_last, t_start, dt_all, dt_prog
      parameter(deg2rad=0.01745329251994d0)


      call ga_sync()
      ndid = ga_nodeid()
      
      
c     ==== Allocation of structure factors and ion yields ====
      if (.not. geom_ncent(geom, natom))
     &     call errquit('WFAT: Getting natom failed', 31, GEOM_ERR)
      call wfat_prep_rate(natom, nqn, nch, n_sg, l_sg, k_sg, nfc,
     &     l_fc, k_fc, nrt, l_rt, k_rt)
      
c     ==== Obtaining info on beta ====
      call wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)
      dbtr = dbt * deg2rad

c     ==== Obtaining info on gamma ====
      call wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)
      dgmr = dgm * deg2rad

c     ==== Obtaining laser field ====
      call wfat_prep_dcfield(rtdb, fieldm)
      
c     ==== Calculate and gather nuclei grid points ====
      call wfat_alloc_gridrange(natom, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natp, nptp)
      call wfat_prep_allgrid(nptp, n_pt, l_pt, k_pt, n_wt, l_wt, k_wt)
      call wfat_collect_mfgrid(geom, gq, natp, nptp, int_mb(k_icp),
     &     int_mb(k_ipp), dbl_mb(k_pt), dbl_mb(k_wt))

c     ==== Calculate and gather channel origin grid points ====
      if (.not. geom_create(geomo, 'geomo'))
     &     call errquit('WFAT: Error in creating geomo.', 25, GEOM_ERR)
      vdum = (/0.0d0, 0.0d0, 0.0d0/)     ! 1)
      tagt = 'dummycent'
      if (.not. rtdb_get(rtdb, 'wfat:ozn', mt_int, 1, ozn))
     &     call errquit('WFAT: Cannot read ozn.', 40, RTDB_ERR)
      if (.not. geom_cart_set(geomo, 1, tagt, vdum, dble(ozn)))
     &     call errquit('WFAT: Error in initializing geometry object.',
     &     32, GEOM_ERR)
c     1) The coordinate of the only center in geomo is set to the zero vector
c     because we will use wfat_collect_mfgrid below to generate the grids
c     in the frame of channel origin. This choice of the coordinate of the
c     center in geomo will not and should not affect other parts of the code.
      call wfat_get_gridrange(geomo, ndid, gq, icpo, ippo, natop, nptop)    ! 3)
c     3) Since geomo contains only one center, we have icpo = natop = 1.
      call wfat_prep_allgrid(nptop, n_pto, l_pto, k_pto, n_wto, l_wto,
     &     k_wto)
      call wfat_collect_mfgrid(geomo, gq, natop, nptop, icpo, ippo, 
     &     dbl_mb(k_pto), dbl_mb(k_wto))    ! 2)
c     2) Although this subroutine returns grids in molecular frame, since
c     the coordinate of the only center in geomo is the zero vector (see
c     note 1)), the output is effectively the grids in the frame of the
c     channel center. So, k_pto contains grid coordinates of the channel center
c     relative to the channel center itself, that's why they are independent of
c     ionizd orbitals.
      if (.not. geom_destroy(geomo))
     &     call errquit('WFAT: Cannot destroy geomo.', 52, GEOM_ERR)
      call wfat_need_cgrid(rtdb, geom, nch, chpos, cgrid)
      call wfat_print_gridpt(rtdb, geom, gq, ozn, nch, cgrid, chpos)
      

c     ==== Allocate new spinid and chpos that can ====
c     ====  take care of wfat1 when it is .true.  ====
      nch1 = 0
      if (wfat1) nch1 = 3 * nch
      if (.not. MA_alloc_get(MT_INT, nch+nch1, 'spinid1', l_si, k_si))
     &     call errquit('WFAT: spinid1 allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, 3*(nch+nch1), 'chpos1',
     &     l_pos, k_pos))
     &     call errquit('WFAT: chpos1 allocation failed', 59, MA_ERR)
      call wfat_get_wfat1_sppos(wfat1, nch, spinid, chpos, int_mb(k_si),
     &     dbl_mb(k_pos))
c     ^^ At this step, only k_si is needed for wfat_collect_vfpsi below.
c     k_pos obtained above is not needed.

c     ==== Calculate and gather V_eff|psi> at nuclear grid points ====
      if (.not. MA_alloc_get(MT_DBL, nptp*(nch+nch1), 'vf_psi',
     &     l_vf, k_vf))
     &     call errquit('WFAT: vf_psi allocation failed', 81, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp*(nch+nch1), 'vs_psi',
     &     l_vs, k_vs))
     &     call errquit('WFAT: vs_psi allocation failed', 81, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptop*(nch+nch1), 'vs_psi_o',
     &     l_vso, k_vso))
     &     call errquit('WFAT: vs_psi_o allocation failed', 81, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp*(nch+nch1), 'vc_psi',
     &     l_vc, k_vc))
     &     call errquit('WFAT: vc_psi allocation failed', 81, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptop*(nch+nch1), 'vc_psi_o',
     &     l_vco, k_vco))
     &     call errquit('WFAT: vc_psi_o allocation failed', 81, MA_ERR)
      if (ismask_e) call errquit('WFAT: ismask_e must not be true '//
     &     'when not using partial wave WFAT.', 60, UERR)
      call wfat_collect_vfpsi(geom, bas, nb, ns, nch, wfat1, vch,
     &     int_mb(k_si), nptp, dbl_mb(k_pt), rdm, smh, .false., ergmask,
     &     amsk, bmsk, dbl_mb(k_vf))
      
c     ==== Calculate and gather the mask function  ====
c     ====         at nuclear grid points          ====
      if (ismask) then
         if (.not. MA_alloc_get(MT_DBL, nptp, 'maskf', l_ma, k_ma))
     &        call errquit('WFAT: maskf allocation failed', 81, MA_ERR)
         if (.not. MA_alloc_get(MT_DBL, nptop*nch, 'maskfo',
     &        l_mao, k_mao))
     &        call errquit('WFAT: maskfo allocation failed', 81, MA_ERR)
         call wfat_collect_maskf(geom, natom, nch, nptp, dbl_mb(k_pt), 
     &        mask_param, dbl_mb(k_ma))
      endif

c     ==== Allocate tau at all grid points and as. coeff. ====
      if (.not. MA_alloc_get(MT_DCPL, nptp*nqn*nch, 'tau0',
     &     l_tu0, k_tu0))
     &     call errquit('WFAT: tau0 allocation failed', 81, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nptop*nqn*nch, 'tau0o',
     &     l_tu0o, k_tu0o))
     &     call errquit('WFAT: tau0o allocation failed', 81, MA_ERR)
      if (wfat1) then
         if (.not. MA_alloc_get(MT_DCPL, nptp*nqn*nch, 'tau1',
     &        l_tu1, k_tu1))
     &        call errquit('WFAT: tau1 allocation failed', 81, MA_ERR)
         if (.not. MA_alloc_get(MT_DCPL, nptop*nqn*nch, 'tau1o',
     &        l_tu1o, k_tu1o))
     &        call errquit('WFAT: tau1o allocation failed', 81, MA_ERR)
      endif
      call wfat_prep_acoef(natom, nqn, nch, wfat1, n_acg, l_acg,
     &     k_acg, n_acj, l_acj, k_acj, n_ack, l_ack, k_ack, n_dum,
     &     l_dum, k_dum)
c     k_acg is a 3D array with the size of natom*nqn*nch.
c     k_acj is a 3D array with the size of natom*nqn*3*nch.
c     k_ack is a 3D array with the size of natom*nqn*nch.

c     ==== Allocate cell functions at all grid points ====
      if (.not. MA_alloc_get(MT_DBL, nptp*nch, 'celf', l_cl, k_cl))
     &     call errquit('WFAT: celf allocation failed', 81, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptop*nch, 'celfo',
     &     l_clo, k_clo))
     &     call errquit('WFAT: celfo allocation failed', 81, MA_ERR)

c     ==== Obtain field factor ====
      call wfat_get_field_factor(rtdb, nch, nqn, nid, mid, zi, fieldm,
     &     erg, dbl_mb(k_fc))
      
c     ==== When WFAT1 is active ====
      if (wfat1) then
         call wfat_prepare_an(maxn, natom, nch, nqn, nsc, l_sc, k_sc, 
     &        nsx, l_sx, k_sx, npbd, l_pbd, k_pbd, n_sj, l_sj, k_sj,
     &        n_sk, l_sk, k_sk, n_on, l_on, k_on)
         if (.not. rtdb_get(rtdb, 'wfat:dthr', mt_dbl, 1, dthr))
     &        call errquit('WFAT: Failed reading dthr from rtdb',
     &        1, RTDB_ERR)
         do i = 1, nch
            do iqn = 1, nqn
               i_on = k_on + 2*nqn*(i-1) + 2*(iqn-1)
               call wfat_calc_o(nid(iqn), mid(iqn), erg(i), zi,
     &              dbl_mb(i_on))
            enddo
         enddo
         if (ndid == 0)
     &        call wfat_print_o(nch, nqn, nid, mid, dbl_mb(k_on))
         if (.not. rtdb_get(rtdb, 'wfat:wfat1_interm', mt_log, 1,
     &        interm)) interm = .false.
      endif


c     ==== Initializing the printings ====
      rtdb_stat = rtdb_parallel(.false.)
      if (ndid == 0)
     &     call wfat_print_field_factor(rtdb, fieldm, nch, nqn, erg,
     &     nid, mid, zi, dbl_mb(k_fc))
      rtdb_stat = rtdb_parallel(.true.)
      if (ndid == 0 .and. pr_srfc) then
         call wfat_initprint_structuref(nch, orbid0, spinid, nqn, nid,
     &        mid)
         if (wfat1) call wfat_initprint_srfc1(nch, orbid0, spinid, nqn,
     &        nid, mid)
      endif
      if (ndid == 0 .and. pr_rate)
     &     call wfat_initprint_rate(nch, orbid0, spinid, nqn, nid, mid)
      if (ndid == 0 .and. pr_tmd)
     &     call wfat_print_kmax(fieldm, nch, erg, nqn, nid, mid)


      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'WFAT computation without '//
     &        'partial wave expansion', 50, .true.)
         write(luout, '("")')
         write(luout, 901)
         write(luout, '("")')
      endif
 901  format(
     &     5x,
     &     'WFAT simulation using explicit integral formula (without ',
     &     'partial wave expansion) has started. Beware that this is ',
     &     /,
     &     5x,
     &     'often a slow computation (so be patient) and may require ',
     &     'large amount of memory. Use the partial wave expansion ',
     &     /,
     &     5x,
     &     'algorithm by setting a nonnegative integer for lmax if ',
     &     'you want fast and equally accurate results.')      

      
      call ga_sync()
      
c     ==== Scan over the orientation angles ====
      call wfat_get_intgl_chkpt(nbt*ngm, 20, pr_chk)
      if (ndid == 0) write(luout, '(5x, "The total number of '//
     &     'orientation angles is ", i9, ".")') nbt*ngm

      
      if (.not. MA_alloc_get(MT_DBL, maxn*(maxn+1+max(nb,nptp)), 
     &     'scratch0', l_sc0, k_sc0))
     &     call errquit('WFAT: scratch0 allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nb*(nch+nch1), 'vch', l_ch, k_ch))
     &     call errquit('WFAT: vch allocation failed', 59, MA_ERR)

      
      t_last = util_cpusec()
      t_start = t_last
      nang = 0
      kk = 1
      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)

         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)

            if (ib /= nbt) then
               gap = .false.
            else
               gap = .true.
            endif

c           ==== Update wfat(0) quantities depending on diagonalization ====
            chpos0 = chpos
            call dcopy(nb*nch, vch(:,:nch), 1, dbl_mb(k_ch), 1)
            call dcopy(nptp*nch, dbl_mb(k_vf), 1, dbl_mb(k_vc), 1)
c           ^^ vch and k_vf must be preserved because they are used to construct
c           new vectors that diagonalize H_int at every orientation.
            call wfat_diag_update_w0(bt, gm, zi, nch, chmmt, nsub, maxn,
     &           dgs, ndgs, dthr, dps, dbl_mb(k_sc0), nb, nptp, 
     &           dbl_mb(k_ch), dbl_mb(k_vc), chpos0)

c           ==== Update wfat(1) quantities depending on diagonalization ====
            if (wfat1) then
               i_ch = k_ch + nb*nch
               call dcopy(nb*3*nch, vch(:,nch+1:), 1, dbl_mb(i_ch), 1)
               i_vc = k_vc + nptp*nch
               call dcopy(nptp*3*nch, dbl_mb(k_vf+nptp*nch), 1, 
     &              dbl_mb(i_vc), 1)
               call wfat_diag_update_w1(bt, gm, zi, nch, nsub, maxn, 
     &              dgs, ndgs, dthr, dps, pbs, dbl_mb(k_sc0), nb, nptp,
     &              dbl_mb(i_ch), dbl_mb(i_vc), dbl_mb(k_pbd))
            endif

            call wfat_get_wfat1_sppos(wfat1, nch, spinid, chpos0, 
     &           int_mb(k_si), dbl_mb(k_pos))
            call wfat_need_cgrid(rtdb, geom, nch, chpos0, cgrid)
            
c           ==== Calculates V_eff|psi(0)> and V_eff|psi(1)> ====
c           ====        at channel center grid points       ====
            call wfat_collect_vfpsi_o(geom, bas, nb, ns, nch, wfat1, 
     &           dbl_mb(k_ch), int_mb(k_si), dbl_mb(k_pos), nptop, 
     &           dbl_mb(k_pto), cgrid, rdm, zi, smh, .false., ergmask,
     &           amsk, bmsk, dbl_mb(k_vco))
            
c           ==== Calculates V_s|psi(0)> and V_s|psi(1)> at all grid points ====
            call wfat_collect_vspsi(geom, bas, nb, nch, wfat1, 
     &           dbl_mb(k_ch), dbl_mb(k_pos), nptp, nptop, dbl_mb(k_pt),
     &           dbl_mb(k_pto), cgrid, zi, smh, .false., ergmask, amsk,
     &           bmsk, dbl_mb(k_vs), dbl_mb(k_vso))

c           ==== Calculates V_c|psi(0)> and V_c|psi(1)> at all grid points ====
            call daxpy(nptp*(nch+nch1), 1.0d0, dbl_mb(k_vs), 1,
     &           dbl_mb(k_vc), 1)
            call daxpy(nptop*(nch+nch1), 1.0d0, dbl_mb(k_vso), 1,
     &           dbl_mb(k_vco), 1)

c           ==== Calculates cell function at all grid points ====
            call wfat_collect_celf(geom, nch, natp, int_mb(k_icp),
     &           int_mb(k_ipp), nptp, nptop, ozn, cgrid, chpos0, 
     &           dbl_mb(k_pt), dbl_mb(k_pto), celm, dbl_mb(k_cl),
     &           dbl_mb(k_clo))

c           ==== Calculates mask function at channel center grid points ====
            if (ismask) call wfat_collect_maskf_o(geom, natom, nch, 
     &           cgrid, nptop, chpos0, dbl_mb(k_pto), mask_param,
     &           dbl_mb(k_mao))

            
c           ==== Calculate tau at all grid points ====
            if (pr_srfc .or. pr_rate .or. pr_tmd) then
               call wfat_collect_tau0(geom, nch, nqn, cgrid, chpos0, 
     &              erg, nid, mid, bt, gm, nptp, nptop, zi, 
     &              dbl_mb(k_pt), dbl_mb(k_pto), dcpl_mb(k_tu0),
     &              dcpl_mb(k_tu0o))
               if (wfat1) then
                  call wfat_collect_tau1(geom, nch, nqn, cgrid, chmmt,
     &                 chpos0, erg, nid, mid, bt, gm, nptp, nptop, zi, 
     &                 dbl_mb(k_pt), dbl_mb(k_pto), wfat1_exc,
     &                 dcpl_mb(k_tu1), dcpl_mb(k_tu1o))
               endif
               if (ndid == 0) then
c                  write(luout, '("sum tu1: ", 2es20.10)')
c     &                 sum(dcpl_mb(k_tu1:k_tu1+nptp*nqn*nch-1))
c                  write(luout, '("sum tu1: ", 2es20.10)')
c     &                 dcpl_mb(k_tu1)
               endif
               
                                             
c              ==== Calculate the structure factors ====
               call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &              nptop, natom, int_mb(k_icp), int_mb(k_ipp), ippo, 
     &              dbl_mb(k_vc), dbl_mb(k_vco), dcpl_mb(k_tu0), 
     &              dcpl_mb(k_tu0o), dbl_mb(k_cl), dbl_mb(k_clo),
     &              dbl_mb(k_wt), dbl_mb(k_wto), ismask, dbl_mb(k_ma),
     &              dbl_mb(k_mao), dcpl_mb(k_acg))
               if (wfat1) then
                  do i = 1, nch
                     i_tu0 = k_tu0 + nptp*nqn*(i-1)
                     i_tu0o = k_tu0o + nptop*nqn*(i-1)
                     do k = 1, 3
                        i_vc = (k_vc+nptp*nch) + nptp*3*(i-1) +
     &                         nptp*(k-1)
                        i_vco = (k_vco+nptop*nch) + nptop*3*(i-1) +
     &                          nptop*(k-1)
                        i_acj = k_acj + natom*nqn*3*(i-1) +
     &                          natom*nqn*(k-1)
                        call wfat_calc_acf_nopw(bt, gm, 1, nqn, natp,
     &                       nptp, nptop, natom, int_mb(k_icp), 
     &                       int_mb(k_ipp), ippo, dbl_mb(i_vc), 
     &                       dbl_mb(i_vco), dcpl_mb(i_tu0),
     &                       dcpl_mb(i_tu0o), dbl_mb(k_cl),
     &                       dbl_mb(k_clo), dbl_mb(k_wt), dbl_mb(k_wto),
     &                       ismask, dbl_mb(k_ma), dbl_mb(k_mao),
     &                       dcpl_mb(i_acj))
                     enddo
                  enddo
                  call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &                 nptop, natom, int_mb(k_icp), int_mb(k_ipp), ippo,
     &                 dbl_mb(k_vc), dbl_mb(k_vco), dcpl_mb(k_tu1),
     &                 dcpl_mb(k_tu1o), dbl_mb(k_cl), dbl_mb(k_clo),
     &                 dbl_mb(k_wt), dbl_mb(k_wto), ismask,
     &                 dbl_mb(k_ma), dbl_mb(k_mao), dcpl_mb(k_ack))      ! 4)
c                 * 4) Unlike k_acj, k_ack is not decomposed into its two
c                      orientation-independent components. This is because
c                      tau1 has been calculated as a whole, not decomposed
c                      into the two components above.
                  call wfat_get_srfc1(bt, gm, natom, nch, nqn, 0, mid, 
     &                 erg, chmmt, dcpl_mb(k_acg), dcpl_mb(k_acj),
     &                 dcpl_mb(k_ack), dcpl_mb(k_dum), nsub, maxn, dgs, 
     &                 ndgs, dthr, dps, pbs, .true., dbl_mb(k_sc),
     &                 dcpl_mb(k_sx), dbl_mb(k_pbd), dcpl_mb(k_sg),
     &                 dcpl_mb(k_sj), dcpl_mb(k_sk))    ! 5)
c                 * 5) k_dum takes the place of the second component of
c                      k_ack in the argument list of wfat_get_srfc1.
c                      Related to note 4).
               else
                  call wfat_get_srfc0(geom, bt, gm, nch, erg, nqn, mid,
     &                 0, chmmt, dcpl_mb(k_acg), .true., dcpl_mb(k_sg))
               endif

c              ==== Print structure factor ====
               if (ndid == 0 .and. pr_srfc) then
                  call wfat_print_structuref(natom, bt, gm, nch, nqn,
     &                 dcpl_mb(k_sg), gap)
                  if (wfat1) call wfat_print_srfc1(natom, bt, gm, nch, 
     &                 nqn, dcpl_mb(k_sj), dcpl_mb(k_sk), gap)
               endif
            endif

            if (pr_rate) then
c              ==== Calculate ionization rate ====
               if (wfat1) then
                  call wfat_calc_wfat1_rate(bt, gm, natom, nqn, nch, 
     &                 nid, mid, fieldm, erg, zi, dcpl_mb(k_sg), 
     &                 dcpl_mb(k_sj), dcpl_mb(k_sk), dbl_mb(k_on), 
     &                 chmmt, dbl_mb(k_pbd), dbl_mb(k_fc), wfat1_exc, 
     &                 dbl_mb(k_rt), rate)
               else
                  call wfat_calc_rate(natom, nch, nqn, dcpl_mb(k_sg), 
     &                 dbl_mb(k_fc), dbl_mb(k_rt), rate)
               endif
               
c              ==== Print ionization rate ====
               if (ndid == 0 .and. pr_rate)
     &              call wfat_print_rate(bt, gm, nch, nqn, dbl_mb(k_rt),
     &              rate, gap)
            endif

c     ==== Print TMD ====
            rtdb_stat = rtdb_parallel(.false.)
            if (ndid == 0 .and. pr_tmd) then
               if (wfat1) then
                  call wfat_calc_tmd1(rtdb, bt, gm, natom, nch, nqn,
     &                 nid, mid, orbid0, spinid, erg, chmmt, 
     &                 dbl_mb(k_pbd), zi, fieldm, wfat1_exc,
     &                 dbl_mb(k_on), dcpl_mb(k_sg), dcpl_mb(k_sj), 
     &                 dcpl_mb(k_sk))
               else
                  call wfat_calc_tmd0(rtdb, bt, gm, natom, nch, nqn,
     &                 nid, mid, orbid0, spinid, erg, zi, fieldm,
     &                 dcpl_mb(k_sg))
               endif
            endif
            rtdb_stat = rtdb_parallel(.true.)

            
            nang = nang + 1
            call ga_sync()
            if (nang == pr_chk(kk)) then
               t_now = util_cpusec()
               dt_prog = t_now - t_last
               t_last = t_now
               if(ndid == 0) call wfat_intgl_progress_nopw(kk, nang,
     &              dt_prog, nbt*ngm)
               kk = kk + 1
            endif
         enddo
      enddo
      dt_all = t_now - t_start
      write(dt_str, '(f12.2)') dt_all
      if (ndid == 0) write(luout, '(5x, "The computation of '//
     &     'ionization quantities takes ", a, " s.")')
     &     trim(adjustl(dt_str))

      
c     ==== Terminating printings ====
      if (ndid==0 .and. pr_srfc) then
         call wfat_endprint_structuref()
         if (wfat1) call wfat_endprint_srfc1()
      endif
      if (ndid==0 .and. pr_rate) call wfat_endprint_rate()
      

c     ==== Deallocations ====
      call wfat_erase_rate(l_sg, l_fc, l_rt)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_erase_allgrid(l_pt, l_wt)
      call wfat_erase_allgrid(l_pto, l_wto)
      call wfat_erase_acoef(wfat1, l_acg, l_acj, l_ack, l_dum)
      if (.not. MA_free_heap(l_vf))
     &     call errquit('WFAT: vf_psi deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_vs))
     &     call errquit('WFAT: vs_psi deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_vso))
     &     call errquit('WFAT: vs_psi_o deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc_psi deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_vco))
     &     call errquit('WFAT: vc_psi_o deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_cl))
     &     call errquit('WFAT: celf deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_clo))
     &     call errquit('WFAT: celfo deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_sc0))
     &     call errquit('WFAT: scratch0 deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_ch))
     &     call errquit('WFAT: vch deallocation failed', 3, MA_ERR)
      if (ismask) then
         if (.not. MA_free_heap(l_ma))
     &        call errquit('WFAT: maskf deallocation failed', 3, MA_ERR)
         if (.not. MA_free_heap(l_mao))
     &        call errquit('WFAT: maskfo deallocation failed', 3,
     &        MA_ERR)
      endif
      if (.not. MA_free_heap(l_tu0))
     &     call errquit('WFAT: tau0 deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_tu0o))
     &     call errquit('WFAT: tau0o deallocation failed', 3, MA_ERR)
      if (wfat1) then
         if (.not. MA_free_heap(l_tu1))
     &        call errquit('WFAT: tau1 deallocation failed', 3, MA_ERR)
         if (.not. MA_free_heap(l_tu1o))
     &        call errquit('WFAT: tau1o deallocation failed', 3, MA_ERR)
         call wfat_erase_an(l_sc, l_sx, l_pbd, l_sj, l_sk, l_on)
      endif
      if (.not. MA_free_heap(l_si))
     &     call errquit('WFAT: spinid1 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_pos))
     &     call errquit('WFAT: chpos1 deallocation failed', 61, MA_ERR)

      call ga_sync()
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, nptop,
     &     natom, centid, ptlim, ptlimo, vc_psi, vc_psi_o, tau, tauo, 
     &     celf, celfo, wt, wto, ismask, maskf, maskfo, acf)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      double precision bt                    ! [input], orientation angle in degree.
      double precision gm                    ! [input], orientation angle in degree.
      integer nch
      integer nqn
      integer natp
      integer nptp
      integer nptop
      integer natom
      integer centid(natp)
      integer ptlim(2,natp)
      integer ptlimo(2)
      double precision vc_psi(nptp,nch)
      double precision vc_psi_o(nptop,nch)
      double complex tau(nptp,nqn,nch)
      double complex tauo(nptop,nqn,nch)
      double precision celf(nptp,nch)
      double precision celfo(nptop,nch)
      double precision wt(nptp)
      double precision wto(nptop)
      logical ismask
      double precision maskf(*)           ! [input], if ismask is .true., its dimension must be at least nptp.
      double precision maskfo(*)          ! [input], if ismask is .true., its dimension must be at least nptop*nch.
      double complex acf(natom,nqn,nch)   ! [output]

      integer nac, nac_p, l_acp, k_acp, l_acop, k_acop, g_ac
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


cnopw      call errquit('WFAT: Make sure all instances of call '//
cnopw     &     'wfat_calc_acf_nopw have used the new list of '//
cnopw     &     'arguments', 54, UERR)

      
      call ga_sync()
      acf = czero               ! 1)
c     1) It is necessary to populate acf with zeros in the beginning because
c        the subroutine wfat_assign_acf_nopw_p below does not modify every
c        element of acf (it only modifies elements that correspond to the centers
c        of the current process). This means that the elements that are not
c        modified will have their whatever initial values at the input of this
c        subroutine, which is undesireable.

      nac = natom * nqn * nch
      if (.not. ga_create(mt_dcpl, nac, 1, 'acf', 0, 0, g_ac))
     &     call errquit('WFAT: Cannot create acf global array.',
     &     1, GA_ERR)
      call ga_zero(g_ac)

c     ==== Calculate the asymptotic coefficients contributions ====
c     ====   from the atoms contained in the current process   ====
      nac_p = natp * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nac_p, 'acf_p', l_acp, k_acp))
     &     call errquit('WFAT: acf_p allocation failed', 59, MA_ERR)
      call wfat_calc_acf_nopw_p(nch, nqn, natp, nptp, ptlim, vc_psi, 
     &     tau, celf, wt, ismask, maskf, dcpl_mb(k_acp))

c     ==== Calculate the asymptotic coefficients contributions ====
c     ====   from the channel origin in the current process    ====
      nac_p = nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nac_p, 'acfo_p', l_acop, k_acop))
     &     call errquit('WFAT: acfo_p allocation failed', 59, MA_ERR)
      call wfat_calc_acf_nopw_p(nch, nqn, 1, nptop, ptlimo, vc_psi_o, 
     &     tauo, celfo, wto, ismask, maskfo, dcpl_mb(k_acop))
      
c     ==== Transfer the values of the asymptotic coefficients in ====
c     ====   the current process to their respective positions   ====
c     ====          in the array that holds all atoms            ====
      call wfat_assign_acf_nopw_p(nch, nqn, natp, natom, centid,
     &     dcpl_mb(k_acp), dcpl_mb(k_acop), acf)    ! 2)
c     2) At this point, acf contains the current process' contribution to the
c        asymptotic coefficients. It will be added with the contributions of
c        all the other processes below in the call to ga_acc.

c     ==== Accumulate the asymptotic coefficients across all processes ====
      call ga_sync()
      call ga_acc(g_ac, 1, nac, 1, 1, acf, nac, 1.0d0)      
      call ga_sync()
      acf = czero    ! For safety
      call ga_get(g_ac, 1, nac, 1, 1, acf, nac)       ! 3)
c     3) Now, acf contains the complete (actual) asymptotic coefficients. At
c        at this point, all MPI processes should hold the same numerical
c        values for acf.

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_acp))
     &     call errquit('WFAT: acf_p deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_acop))
     &     call errquit('WFAT: acfo_p deallocation failed', 61, MA_ERR)
      if (.not. ga_destroy(g_ac))
     &     call errquit('WFAT: Cannot destroy acf.', 1, GA_ERR)
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_assign_acf_nopw_p(nch, nqn, natp, natom, centid,
     &     acf_p, acfo_p, acf)

      implicit none
      integer nch
      integer nqn
      integer natp
      integer natom
      integer centid(natp)
      double complex acf_p(natp,nqn,nch)
      double complex acfo_p(nqn,nch)
      double complex acf(natom,nqn,nch)

      integer i, j, icp, ic

      do i = 1, nch
         do j = 1, nqn
            do icp = 1, natp
               ic = centid(icp)
               acf(ic,j,i) = acf_p(icp,j,i)
            enddo
            acf(:,j,i) = acf(:,j,i) + acfo_p(j,i)/natom
         enddo
      enddo

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_calc_acf_nopw_p(nch, nqn, natp, nptp, ptlim, 
     &     vc_psi, tau0, celf, wt, ismask, maskf, acf)

      implicit none
#include "stdio.fh"
#include "global.fh"
      integer nch
      integer nqn
      integer natp
      integer nptp
      integer ptlim(2,natp)
      double precision vc_psi(nptp,nch)
      double complex tau0(nptp,nqn,nch)
      double precision celf(nptp)
      double precision wt(nptp)
      logical ismask
      double precision maskf(*)          ! [input], if ismask is .true., its dimension must be at least nptp.
      double complex acf(natp,nqn,nch)  ! [output]

      integer i, j, ic, k, kk
      double precision mfun
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


      acf = czero
      do i = 1, nch
         do j = 1, nqn

            kk = 0
            do ic = 1, natp
               do k = ptlim(1,ic), ptlim(2,ic)
                  kk = kk + 1
                  if (ismask) then
                     mfun = maskf(kk)
                  else
                     mfun = 1.0d0
                  endif
                  
                  acf(ic,j,i) = acf(ic,j,i) +
     &                          conjg(tau0(kk,j,i)) * vc_psi(kk,i) * 
     &                          mfun * celf(kk) * wt(kk)

cfordebug                  if (ga_nodeid() == 0)
cfordebug     &                 write(luout, '("@grid: ", i8, 2es20.10, '//
cfordebug     &                 '3es20.10)')
cfordebug     &                 kk, tau0(kk,j,i), vc_psi(kk,i), celf(kk), wt(kk)
               enddo
            enddo
            
         enddo
      enddo

      end
c=======================================================
        
      
