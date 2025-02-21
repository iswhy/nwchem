c=======================================================
c=======================================================
      subroutine wfat_calc_ti_nopw(rtdb, geom, bas, pr_srfc, pr_rate,
     &     fieldm, nb, ns, rdm, nch, orbid0, vch, spinid, erg, chmmt, 
     &     chpos, nqn, nid, mid, zi, gq, celm, smh, ismask_e, ergmask,
     &     amsk, bmsk, ismask, mask_param, wfat1, nsub, maxn, dgs, ndgs,
     &     dps, pbs, wfat1_exc)

      implicit none
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
      double precision fieldm
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
      
      logical gap, interm
      integer l_sg, l_fc, l_rt, l_icp, l_ipp, l_pt, l_wt, l_vc, l_cl
      integer l_ma, l_tu0, l_tu1, l_acg, l_acj, l_ack, l_dum, l_sc, l_sx
      integer l_pbd, l_sj, l_sk, l_on, l_si, l_pos
      integer k_sg, k_fc, k_rt, k_icp, k_ipp, k_pt, k_wt, k_vc, k_cl
      integer k_ma, k_tu0, k_tu1, k_acg, k_acj, k_ack, k_dum, k_sc, k_sx
      integer k_pbd, k_sj, k_sk, k_on, k_si, k_pos
      integer natom, n_sg, nfc, nrt, n_pt, n_wt, n_acg, n_acj, n_ack
      integer n_dum, nsc, nsx, npbd, n_sj, n_sk, n_on, nch1
      integer ndid, nbt, ngm, natp, nptp, i, k, iqn, i_on, ig, ib, i_tu0
      integer i_vc, i_acj
      double precision bt_start, bt_last, dbt, dbtr, gm_start, gm_last
      double precision dgm, dgmr, dthr, bt, gm, rate
      double precision deg2rad
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

c     ==== Calculate and gather grid points ====
      call wfat_alloc_gridrange(natom, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natp, nptp)
      call wfat_prep_allgrid(nptp, n_pt, l_pt, k_pt, n_wt, l_wt, k_wt)
      call wfat_collect_mfgrid(geom, gq, natp, nptp, int_mb(k_icp),
     &     int_mb(k_ipp), dbl_mb(k_pt), dbl_mb(k_wt))

c     ==== Construct new spinid and chpos that can ====
c     ====  take care of wfat1 when it is .true.   ====
      nch1 = 0
      if (wfat1) nch1 = 3 * nch
      if (.not. MA_alloc_get(MT_INT, nch+nch1, 'spinid1', l_si, k_si))
     &     call errquit('WFAT: spinid1 allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, 3*(nch+nch1), 'chpos1',
     &     l_pos, k_pos))
     &     call errquit('WFAT: chpos1 allocation failed', 59, MA_ERR)
      call wfat_get_wfat1_sppos(wfat1, nch, spinid, chpos, int_mb(k_si),
     &     dbl_mb(k_pos))

c     ==== Calculate and gather V_c|psi> at all grid points ====
      if (.not. MA_alloc_get(MT_DBL, nptp*(nch+nch1), 'vc_psi',
     &     l_vc, k_vc))
     &     call errquit('WFAT: vc_psi allocation failed', 81, MA_ERR)
      call wfat_collect_vcpsi(geom, bas, nb, ns, nch+nch1, nptp, vch,
     &     int_mb(k_si), dbl_mb(k_pos), dbl_mb(k_pt), rdm, zi, smh,
     &     ismask_e, ergmask, amsk, bmsk, dbl_mb(k_vc))

c     ==== Calculate and gather the cell ====
c     ====  function at all grid points  ====
      if (.not. MA_alloc_get(MT_DBL, nptp, 'celf', l_cl, k_cl))
     &     call errquit('WFAT: celf allocation failed', 81, MA_ERR)
      call wfat_collect_celf(geom, natp, nptp, int_mb(k_icp), 
     &     int_mb(k_ipp), celm, dbl_mb(k_pt), dbl_mb(k_cl))

c     ==== Calculate and gather the mask ====
c     ====  function at all grid points  ====
      if (ismask) then
         if (.not. MA_alloc_get(MT_DBL, nptp, 'maskf', l_ma, k_ma))
     &        call errquit('WFAT: maskf allocation failed', 81, MA_ERR)
         call wfat_collect_maskf(geom, natom, nptp, dbl_mb(k_pt),
     &        mask_param, dbl_mb(k_ma))
      endif

c     ==== Allocate tau at all grid points and as. coeff. ====
      if (.not. MA_alloc_get(MT_DCPL, nptp*nqn*nch, 'tau0',
     &     l_tu0, k_tu0))
     &     call errquit('WFAT: tau0 allocation failed', 81, MA_ERR)
      if (wfat1) then
         if (.not. MA_alloc_get(MT_DCPL, nptp*nqn*nch, 'tau1',
     &        l_tu1, k_tu1))
     &        call errquit('WFAT: tau1 allocation failed', 81, MA_ERR)
      endif
      call wfat_prep_acoef(natom, nqn, nch, wfat1, n_acg, l_acg,
     &     k_acg, n_acj, l_acj, k_acj, n_ack, l_ack, k_ack, n_dum,
     &     l_dum, k_dum)

      
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
      if (ndid == 0)
     &     call wfat_print_field_factor(rtdb, fieldm, nch, nqn, erg, 
     &     nid, mid, zi, dbl_mb(k_fc))
      if (ndid == 0 .and. pr_srfc) then
         call wfat_initprint_structuref(nch, orbid0, spinid, nqn, nid,
     &        mid)
         if (wfat1) call wfat_initprint_srfc1(nch, orbid0, spinid, nqn,
     &        nid, mid)
      endif
      if (ndid == 0 .and. pr_rate)
     &     call wfat_initprint_rate(nch, orbid0, spinid, nqn, nid, mid)
      call ga_sync()

      
      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'WFAT computation without '//
     &        'partial wave expansion', 50, .true.)
         write(luout, '("")')
         write(luout, 901)
         write(luout, '("")')
      endif
 901  format(
     &     5x, 'Calculation of some WFAT quantities using explicit ',
     &     'integral formula (without partial wave expansion) has ', /,
     &     5x, 'started. Beware that this is often a slow computation,',
     &     ' so be patient. Use a partial wave expansion algorithm ', /,
     &     5x, 'if you want fast results by specifying a positive ',
     &     'value for the input keyword lmax.')

      
c     ==== Scan over the orientation angles ====
      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)

         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)

            if (ib /= nbt) then
               gap = .false.
            else
               gap = .true.
            endif
            
c           ==== Calculate tau ====
            if (pr_srfc .or. pr_rate) then
               call wfat_collect_tau0(geom, nch, nqn, chpos, erg, nid, 
     &              mid, bt, gm, nptp, zi, dbl_mb(k_pt), dcpl_mb(k_tu0))
               if (wfat1) then
                  call wfat_collect_tau1(geom, nch, nqn, chmmt, chpos, 
     &                 erg, nid, mid, bt, gm, nptp, zi, dbl_mb(k_pt),
     &                 wfat1_exc, dcpl_mb(k_tu1))
               endif
               if (ndid == 0) then
c                  write(luout, '("sum tu1: ", 2es20.10)')
c     &                 sum(dcpl_mb(k_tu1:k_tu1+nptp*nqn*nch-1))
c                  write(luout, '("sum tu1: ", 2es20.10)')
c     &                 dcpl_mb(k_tu1)
               endif
               
                                             
c              ==== Calculate the structure factors ====
               call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &              natom, int_mb(k_icp), int_mb(k_ipp), dbl_mb(k_vc), 
     &              dcpl_mb(k_tu0), dbl_mb(k_cl), dbl_mb(k_wt), ismask,
     &              dbl_mb(k_ma), dcpl_mb(k_acg))               
               if (wfat1) then
                  do i = 1, nch
                     i_tu0 = k_tu0 + nptp*nqn*(i-1)
                     do k = 1, 3
                        i_vc = (k_vc+nptp*nch) + nptp*3*(i-1) +
     &                         nptp*(k-1)
                        i_acj = k_acj + natom*nqn*3*(i-1) +
     &                          natom*nqn*(k-1)
                        call wfat_calc_acf_nopw(bt, gm, 1, nqn, natp,
     &                       nptp, natom, int_mb(k_icp), int_mb(k_ipp), 
     &                       dbl_mb(i_vc), dcpl_mb(i_tu0), 
     &                       dbl_mb(k_cl), dbl_mb(k_wt), ismask, 
     &                       dbl_mb(k_ma), dcpl_mb(i_acj))
                     enddo
                  enddo
                  call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &                 natom, int_mb(k_icp), int_mb(k_ipp), 
     &                 dbl_mb(k_vc), dcpl_mb(k_tu1), dbl_mb(k_cl), 
     &                 dbl_mb(k_wt), ismask, dbl_mb(k_ma),
     &                 dcpl_mb(k_ack))
                  call wfat_get_srfc1(bt, gm, natom, nch, nqn, 0, mid, 
     &                 erg, chmmt, dcpl_mb(k_acg), dcpl_mb(k_acj),
     &                 dcpl_mb(k_ack), dcpl_mb(k_dum), nsub, maxn, dgs, 
     &                 ndgs, dthr, dps, pbs, .true., dbl_mb(k_sc),
     &                 dcpl_mb(k_sx), dbl_mb(k_pbd), dcpl_mb(k_sg),
     &                 dcpl_mb(k_sj), dcpl_mb(k_sk))
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
     &                 interm, dbl_mb(k_rt), rate)
               else
                  call wfat_calc_rate(natom, nch, nqn, dcpl_mb(k_sg), 
     &                 dbl_mb(k_fc), dbl_mb(k_rt), rate)
               endif
               
c              ==== Print ionization rate ====
               if (ndid == 0 .and. pr_rate)
     &              call wfat_print_rate(bt, gm, nch, nqn, dbl_mb(k_rt),
     &              rate, gap)
            endif
            
         enddo
      enddo

      
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
      call wfat_erase_acoef(wfat1, l_acg, l_acj, l_ack, l_dum)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc_psi deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_cl))
     &     call errquit('WFAT: celf deallocation failed', 3, MA_ERR)
      if (ismask) then
         if (.not. MA_free_heap(l_ma))
     &        call errquit('WFAT: maskf deallocation failed', 3, MA_ERR)
      endif
      if (.not. MA_free_heap(l_tu0))
     &     call errquit('WFAT: tau0 deallocation failed', 3, MA_ERR)
      if (wfat1) then
         if (.not. MA_free_heap(l_tu1))
     &        call errquit('WFAT: tau1 deallocation failed', 3, MA_ERR)
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
      subroutine wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &     natom, centid, ptlim, vc_psi, tau0, celf, wt, ismask, maskf,
     &     acf)

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
      integer natom
      integer centid(natp)
      integer ptlim(2,natp)
      double precision vc_psi(nptp,nch)
      double complex tau0(nptp,nqn,nch)
      double precision celf(nptp)
      double precision wt(nptp)
      logical ismask
      double precision maskf(*)      ! [input], if ismask is .true., its dimension must be at least nptp.
      double complex acf(natom,nqn,nch)   ! [output]

      integer nac, nac_p, l_acp, k_acp, g_ac
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


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

c     ==== Calculate the asymptotic coefficients for the ====
c     ====    atoms contained in the current process     ====
      nac_p = natp * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nac_p, 'acf_p', l_acp, k_acp))
     &     call errquit('WFAT: acf_p allocation failed', 59, MA_ERR)
      call wfat_calc_acf_nopw_p(nch, nqn, natp, nptp, ptlim, 
     &     vc_psi, tau0, celf, wt, ismask, maskf, dcpl_mb(k_acp))

c     ==== Transfer the values of the asymptotic coefficients in ====
c     ====   the current process to their respective positions   ====
c     ====          in the array that holds all atoms            ====
      call wfat_assign_acf_nopw_p(nch, nqn, natp, natom, centid,
     &     dcpl_mb(k_acp), acf)    ! 2)
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
      if (.not. ga_destroy(g_ac))
     &     call errquit('WFAT: Cannot destroy acf.', 1, GA_ERR)
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_assign_acf_nopw_p(nch, nqn, natp, natom, centid,
     &     acf_p, acf)

      implicit none
      integer nch
      integer nqn
      integer natp
      integer natom
      integer centid(natp)
      double complex acf_p(natp,nqn,nch)
      double complex acf(natom,nqn,nch)

      integer i, j, icp, ic

      do i = 1, nch
         do j = 1, nqn
            do icp = 1, natp
               ic = centid(icp)
               acf(ic,j,i) = acf_p(icp,j,i)
            enddo
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

      
c=======================================================
c=======================================================
      subroutine wfat_collect_celf(geom, natp, nptp, centid, ptlim,
     &     celm, v_m, celf)

      implicit none
#include "errquit.fh"
      integer geom
      integer natp
      integer nptp
      integer centid(natp)
      integer ptlim(2,natp)
      integer celm
      double precision v_m(3,nptp)
      double precision celf(nptp)

      logical t1
      integer icp, ic, i, j
      double precision v_c(3)

      j = 0
      do icp = 1, natp
         ic = centid(icp)
         do i = ptlim(1,icp), ptlim(2,icp)
            j = j + 1
            call geom_cent2cent(geom, 0, v_m(:,j), ic, v_c, t1)
            call wfat_celf_1pt(geom, ic, v_c, celm, celf(j))
         enddo
      enddo

      if (j /= nptp) call errquit('WFAT: The requirement j=nptp is '//
     &     'not satisfied. This is a bug, check the wfat_collect_celf'//
     &     ' subroutine file.', 99, UERR)
      
      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_collect_maskf(geom, natom, nptp, v_m, mask_param,
     &     maskf)

      implicit none
#include "errquit.fh"
      integer geom
      integer natom
      integer nptp
      double precision v_m(3,nptp)
      double precision mask_param(2)
      double precision maskf(nptp)

      integer i
      double complex maskf_

      do i = 1, nptp
         maskf_ = cmplx(1.0d0, 0.0d0, kind=kind(1.0d0))
         call wfat_mask_intd(geom, natom, 1, 1, 1, v_m(:,i), 
     &        mask_param, maskf_)
         call errquit('ERROR: No-partial-wave algorithm cannot be '//
     &        'used at the moment. First, change the dimension of v_m'//
     &        ' to also include nch.', 55, UERR)
         maskf(i) = real(maskf_)
      enddo

      end
c=======================================================
      
