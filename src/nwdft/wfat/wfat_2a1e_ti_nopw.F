c=======================================================
c=======================================================
      subroutine wfat_2a1e_ti_nopw(rtdb, geom, nqn, nid, mid, gq, zi,
     &     celm, sphi, mch, nr, nt, ke, vcf, erg, chmmt, chpos, orbid0,
     &     spinid, pottype, nel, potparam, pr_srfc, pr_rate, fieldm, 
     &     wfat1, nr_m, nt_m, nr_p, nt_p, cfx, cfy, cfz, pbs, wfat1_exc,
     &     fedvr, nfe, ndvr, maxn, xq, wq, ismask, mask_param)

      implicit none
#include "global.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb                    ! [input],
      integer geom                    ! [input],
      integer nqn                     ! [input],
      integer nid(nqn)                ! [input],
      integer mid(nqn)                ! [input],
      integer gq                      ! [input],
      double precision zi             ! [input],
      integer celm                    ! [input],
      integer sphi                    ! [input],
      integer mch                     ! [input],
      integer nr                      ! [input],
      integer nt                      ! [input],
      double precision ke             ! [input],
      double precision vcf(nr*nt)     ! [input], the minimum number of columns is nch if wfat1 is .false., otherwise 4*nch.
      double precision erg            ! [input],
      double precision chmmt(3)       ! [input],
      double precision chpos(3)       ! [input],
      integer orbid0                  ! [input], orbital ID of ionized channels.
      integer spinid                  ! [input], spin ID of ionized channels.
      character*3 pottype             ! [input],
      integer nel                     ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)    ! [input], 1)
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_rate                 ! [input], if .true. then print ionization rate.
      double precision fieldm         ! [input],
      logical wfat1                   ! [input],
      integer nr_m                    ! [input], only used when wfat1=.true.
      integer nt_m                    ! [input], only used when wfat1=.true.
      integer nr_p                    ! [input], only used when wfat1=.true.
      integer nt_p                    ! [input], only used when wfat1=.true.
      double precision cfx(1:*)       ! [input], only used when wfat1=.true.
      double precision cfy(1:*)       ! [input], only used when wfat1=.true.
      double precision cfz(1:*)       ! [input], only used when wfat1=.true.
      double precision pbs(3)         ! [input]. Only used when wfat1=.true.. The diagonal elements of polarizability tensor (matrix), i.e. the xx, yy, zz components.
      logical wfat1_exc(1:*)          ! [input], only used when wfat1=.true.
      logical fedvr                   ! [input],
      integer nfe                     ! [input], only used when fedvr=.true.
      integer ndvr(*)                 ! [input], only used when fedvr=.true., size is nfe
      integer maxn                    ! [input], only used when fedvr=.true.
      double precision xq(*)          ! [input], only used when fedvr=.true., size is maxn*nfe
      double precision wq(*)          ! [input], only used when fedvr=.true., size is maxn*nfe
      logical ismask                  ! [input],
      double precision mask_param(2)  ! [input]. Only used when ismask=.true.. Must always be supplied regardless.

      logical gap, interm
      integer l_sg, l_fc, l_rt, l_icp, l_ipp, l_pt, l_wt, l_vc, l_cl
      integer l_ma, l_tu0, l_tu1, l_acg, l_acj, l_ack, l_dum, l_sc, l_sx
      integer l_pbd, l_sj, l_sk, l_on
      integer k_sg, k_fc, k_rt, k_icp, k_ipp, k_pt, k_wt, k_vc, k_cl
      integer k_ma, k_tu0, k_tu1, k_acg, k_acj, k_ack, k_dum, k_sc, k_sx
      integer k_pbd, k_sj, k_sk, k_on
      integer natom, n_sg, nfc, nrt, n_pt, n_wt, n_acg, n_acj, n_ack
      integer n_dum, nsc, nsx, npbd, n_sj, n_sk, n_on
      integer ndid, nbt, ngm, natp, nptp, i, k, iqn, i_on, ig, ib, i_tu0
      integer i_vc, i_vcc, i_acj, nsub, maxnsub, dgs, ndgs
      double precision bt_start, bt_last, dbt, dbtr, gm_start, gm_last
      double precision dgm, dgmr, dthr, bt, gm, rate
      double precision deg2rad
      parameter(deg2rad=0.01745329251994d0)

      integer nch, nch1
      parameter(nch=1, nch1=3, nsub=1, maxnsub=1, dgs=1, ndgs=1)
      

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

c     ==== Calculate and gather V_c|psi> at all grid points ====
      if (.not. MA_alloc_get(MT_DBL, nptp*(nch+nch1), 'vc_psi',
     &     l_vc, k_vc))
     &     call errquit('WFAT: vc_psi allocation failed', 81, MA_ERR)
      i_vc = k_vc
      if (wfat1) then
         i_vcc = k_vc + nptp
         if (fedvr) then
            call wfat_2a1efd_collect_vcpsi_ord1(geom, sphi, mch, nfe,
     &           ndvr, maxn, xq, wq, nr, nt, nt_m, nt_p, nptp, pottype,
     &           nel, potparam, zi, vcf, cfx, cfy, cfz, chpos, 
     &           dbl_mb(k_pt), dbl_mb(i_vc), dbl_mb(i_vcc))
         else
            call wfat_2a1e_collect_vcpsi_ord1(geom, sphi, mch, nr, nt,
     &           nr_m, nt_m, nr_p, nt_p, nptp, ke, pottype, nel, 
     &           potparam, zi, vcf, cfx, cfy, cfz, chpos, dbl_mb(k_pt),
     &           dbl_mb(i_vc), dbl_mb(i_vcc))
         endif
      else
         if (fedvr) then
            call wfat_2a1efd_collect_vcpsi(geom, sphi, mch, nfe, ndvr, 
     &           maxn, xq, wq, nr, nt, nptp, pottype, nel, potparam, zi,
     &           vcf, chpos, dbl_mb(k_pt), dbl_mb(i_vc))
         else
            call wfat_2a1e_collect_vcpsi(geom, sphi, mch, nr, nt, nptp, 
     &           ke, pottype, nel, potparam, zi, vcf, chpos,
     &           dbl_mb(k_pt), dbl_mb(i_vc))
         endif
      endif

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
         call wfat_prepare_an(maxnsub, natom, nch, nqn, nsc, l_sc, k_sc,
     &        nsx, l_sx, k_sx, npbd, l_pbd, k_pbd, n_sj, l_sj, k_sj,
     &        n_sk, l_sk, k_sk, n_on, l_on, k_on)
         if (.not. rtdb_get(rtdb, 'wfat:dthr', mt_dbl, 1, dthr))
     &        call errquit('WFAT: Failed reading dthr from rtdb',
     &        1, RTDB_ERR)
         do iqn = 1, nqn
            i_on = k_on + 2*nqn*(1-1) + 2*(iqn-1)
            call wfat_calc_o(nid(iqn), mid(iqn), erg, zi,
     &           dbl_mb(i_on))
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
     &                 dcpl_mb(k_ack), dcpl_mb(k_dum), nsub, maxnsub, 
     &                 dgs, ndgs, dthr, chmmt, pbs, .true., 
     &                 dbl_mb(k_sc), dcpl_mb(k_sx), dbl_mb(k_pbd), 
     &                 dcpl_mb(k_sg), dcpl_mb(k_sj), dcpl_mb(k_sk))
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
      call ga_sync()
      
      end
c=======================================================
