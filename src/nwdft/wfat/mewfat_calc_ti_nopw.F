c========================================================
c========================================================
      subroutine mewfat_calc_ti_nopw(rtdb, geom, bas, pr_srfc, pr_rate,
     &     fieldm, nb, sx, sy, n0, n1, ocmo0, ocmo1, pcf, rcf, qcf, scf,
     &     nch, chnid0, ionerg, chmmt, chpos, nqn, nid, mid, zi, gq, 
     &     celm, smh, ismask_e, ergmask, amsk, bmsk, ismask, mask_param)

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
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)             ! [input], no. of electrons in cation.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients returned by mewfat_calc_detcoef. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient returned by mewfat_calc_detcoef.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient returned by mewfat_calc_detcoef. The 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients returned by mewfat_calc_detcoef. Its rows index neutral MO's while its columns index cation MO's.
      integer nch               ! [input], number of ionized channels.
      integer chnid0(nch)       ! [input], orbital ID of ionized channels.
      double precision ionerg(nch)       ! [input]
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

      character*5 si_type
      logical gap, isdft, rs
      integer l_sg, l_fc, l_rt, l_icp, l_ipp, l_pt, l_wt, l_vc, l_cl
      integer l_ma, l_tu0, l_acg
      integer k_sg, k_fc, k_rt, k_icp, k_ipp, k_pt, k_wt, k_vc, k_cl
      integer k_ma, k_tu0, k_acg
      integer natom, n_sg, nfc, nrt, n_pt, n_wt
      integer ndid, nbt, ngm, natp, nptp, i, k, iqn, ig, ib
      double precision bt_start, bt_last, dbt, dbtr, gm_start, gm_last
      double precision nall, dgm, dgmr, dthr, bt, gm, rate
      double precision deg2rad
      parameter(deg2rad=0.01745329251994d0)

      logical nwxc_has_cam, nwxc_has_xc
      external nwxc_has_cam, nwxc_has_xc
      

      if (nch /= 1) call errquit('MEWFAT: At the moment, the argument'//
     &     ' nch of mewfat_calc_ti_nopw can only be equal to 1.',
     &     32, UERR)
      
      
      call ga_sync()
      ndid = ga_nodeid()
      nall = dble(sum(n0))

      
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
      if (.not. MA_alloc_get(MT_DBL, nptp*nch, 'vc_psi', l_vc, k_vc))
     &     call errquit('WFAT: vc_psi allocation failed', 81, MA_ERR)
      rs = nwxc_has_cam()
      isdft = nwxc_has_xc()
      if (isdft) then
         if (.not. rtdb_cget(rtdb, 'mewfat:si', 1, si_type))
     &        call errquit('MEWFAT: Failed reading SI from RTDB.', 17,
     &        RTDB_ERR)
         call mewfat_collect_dft_vcpsi(geom, bas, nb, sx, sy, n0, n1,
     &     nch, nptp, dbl_mb(k_pt), ocmo0, ocmo1, pcf, rcf, qcf, scf, 
     &     zi, chpos, smh, rs, si_type, dbl_mb(k_vc))
      else
         call mewfat_collect_vcpsi(geom, bas, nb, sx, sy, n0, n1,
     &     nch, nptp, dbl_mb(k_pt), ocmo0, ocmo1, pcf, rcf, qcf, scf, 
     &     zi, chpos, smh, dbl_mb(k_vc))
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

c     ==== Allocate tau and as. coeff. at all grid points ====
      if (.not. MA_alloc_get(MT_DCPL, nptp*nqn*nch, 'tau0',
     &     l_tu0, k_tu0))
     &     call errquit('WFAT: tau0 allocation failed', 81, MA_ERR)

c     ==== Obtain field factor ====
      call wfat_get_field_factor(rtdb, nch, nqn, nid, mid, zi, fieldm,
     &     ionerg, dbl_mb(k_fc))
      
c     ==== Initializing the printings ====
      if (ndid == 0)
     &     call wfat_print_field_factor(rtdb, fieldm, nch, nqn, ionerg, 
     &     nid, mid, zi, dbl_mb(k_fc))
      if (ndid == 0 .and. pr_srfc) then
         call wfat_initprint_structuref(nch, chnid0, sx, nqn, nid,
     &        mid)
      endif
      if (ndid == 0 .and. pr_rate)
     &     call wfat_initprint_rate(nch, chnid0, sx, nqn, nid, mid)
      call ga_sync()


      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'WFAT computation of without'//
     &        ' partial wave expansion', 50, .true.)
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
      if (.not. MA_alloc_get(MT_DCPL, natom*nqn*nch, 'acg',
     &     l_acg, k_acg))
     &     call errquit('WFAT: acg allocation failed', 81, MA_ERR)
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
               call wfat_collect_tau0(geom, nch, nqn, chpos, ionerg, 
     &              nid, mid, bt, gm, nptp, zi, dbl_mb(k_pt),
     &              dcpl_mb(k_tu0))
               
c              ==== Calculate the structure factors ====
               call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &              natom, int_mb(k_icp), int_mb(k_ipp), dbl_mb(k_vc), 
     &              dcpl_mb(k_tu0), dbl_mb(k_cl), dbl_mb(k_wt), ismask,
     &              dbl_mb(k_ma), dcpl_mb(k_acg))
               dcpl_mb(k_acg:k_acg+natom*nqn*nch-1) =
     &              sqrt(nall) * dcpl_mb(k_acg:k_acg+natom*nqn*nch-1)
               call wfat_get_srfc0(geom, bt, gm, nch, ionerg, nqn, mid, 
     &              0, chmmt, dcpl_mb(k_acg), .true., dcpl_mb(k_sg))
               
c              ==== Print structure factor ====
               if (ndid == 0 .and. pr_srfc) then
                  call wfat_print_structuref(natom, bt, gm, nch, nqn,
     &                 dcpl_mb(k_sg), gap)
               endif
            endif

            if (pr_rate) then
c              ==== Calculate ionization rate ====
               call wfat_calc_rate(natom, nch, nqn, dcpl_mb(k_sg), 
     &              dbl_mb(k_fc), dbl_mb(k_rt), rate)

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
      endif
      if (ndid==0 .and. pr_rate) call wfat_endprint_rate()

      
c     ==== Deallocations ====
      call wfat_erase_rate(l_sg, l_fc, l_rt)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_erase_allgrid(l_pt, l_wt)
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

      call ga_sync()
      
      end
c========================================================
      
