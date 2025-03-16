c========================================================
c========================================================
      subroutine mewfat_calc_ti_nopw(rtdb, geom, bas, pr_srfc, pr_rate,
     &     fieldm, nb, sx, sy, n0, n1, ocmo0, ocmo1, pcf, rcf, qcf, scf,
     &     nch, chnid0, ionerg, chmmt, chpos, nqn, nid, mid, zi, gq, 
     &     celm, smh, ismask_e, ergmask, amsk, bmsk, ismask, mask_param,
     &     wphase)

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
      double complex wphase
      
      character*16 tagt
      character*5 si_type
      logical gap, isdft, rs, cgrid(nch)
      integer l_sg, l_fc, l_rt, l_icp, l_ipp, l_pt, l_wt, l_vc, l_cl
      integer l_ma, l_tu0, l_acg, l_pto, l_wto, l_vco, l_clo, l_mao
      integer l_tu0o
      integer k_sg, k_fc, k_rt, k_icp, k_ipp, k_pt, k_wt, k_vc, k_cl
      integer k_ma, k_tu0, k_acg, k_pto, k_wto, k_vco, k_clo, k_mao
      integer k_tu0o
      integer natom, n_sg, nfc, nrt, n_pt, n_wt, n_pto, n_wto
      integer ndid, nbt, ngm, natp, natop, nptp, nptop, i, k, iqn, ig
      integer ib, icpo, ippo(3), ozn, geomo
      double precision bt_start, bt_last, dbt, dbtr, gm_start, gm_last
      double precision nall, dgm, dgmr, dthr, bt, gm, rate, vdum(3)
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
c     channel center.
      if (.not. geom_destroy(geomo))
     &     call errquit('WFAT: Cannot destroy geomo.', 52, GEOM_ERR)
      call wfat_need_cgrid(rtdb, geom, nch, chpos, cgrid)
      call wfat_print_gridpt(rtdb, geom, gq, ozn, nch, cgrid, chpos)
      
      
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
     &        nch, nptp, nptop, dbl_mb(k_pt), dbl_mb(k_pto), cgrid, 
     &        ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, chpos, smh, rs, 
     &        si_type, dbl_mb(k_vc), dbl_mb(k_vco))
      else
         call mewfat_collect_hf_vcpsi(geom, bas, nb, sx, sy, n0, n1,
     &        nch, nptp, nptop, dbl_mb(k_pt), dbl_mb(k_pto), cgrid, 
     &        ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, chpos, smh, 
     &        dbl_mb(k_vc), dbl_mb(k_vco))
      endif

c     ==== Calculate and gather the cell ====
c     ====  function at all grid points  ====
      if (.not. MA_alloc_get(MT_DBL, nptp, 'celf', l_cl, k_cl))
     &     call errquit('WFAT: celf allocation failed', 98, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptop*nch, 'celfo',
     &     l_clo, k_clo))
     &     call errquit('WFAT: celfo allocation failed', 98, MA_ERR)
      call wfat_collect_celf(geom, nch, natp, int_mb(k_icp),
     &     int_mb(k_ipp), nptp, nptop, ozn, cgrid, chpos, dbl_mb(k_pt),
     &     dbl_mb(k_pto), celm, dbl_mb(k_cl), dbl_mb(k_clo))

c     ==== Calculate and gather the mask ====
c     ====  function at all grid points  ====
      if (ismask) then
         if (.not. MA_alloc_get(MT_DBL, nptp, 'maskf', l_ma, k_ma))
     &        call errquit('WFAT: maskf allocation failed', 98, MA_ERR)
         if (.not. MA_alloc_get(MT_DBL, nptop*nch, 'maskfo',
     &        l_mao, k_mao))
     &        call errquit('WFAT: maskfo allocation failed', 98, MA_ERR)
         call wfat_collect_maskf(geom, natom, nch, cgrid, nptp, nptop,
     &        chpos, dbl_mb(k_pt), dbl_mb(k_pto), mask_param, 
     &        dbl_mb(k_ma), dbl_mb(k_mao))
      endif

c     ==== Allocate tau and as. coeff. at all grid points ====
      if (.not. MA_alloc_get(MT_DCPL, nptp*nqn*nch, 'tau0',
     &     l_tu0, k_tu0))
     &     call errquit('WFAT: tau0 allocation failed', 98, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nptop*nqn*nch, 'tau0o',
     &     l_tu0o, k_tu0o))
     &     call errquit('WFAT: tau0o allocation failed', 98, MA_ERR)

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
               call wfat_collect_tau0(geom, nch, nqn, cgrid, chpos,
     &              ionerg, nid, mid, bt, gm, nptp, nptop, zi, 
     &              dbl_mb(k_pt), dbl_mb(k_pto), dcpl_mb(k_tu0),
     &              dcpl_mb(k_tu0o))
               
c              ==== Calculate the structure factors ====
               call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &              nptop, natom, int_mb(k_icp), int_mb(k_ipp), ippo, 
     &              dbl_mb(k_vc), dbl_mb(k_vco), dcpl_mb(k_tu0), 
     &              dcpl_mb(k_tu0o), dbl_mb(k_cl), dbl_mb(k_clo),
     &              dbl_mb(k_wt), dbl_mb(k_wto), ismask, dbl_mb(k_ma),
     &              dbl_mb(k_mao), dcpl_mb(k_acg))
               dcpl_mb(k_acg:k_acg+natom*nqn*nch-1) = wphase *
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
      call wfat_erase_allgrid(l_pto, l_wto)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc_psi deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_vco))
     &     call errquit('WFAT: vc_psi_o deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_cl))
     &     call errquit('WFAT: celf deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_clo))
     &     call errquit('WFAT: celfo deallocation failed', 3, MA_ERR)
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

      call ga_sync()
      
      end
c========================================================
      
