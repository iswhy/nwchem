c=======================================================
c=======================================================
      subroutine wfat_sometests(rtdb, zi, nch, erg, lmax_, chmmt)

      implicit none
#include "rtdb.fh"      
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb
      double precision zi
      integer nch
      double precision erg(nch)
      integer lmax_
      double precision chmmt(3, nch)

      integer ndid
      integer l_chgtest(2), l_r1test(2), l_qtest(2), l_dvrttest(3)
      integer l_phi1test(2), l_tautest(2), l_otest(2), l_fbrrtest(3)
      integer l_altest(2), l_d2lagtest(2), l_2a1etest(7), l_dvrrtest(3)
      integer l_3jtest(4), l_3j0test(2), l_gntest(4)
      double precision d_chgutest(2), d_erfctest(5)


      ndid = ga_nodeid()

      if (rtdb_get(rtdb, 'wfat:test_chg', mt_int, 2, l_chgtest)) then
         if (ndid==0) call wfat_test_chg(l_chgtest(1), l_chgtest(2),
     &        zi, erg(1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_chgu', mt_dbl, 2, d_chgutest)) then
         if (ndid==0) call wfat_test_chgu(d_chgutest(1), d_chgutest(2))
      endif
      if (rtdb_get(rtdb, 'wfat:test_r1', mt_int, 2, l_r1test)) then
         if (ndid==0) call wfat_test_r1(l_r1test(1), l_r1test(2),
     &        erg(1), zi, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_q', mt_int, 2, l_qtest)) then
         if (ndid==0) call wfat_test_q(l_qtest(1), l_qtest(2),
     &        erg(1), zi, lmax_, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_phi1', mt_int, 2, l_phi1test)) then
         if (ndid==0) call wfat_test_phi1xi(l_phi1test(1), 
     &        l_phi1test(2), erg(1), zi, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_tau', mt_int, 2, l_tautest)) then
         if (ndid==0) call wfat_test_tau(l_tautest(1), l_tautest(2),
     &        erg(1), zi, lmax_, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test2_tau', mt_int, 2, l_tautest)) then
         if (ndid==0) call wfat_test2_tau(l_tautest(1), l_tautest(2),
     &        erg(1), zi, lmax_, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_o', mt_int, 2, l_otest)) then
         if (ndid==0) call wfat_test_o(l_otest(1), l_otest(2),
     &        erg(1), zi)
      endif
      if (rtdb_get(rtdb, 'wfat:test_erfc', mt_dbl, 5, d_erfctest)) then
         if (ndid==0)
     &        call wfat_test_erfc(d_erfctest(1),
     &                            d_erfctest(2),
     &                            d_erfctest(3),
     &                            d_erfctest(4),
     &                            idnint(d_erfctest(5)))
      endif
      if (rtdb_get(rtdb, 'wfat:test_assocleg', mt_int, 2, l_altest))
     &     then
         if (ndid==0) call wfat_test_assocleg(l_altest(1), l_altest(2))
      endif
      if (rtdb_get(rtdb, 'wfat:test_d2lag', mt_int, 2, l_d2lagtest))
     &     then
         if (ndid==0)
     &        call wfat_test_d2lag(l_d2lagtest(1), l_d2lagtest(2))
      endif
      if (rtdb_get(rtdb, 'wfat:test_fbrr', mt_int, 3, l_fbrrtest))
     &     then
         if (ndid==0)
     &        call wfat_test_fbrr(l_fbrrtest(1), l_fbrrtest(2),
     &        l_fbrrtest(3))
      endif
      if (rtdb_get(rtdb, 'wfat:test_dvrr', mt_int, 3, l_dvrrtest))
     &     then
         if (ndid==0)
     &        call wfat_test_dvrr(l_dvrrtest(1), l_dvrrtest(2),
     &        l_dvrrtest(3))
      endif
      if (rtdb_get(rtdb, 'wfat:test_dvrt', mt_int, 3, l_dvrttest))
     &     then
         if (ndid==0)
     &        call wfat_test_dvrt(l_dvrttest(1), l_dvrttest(2),
     &        l_dvrttest(3))
      endif
      if (rtdb_get(rtdb, 'wfat:test_3j', mt_int, 4, l_3jtest))
     &     then
         if (ndid==0)
     &        call wfat_test_3j(l_3jtest(1), l_3jtest(2), l_3jtest(3),
     &        l_3jtest(4))
      endif
      if (rtdb_get(rtdb, 'wfat:test_3j000', mt_int, 2, l_3j0test))
     &     then
         if (ndid==0)
     &        call wfat_test_3j000(l_3j0test(1), l_3j0test(2))
      endif
      if (rtdb_get(rtdb, 'wfat:test_gaunt', mt_int, 4, l_gntest))
     &     then
         if (ndid==0)
     &        call wfat_test_gaunt(l_gntest(1), l_gntest(2),
     &        l_gntest(3), l_gntest(4))
      endif
      if (rtdb_get(rtdb, 'wfat:test_2a1e', mt_int, 7, l_2a1etest))
     &     then
         if (ndid==0)
     &        call wfat_test_2a1e(l_2a1etest(1), l_2a1etest(2),
     &        l_2a1etest(3), l_2a1etest(4), l_2a1etest(5),
     &        l_2a1etest(6), l_2a1etest(7))
      endif
      

      call ga_sync()
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_test()

      implicit none
#include "global.fh"      
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
c#include "cscf.fh"
#include "wfat.fh"
#include "wfat_qm_param_test.fh"

         
      integer nch, nqn, gqual
      parameter ( nch=2, nqn=3, gqual=2 )
      integer l_in, l_cmm, l_nd, l_sf, l_fc, l_iy
      integer k_in, k_cmm, k_nd, k_sf, k_fc, k_iy


      write(luout, '("------------")')
      write(luout, '("Testing WFAT")')
      write(luout, '("------------")')
      call wfat_init_test()



      
c      write(luout, '("movecs_in = ", a256)') trim(adjustl(movecs_in))
c      write(luout, '("movecs_out = ", a256)') trim(adjustl(movecs_out))
c
c     ==== Checking numerical density ====
      if (.not. MA_alloc_get(MT_DBL, 2, 'numden', l_nd, k_nd))
     &     call errquit('WFAT: numden allocation failed', 59, MA_ERR)
      dbl_mb(k_nd : k_nd+1) = 0.0d0
      call wfat_numdens(bas_hd, geom_hd, nspin0, gqual,
     &     dbl_mb(k_nd:k_nd+1))
      write(luout, '("Numerically integrated density in WFAT test = "'//
     &             ', 2f12.6)') dbl_mb(k_nd : k_nd+1)

      if (.not. MA_free_heap(l_nd))
     &     call errquit('WFAT: numden deallocation failed',
     &     3, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_test_mo()
      
      implicit none
#include "bas.fh"      
#include "cdft.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer nm0, nm1, nm2, nm3
      integer pm0, pm1, pm2, pm3
      integer nscr, nmov, nbas, is, io, l_s, k_s, l_mo, k_mo
      double precision vpt(3)
      

      if (.not. bas_numbf(ao_bas_han, nbas))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)

      nscr = nbas * (20 + nbas)
      if (.not. MA_alloc_get(mt_dbl, nscr, 'scr', l_s, k_s))
     &     call errquit('WFAT: Allocating scr failed', 53, MA_ERR)
      nmov = 20 * nbas
      if (.not. MA_alloc_get(mt_dbl, nmov, 'mov', l_mo, k_mo))
     &     call errquit('WFAT: Allocating mov failed', 53, MA_ERR)

      nm0 =  1 * nbas
      nm1 =  3 * nbas
      nm2 =  6 * nbas
      nm3 = 10 * nbas

      pm0 = k_mo          !! molecular orbitals
      pm1 = pm0 + nm0     !! 1st derivative of molecular orbitals
      pm2 = pm1 + nm1     !! 2nd derivative of molecular orbitals
      pm3 = pm2 + nm2     !! 3rd derivative of molecular orbitals

      vpt = (/ 0.4d0, 1.0d0, -2.3d0 /)
      do is = 1, ipol
         call wfat_eval_mo_1pt(geom, ao_bas_han, is, 0, vpt,
     &        dbl_mb(k_s), dbl_mb(pm0), dbl_mb(pm1), dbl_mb(pm2),
     &        dbl_mb(pm3))
         do io = 1, nbas
            write(luout, '("is/MO value = ", i1, i5, e14.5)') is, io,
     &           dbl_mb(pm0+io-1)
         enddo
      enddo

      if (.not. ma_free_heap(l_s))
     &     call errquit('WFAT: Deallocating scr failed', 2, MA_ERR)
      if (.not. ma_free_heap(l_mo))
     &     call errquit('WFAT: Deallocating mov failed', 2, MA_ERR)

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_test_fock_2e(rtdb)

      implicit none
#include "global.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "stdio.fh"      
      integer rtdb
      
      logical oskel_, rtdb_stat
      integer g_rd(1), g_j(1), nb
      double precision jfac(1), kfac(1), tol2e_


c      rtdb_stat = rtdb_parallel(.false.)


      if (.not. bas_numbf(ao_bas_han, nb))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)

      if (.not. rtdb_get(rtdb, 'dft:skeleton', mt_log, 1, oskel_))
     &     call errquit('WFAT: Failed reading oskel from RTDB.', 1,
     &     RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'dft:tol2e', mt_dbl, 1, tol2e_))
     &     call errquit('WFAT: Failed reading tol2e from RTDB.', 1,
     &     RTDB_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, 'rdm_mat', 0, 0, g_rd(1)))
     &     call errquit('WFAT: Cannot create ''hartree_mat'' global '//
     &     'array.', 1, GA_ERR)
      call ga_zero(g_rd(1))
      call ga_add_constant(g_rd(1), 5.0d0)
      
      if (.not. ga_create(mt_dbl, nb, nb, 'hartree_mat', 0, 0, g_j(1)))
     &     call errquit('WFAT: Cannot create ''hartree_mat'' global '//
     &     'array.', 1, GA_ERR)
      call ga_zero(g_j(1))

      jfac(1) = 1.0d0
      kfac(1) = 0.0d0
      call fock_2e(geom, ao_bas_han, 1, jfac, kfac, tol2e_, oskel_,
     &     g_rd, g_j, .false.)

      write(luout, '("Test J matrix : ")')
      call ga_print(g_j(1))


      if (.not. ga_destroy(g_rd(1)))
     &     call errquit('WFAT: Cannot destroy rdm_mat.', 1, GA_ERR)
      if (.not. ga_destroy(g_j(1)))
     &     call errquit('WFAT: Cannot destroy hartree_mat.', 1, GA_ERR)

c      rtdb_stat = rtdb_parallel(.true.)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_test_fock(rtdb, geom, bas, nb, ns, ne, g_mo,
     &           chg, mult, tol2e, oskel, ixcfit, iolgc)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb            ! [input]
      integer geom            ! [input]
      integer bas             ! [input]
      integer nb              ! [input]
      integer ns              ! [input]
      integer ne(ns)          ! [input]
      integer g_mo(ns)        ! [input]
      double precision chg    ! [input]
      integer mult            ! [input]
      double precision tol2e  ! [input]
      logical oskel           ! [input]
      integer ixcfit          ! [input]
      logical iolgc           ! [input]

      character*9 hfx_name(2), vxcname(2)      
      integer g_fk(ns), g_tmp, g_ev
      integer is, ndid

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      
      ndid = ga_nodeid()
      if (ndid==0) write(luout, '("Testing Fock building...")')
      hfx_name = (/ 'hfx_alpha', 'hfx_beta' /)

      do is = 1, ns
         g_fk(is) = ga_create_atom_blocked(geom, bas, hfx_name(is))
         call ga_zero(g_fk(is))
      enddo
      call wfat_get_fock(rtdb, geom, bas, nb, ns, ne, g_mo, chg, mult, 
     &     tol2e, oskel, ixcfit, iolgc, .true., g_fk)
      

c     ==== Calculate (total) Fock matrices and print their eigvals ====
      g_tmp = ga_create_atom_blocked(geom, bas, 'tmp')
      g_ev = ga_create_atom_blocked(geom, bas, 'eigval')
      do is = 1, ns
         call ga_zero(g_tmp)
         call ga_zero(g_ev)
         call ga_dgemm('n', 'n', nb, nb, nb, 1.0d0, g_fk(is), g_mo(is),
     &        0.0d0, g_tmp)
         call ga_dgemm('t', 'n', nb, nb, nb, 1.0d0, g_mo(is), g_tmp,
     &        0.0d0, g_ev)

c        ==== Printing eigenvalues ====
         if (ndid == 0) then
            write(luout, '("Printing eigenvalues for spin: ", i1)') is
         endif
         call ga_sync()
         call ga_print(g_ev)
      enddo


c     ==== Deallocations ====
      do is = 1, ns
         if (.not. ga_destroy(g_fk(is)))
     &        call errquit('wfat_test_hamiltonian: Cannot destroy '//
     &        hfx_name(is), 1 , GA_ERR)
      enddo
      if (.not. ga_destroy(g_tmp))
     &     call errquit('wfat_test_hamiltonian: Cannot destroy tmp.',
     &     1 , GA_ERR)
      if (.not. ga_destroy(g_ev))
     &     call errquit('wfat_test_hamiltonian: Cannot destroy '//
     &     'eigval.', 1 , GA_ERR)
      
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_test_numovl(geom, bas, nz, zp, nx, xlims,
     &           ny, ylims, g_rh, oskel)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "bas.fh"
#include "stdio.fh"
#include "global.fh"      
      integer geom                    ! [input], geometry handle using which bas handle was created.
      integer bas                     ! [input], basis handle in which vcf is represented.
      integer nz                      ! [input], number of z points to evaluate the integrated density.
      double precision zp(nz)         ! [input], the z points at which to evaluate the integrated density.
      integer nx                      ! [input], number of integration points along x.
      double precision xlims(2)       ! [input], integration limits along x.
      integer ny                      ! [input], number of integration points along y.
      double precision ylims(2)       ! [input], integration limits along y.
      integer g_rh                    ! [output], the xy-integrated density.
c                                       MUST BE declared such that IT HAS ONLY ONE COLUMN.
      logical oskel
      
      integer nb, nz0
      integer ndid, iz_in, iz_fin, iz, i, col1, col2
      integer g_ov, g_ecov, l_ov, k_ov, l_ov0, k_ov0
      double precision dz

      integer ga_create_atom_blocked
      external ga_create_atom_blocked


      call ga_sync()

      ndid = ga_nodeid()
      if (.not. bas_numbf(bas, nb))
     &     call errquit('WFAT: Getting nb failed in '//
     &     'wfat_vecdens_numinteg_xy', 1, BASIS_ERR)
      
      

c     ==== Obtain range of indices of zp for the present process ====
      call ga_distribution(g_rh, ndid, iz_in, iz_fin, col1, col2)
      if ((col1 /= 1) .or. (col1 /= col2))
     &     call errquit('WFAT: Argument g_rh of wfat_vecdens_numinteg'//
     &     '_xy must be declared such that it has one column.', 1, UERR)
      nz0 = iz_fin - iz_in + 1
      dz = (zp(iz_fin) - zp(iz_in)) / (nz0 - 1)


c     ==== Calculate numerically integrated overlap matrix ====
      g_ov = ga_create_atom_blocked(geom, bas, 'ga_ovl')
      call ga_zero(g_ov)
      if (.not. ma_alloc_get(mt_dbl, nb*nb, 'overlap_z', l_ov0, k_ov0))
     &     call errquit('WFAT: Cannot allocate overlap_z', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nb*nb, 'overlap', l_ov, k_ov))
     &     call errquit('WFAT: Cannot allocate overlap', 1, MA_ERR)
      dbl_mb(k_ov:k_ov+nb*nb-1) = 0.0d0
      
      do iz = iz_in, iz_fin
         call wfat_overlap_numinteg_xy(geom, bas, nx, xlims, ny, ylims, 
     &        zp(iz), dbl_mb(k_ov0))
         dbl_mb(k_ov:k_ov+nb*nb-1) = dbl_mb(k_ov:k_ov+nb*nb-1) +
     &                               dbl_mb(k_ov0:k_ov0+nb*nb-1) * dz
      enddo
      call ga_acc(g_ov, 1, nb, 1, nb, dbl_mb(k_ov), nb, 1.0d0)

      
c     ==== Calculate analytically integrated overlap matrix ====
      g_ecov = ga_create_atom_blocked(geom, bas, 'ga_ecovl')
      call ga_zero(g_ecov)
      call int_1e_ga(bas, bas, g_ecov, 'overlap', oskel)

      
c     ==== Printing ====
      if (ndid == 0) write(luout, '("Numerical overlap :")')
      call ga_print(g_ov)

      if (ndid == 0) write(luout, '("Analytical overlap :")')
      call ga_print(g_ecov)
      
      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_ov0))
     &     call errquit('WFAT: Cannot deallocate overlap_z', 1, MA_ERR)
      if (.not. ma_free_heap(l_ov))
     &     call errquit('WFAT: Cannot deallocate overlap', 1, MA_ERR)
      if (.not. ga_destroy(g_ov))
     &     call errquit('WFAT: Cannot destroy ga_ovl', 1, MA_ERR)
      if (.not. ga_destroy(g_ecov))
     &     call errquit('WFAT: Cannot destroy ga_ecovl', 1, MA_ERR)

      call ga_sync()
      

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_test_numerfexch(rtdb)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "errquit.fh"
#include "geom.fh"
#include "wfat.fh"      
#include "wfat_qm_param.fh"
      integer rtdb
      
      character*16 tag
      character*4 xclvl
      logical rs, srhfx, t1
      integer l_rd, l_bv, l_bv1, l_bv2, l_cf, l_hf0, l_hf, l_xc0, l_xc 
      integer l_sc, l_pt, l_wt, l_icp, l_ipp
      integer k_rd, k_bv, k_bv1, k_bv2, k_cf, k_hf0, k_hf, k_xc0, k_xc
      integer k_sc, k_pt, k_wt, k_icp, k_ipp
      integer nbas, natom, nrd, ncf, nhf, nxc, nsc, max_npt, npt, natomp
      integer nptp, dbas, ic, icent, ipt_i, ipt_f, ipt, ipt0, is, ib_l
      integer ib_k, norb, sid(2), i_hf0, i_xc0, i_hf, i_xc
      integer gqual, ndid, imf
      integer g_veff(2), g_hf(2), g_xc(2)
      double precision chfx, alpha, beta, gamma, v_c(3), v_m(3), celf
      double precision znuc, atp(3), vhfpsi(nspin0), vxcpsi(nspin0)
      double precision dummy3
      parameter (imf=0)

      logical xc_gotxc, xc_chkgrad, xc_chktau, xc_hascam, xc_hassrhf
      logical nwxc_has_cam
      integer ga_create_atom_blocked
      double precision xc_hfexch, nwxc_wght_hfx
      external xc_gotxc, xc_chkgrad, xc_chktau, xc_hascam, xc_hassrhf
      external nwxc_has_cam, ga_create_atom_blocked, xc_hfexch
      external nwxc_wght_hfx


      call ga_sync()
      ndid = ga_nodeid()
      
      norb = nspin0
      sid = (/ 1, 2 /)


      
      
c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom_hd, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
c     ==== Get number of bases ====
      if (.not. bas_numbf(bas_hd, nbas))
     &     call errquit('WFAT: Getting nbas failed in ', 53, UERR)
      
      
      if (nbas > 30) then
         if (ndid==0)
     &        write(luout, '("ATTENTION: To run wfat_test_numerfexch '//
     &        'nbas must be less than or  equal to 30. Returning ...")')
         return
      endif



      
      do is = 1, nspin0
         g_veff(is) = ga_create_atom_blocked(geom_hd, bas_hd, 'vef_mat')
      enddo
cc      call wfat_get_dft_veffmat(rtdb, geom_hd, bas_hd, nbas, nspin0,
cc     &     nel, g_mov, tol2e, oskel, ixcfit, iolgc, g_veff)

      call wfat_test_fock(rtdb, geom_hd, bas_hd, nbas, nspin0, nel, 
     &     g_mov, charge, multip, tol2e, oskel, ixcfit, iolgc)
      
      
      

c     ==== Local DFT functional level ====
      if (ndid==0) write(luout, '("herek1")')
      dbas = 0
      xclvl = 'none'
      if (xc_gotxc()) then
         dbas = 0
         xclvl = 'lda'
      endif
      if (xc_chkgrad()) then
         dbas = 2
         xclvl = 'gga'
      endif      
      if (xc_chktau()) then
         dbas = 2
         xclvl = 'mgga'
      endif
      if (xc_hascam()) then
cOLD         call errquit('WFAT: At the moment, WFAT cannot be used with '//
cOLD     &                'CAM functionals.', 1, UERR)
      endif
      if (xc_hassrhf()) then
         call errquit('WFAT: At the moment, WFAT cannot be used with '//
     &                'short-range Hartree-Fock functionals.', 1, UERR)
      endif


c     ==== Exact exchange parameters ====
      if (ndid==0) write(luout, '("herek2")')
      chfx = nwxc_wght_hfx()
      rs = nwxc_has_cam()
      if (rs) then
         call nwxc_get_cam(alpha, beta, gamma, srhfx)
         if (srhfx) call errquit('WFAT: Use of range-separated '//
     &        'exchange in WFAT at the moment does not allow '//
     &        'short range HF.', 1, UERR)
      else
         alpha = 1.0d0
         beta = 0.0d0
         gamma = 0.0d0
         srhfx = .false.
      endif
      if (ndid==0) then
         write(luout, '("a,b,c = ", 3es14.6)') alpha, beta, gamma
      endif
      

c     ==== Get the reduced density matrix ====
      nrd = nbas * nbas * nspin0
      if (.not. MA_alloc_get(MT_DBL, nrd, 'rdm', l_rd, k_rd))
     &     call errquit('WFAT: rdm allocation failed', 59, MA_ERR)
      call wfat_get_rdm(dbl_mb(k_rd:k_rd+nrd-1))


c     ==== Allocate vectors and matrrices ====
      if (.not. ma_alloc_get(mt_dbl, nbas, 'bas', l_bv, k_bv))
     &     call errquit('WFAT: Cannot allocate bas', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, 3*nbas, 'dbas1', l_bv1, k_bv1))
     &     call errquit('WFAT: Cannot allocate dbas1', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, 6*nbas, 'dbas2', l_bv2, k_bv2))
     &     call errquit('WFAT: Cannot allocate dbas2', 1, MA_ERR)
      ncf = nspin0 * nbas
      if (.not. ma_alloc_get(mt_dbl, ncf, 'vcoef', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate vcoef', 1, MA_ERR)
      nhf = nspin0 * nbas**2
      if (.not. ma_alloc_get(mt_dbl, nhf, 'hf_mat', l_hf, k_hf))
     &     call errquit('WFAT: Cannot allocate hf_mat', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nhf, 'hf_mat', l_hf0, k_hf0))
     &     call errquit('WFAT: Cannot allocate hf_mat0', 1, MA_ERR)
      nxc = nspin0 * nbas**2
      if (.not. ma_alloc_get(mt_dbl, nxc, 'xc_mat', l_xc, k_xc))
     &     call errquit('WFAT: Cannot allocate xc_mat', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nxc, 'xc_mat', l_xc0, k_xc0))
     &     call errquit('WFAT: Cannot allocate xc_mat0', 1, MA_ERR)
      nsc = nbas * (3*nbas + 1 + nspin0)
      if (.not. ma_alloc_get(mt_dbl, nsc, 'scr', l_sc, k_sc))
     &     call errquit('WFAT: Cannot allocate scracth', 1, MA_ERR)
      do is = 1, nspin0
         g_hf(is) = ga_create_atom_blocked(geom_hd, bas_hd, 'hf_mat')
         g_xc(is) = ga_create_atom_blocked(geom_hd, bas_hd, 'xc_mat')
         call ga_zero(g_hf(is))
         call ga_zero(g_xc(is))
      enddo
      

c     ==== Allocate grid-related arrays ====
      if (.not. rtdb_get(rtdb, 'wfat:gqual', mt_int, 1, gqual))
     &     call errquit('WFAT: Failed reading gqual from rtdb',
     &     1, INPUT_ERR)
      
      call wfat_get_gridmaxpt(gqual, max_npt)
c     == Quadrature points ==
      if (.not. MA_alloc_get(MT_DBL, 3*max_npt, 'gpt', l_pt, k_pt))
     &     call errquit('WFAT: gpt allocation failed', 59, MA_ERR)
c     == Quadrature weights ==
      if (.not. MA_alloc_get(MT_DBL, max_npt, 'gwt', l_wt, k_wt))
     &     call errquit('WFAT: gwt allocation failed', 59, MA_ERR)
c     == Process-dependent grid range ==      
      if (.not. ma_alloc_get(MT_INT, natom, 'icent_p', l_icp, k_icp))
     &     call errquit('WFAT: icent_p allocation failed.', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_INT, 2*natom, 'ipt_p', l_ipp, k_ipp))
     &     call errquit('WFAT: ipt_p allocation failed.', 1, MA_ERR)
      call wfat_get_gridrange(geom_hd, ndid, gqual, int_mb(k_icp),
     &                        int_mb(k_ipp), natomp, nptp)
      
      
c     ==== BEGIN QUADRATURE SUMMATION ====
c     ==== Loop over atoms ====
      if (ndid==0) write(luout, '("herek3")')
      call ga_sync()
      dbl_mb(k_hf:k_hf+nhf-1) = 0.0d0
      dbl_mb(k_xc:k_xc+nxc-1) = 0.0d0
      do ic = 1, natomp
         icent = int_mb(k_icp + ic - 1)
         if (.not. geom_cent_get(geom_hd, icent, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)

         call wfat_get_gridpt(int(znuc), gqual, max_npt, npt,
     &                        dbl_mb(k_pt : k_pt+3*max_npt-1),
     &                        dbl_mb(k_wt : k_wt+max_npt-1))

         ipt_i = int_mb(k_ipp + (ic-1)*2)
         ipt_f = int_mb(k_ipp + 1 + (ic-1)*2)
c        ==== Loop over all quadrature points around the current atom ====
c         if (ndid==0) write(luout, '("herek4")')
         do ipt = ipt_i, ipt_f
cc            if (ndid==0) write(luout, '(4i6)') ipt, ipt_i, ipt_f, icent
            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            call geom_cent2cent(geom_hd, icent, v_c, imf, v_m, t1)
cc            write(luout, '("        v_m = ", 3e14.5)') v_m

            
c           ==== Evaluate basis at v_m ====
c            if (ndid==0) write(luout, '("herek5")')
            call wfat_eval_basis_npt(geom_hd, bas_hd, dbas, 1, v_m,
     &           dbl_mb(k_bv), dbl_mb(k_bv1), dbl_mb(k_bv2),
     &           dummy3)

            
c           ==== Calculate (J+K)|bas_l> and Vxc|bas_l> ====
            do ib_l = 1, nbas     ! ib_l is column index
               dbl_mb(k_cf:k_cf+ncf-1) = 0.0d0
               do is = 1, nspin0
                  dbl_mb(k_cf + (is-1)*nbas + ib_l - 1) = 1.0d0
               enddo

c               if (ndid==0) write(luout, '("herek6")')
               call wfat_vhfpsi_1ptbas(bas_hd, nbas, nspin0, v_m,
     &              dbl_mb(k_bv), chfx, rs, alpha, beta, gamma, norb, 
     &              dbl_mb(k_cf), sid, dbl_mb(k_rd), dbl_mb(k_sc),
     &              vhfpsi)

c               if (ndid==0) write(luout, '("herek7")')
               call wfat_vxcpsi_1ptbas(nbas, nspin0, dbl_mb(k_bv), 
     &              dbl_mb(k_bv1), dbl_mb(k_bv2), norb, dbl_mb(k_cf), 
     &              sid, dbl_mb(k_rd), xclvl, vxcpsi)
               
               do is = 1, nspin0
                  do ib_k = 1, nbas    ! ib_k is row index

c                     if (ndid==0) write(luout, '(2i6)') ib_k, ib_l
                     i_hf0 = k_hf0 + (is-1)*nbas**2 + (ib_l-1)*nbas +
     &                       ib_k - 1
                     dbl_mb(i_hf0) = dbl_mb(k_bv+ib_k-1) * vhfpsi(is)
                     i_xc0 = k_xc0 + (is-1)*nbas**2 + (ib_l-1)*nbas +
     &                       ib_k - 1
                     dbl_mb(i_xc0) = dbl_mb(k_bv+ib_k-1) * vxcpsi(is)
                  enddo
               enddo
            enddo
            
            
c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom_hd, icent, v_c, wf_msm, celf)
            if (abs(celf) > 1.0d-12) then
               dbl_mb(k_hf:k_hf+nhf-1) = dbl_mb(k_hf:k_hf+nhf-1) +
     &                                   dbl_mb(k_hf0:k_hf0+nhf-1) *
     &                                   celf * dbl_mb(k_wt+ipt-1)

               dbl_mb(k_xc:k_xc+nxc-1) = dbl_mb(k_xc:k_xc+nxc-1) +
     &                                   dbl_mb(k_xc0:k_xc0+nxc-1) *
     &                                   celf * dbl_mb(k_wt+ipt-1)
            endif

         enddo
      enddo


c     ==== Accumulate to g_hf and g_xc ====
      do is = 1, nspin0
         i_hf = k_hf + (is-1)*nbas**2
         call ga_acc(g_hf(is), 1, nbas, 1, nbas, dbl_mb(i_hf),
     &        nbas, 1.0d0)
         call ga_sync()

         i_xc = k_xc + (is-1)*nbas**2
         call ga_acc(g_xc(is), 1, nbas, 1, nbas, dbl_mb(i_xc),
     &        nbas, 1.0d0)
         call ga_sync()
      enddo
      

c     ==== Printing ====
      do is = 1, nspin0
         if (ndid==0)
     &        write(luout, '("HF e-e matrix, is = ", i1)') is
         call ga_print(g_hf(is))
      enddo

      do is = 1, nspin0
         if (ndid==0)
     &        write(luout, '("XC potential matrix, is = ", i1)') is
         call ga_print(g_xc(is))
      enddo      
      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_rd))
     &     call errquit('WFAT: Cannot deallocate rdm', 1, MA_ERR)
      if (.not. ma_free_heap(l_bv))
     &     call errquit('WFAT: Cannot deallocate bas', 1, MA_ERR)
      if (.not. ma_free_heap(l_bv2))
     &     call errquit('WFAT: Cannot deallocate dbas2', 1, MA_ERR)
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Cannot deallocate vcoef', 1, MA_ERR)
      if (.not. ma_free_heap(l_hf))
     &     call errquit('WFAT: Cannot deallocate hf_mat', 1, MA_ERR)
      if (.not. ma_free_heap(l_hf0))
     &     call errquit('WFAT: Cannot deallocate hf_mat0', 1, MA_ERR)
      if (.not. ma_free_heap(l_xc))
     &     call errquit('WFAT: Cannot deallocate xc_mat', 1, MA_ERR)
      if (.not. ma_free_heap(l_xc0))
     &     call errquit('WFAT: Cannot deallocate xc_mat0', 1, MA_ERR)
      if (.not. ma_free_heap(l_sc))
     &     call errquit('WFAT: Cannot deallocate scracth', 1, MA_ERR)
      do is = 1, nspin0
         if (.not. ga_destroy(g_hf(is)))
     &        call errquit('WFAT: Cannot destroy hf_mat', 1, GA_ERR)
         if (.not. ga_destroy(g_xc(is)))
     &        call errquit('WFAT: Cannot destroy xc_mat', 1, GA_ERR)
      enddo
      if (.not. MA_free_heap(l_pt))
     &     call errquit('WFAT: gpt deallocation failed', 59, MA_ERR)
      if (.not. MA_free_heap(l_wt))
     &     call errquit('WFAT: gwt deallocation failed', 59, MA_ERR)
      if (.not. ma_free_heap(l_icp))
     &     call errquit('WFAT: icent_p deallocation failed.', 1, MA_ERR)
      if (.not. ma_free_heap(l_ipp))
     &     call errquit('WFAT: ipt_p deallocation failed.', 1, MA_ERR)
      
      end
c=======================================================
      
