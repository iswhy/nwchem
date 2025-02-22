c=======================================================
c=======================================================
      subroutine wfat_dft_veffvec_coef(rtdb, geom, bas, nb, ns, g_rd, 
     &           nv, vin, sid, tol2e, oskel, ixcfit, iolgc,
     &           numnel, exc, vout)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer rtdb                    ! [input], RTDB handle.                                     
      integer geom                    ! [input], geometry handle.                                 
      integer bas                     ! [input], basis set handle.                                
      integer nb                      ! [input], no. of bases associated with the argument 'bas'. 
      integer ns                      ! [input], no. of spin channels in the system.              
      integer g_rd(2)                 ! [input], GA handle to an RDM matrix. If ns=1, then g_rd(1) is the alpha RDM while g_rd(2) is unused.
      integer nv                      ! [input], no. of vectors in 'vin'.                          
      double precision vin(nb, nv)    ! [input], input vectors to be acted on by Vxc matrix.       
      integer sid(nv)                 ! [input], spin id (i.e. '1' or '2') of the vectors in 'vin'.
      double precision tol2e          ! [input], tolerance for elements of 2e integrals.
      logical oskel                   ! [input], toggle skeleton Fock matrix.
      integer ixcfit                  ! [input], if 1 then do density fitting to calculate exchange potential.
      logical iolgc                   ! [input], if .true. then does not use disk for exchange fitting, otherwise, use disk.
      double precision numnel         ! [output], numerically calculated no. of electrons.
      double precision exc(2)         ! [output], XC energies.                            
      double precision vout(nb, nv)   ! [output], coefficients of veff|vin>.

      integer l_vnv, k_vnv, i_vnv, nvnv, i
      

c     ==== Nuclear attraction potential contribution ====
      nvnv = nb * nv
      if (.not. ma_alloc_get(mt_dbl, nvnv, 'vnucvin', l_vnv, k_vnv))
     &     call errquit('WFAT: Cannot allocate vnucvin.', 1, MA_ERR)
      call wfat_vnucvec_coef(bas, nb, nv, vin, dbl_mb(k_vnv))


c     ==== Electron-electron interaction potential contribution ====
      call wfat_dft_veevec_coef(rtdb, geom, bas, nb, ns, g_rd, nv, vin,
     &                          sid, tol2e, oskel, ixcfit, iolgc,
     &                          numnel, exc, vout)


c     ==== Calculate vout = veff|vin> ====
      do i = 1, nv
         i_vnv = k_vnv + (i-1)*nb
         vout(:,i) = dbl_mb(i_vnv:i_vnv+nb-1) + vout(:,i)
      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_vnv))
     &     call errquit('WFAT: Cannot deallocate vnucvin.', 1, MA_ERR)
      
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_vnucvec_coef(bas, nb, nv, vin, vout)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"      
      integer bas
      integer nb
      integer nv
      double precision vin(nb,nv)
      double precision vout(nb,nv)

      integer l_pe, k_pe
      integer npe


      npe = nb * nb
      if (.not. ma_alloc_get(MT_DBL, npe, 'vnuc', l_pe, k_pe))
     &     call errquit('WFAT: Cannot allocate vnuc in '//
     &     'wfat_vnucvec_coef.', 1, MA_ERR)
      
      call int_1epe_all(bas, nb, dbl_mb(k_pe))
      call wfat_sopvec_coef(bas, nb, dbl_mb(k_pe), nv, vin, vout)

      if (.not. ma_free_heap(l_pe))
     &     call errquit('WFAT: Cannot deallocate vnuc in '//
     &     'wfat_vnucvec_coef.', 1, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dft_veevec_coef(rtdb, geom, bas, nb, ns, g_rd, nv, 
     &     vin, sid, tol2e, oskel, ixcfit, iolgc,
     &     numnel, exc, vout)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"      
      integer rtdb                    ! [input], RTDB handle.                                     
      integer geom                    ! [input], geometry handle.                                 
      integer bas                     ! [input], basis set handle.                                
      integer nb                      ! [input], no. of bases associated with the argument 'bas'. 
      integer ns                      ! [input], no. of spin channels in the system.              
      integer g_rd(2)                 ! [input], GA handle to an RDM matrix. If ns=1, then g_rd(1) is the alpha RDM while g_rd(2) is unused.
      integer nv                      ! [input], no. of vectors in 'vin'.                          
      double precision vin(nb, nv)    ! [input], input vectors to be acted on by Vxc matrix.       
      integer sid(nv)                 ! [input], spin id (i.e. '1' or '2') of the vectors in 'vin'.
      double precision tol2e          ! [input], tolerance for elements of 2e integrals.
      logical oskel                   ! [input], toggle skeleton Fock matrix.
      integer ixcfit                  ! [input], if 1 then do density fitting to calculate exchange potential.
      logical iolgc                   ! [input], if .true. then does not use disk for exchange fitting, otherwise, use disk.
      double precision numnel         ! [output], numerically calculated no. of electrons.
      double precision exc(2)         ! [output], XC energies.                            
      double precision vout(nb, nv)   ! [output], coefficients of vee|vin>.

      integer i, nhf, l_hf, k_hf, i_hf
      double precision chfx

      double precision xc_hfexch
      external xc_hfexch


      
      if (ns /= 1 .and. ns /= 2) then
         call errquit('WFAT: Argument ''ns'' of wfat_dft_veevec_coef '//
     &        'must be either 1 or 2.', 1, UERR)
      endif
      

c     ==== Hartree-Fock potential (J + K) contribution ====
      chfx = xc_hfexch()
      nhf = nb * nv
      if (.not. ma_alloc_get(MT_DBL, nhf, 'vhfvec', l_hf, k_hf))
     &     call errquit('WFAT: Cannot allocate vhfvec.', 1, MA_ERR)
      call wfat_vhfvec_coef(geom, bas, nb, ns, g_rd(1:ns), nv, vin,
     &     sid, chfx, tol2e, oskel, dbl_mb(k_hf))


c     ==== DFT XC potential contribution ====
      call wfat_vxcvec_coef(rtdb, geom, bas, nb, ns, g_rd(1:2), nv, vin,
     &     sid, ixcfit, iolgc, numnel, exc, vout)


c     ==== Calculate vout = vee|vin> ====
      do i = 1, nv
         i_hf = k_hf + (i-1)*nb
         vout(:,i) = dbl_mb(i_hf:i_hf+nb-1) + vout(:,i)
      enddo

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_hf))
     &     call errquit('WFAT: Cannot deallocate vhfvec.', 1, MA_ERR)
      
      end
c=======================================================

      
c=======================================================
c=======================================================      
      subroutine wfat_vhfvec_coef(geom, bas, nb, ns, g_rd, nv, vin,
     &           sid, chfx, tol2e, oskel, vout)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"      
      integer geom                  ! [input], geometry handle.                                 
      integer bas                   ! [input], basis set handle.                                
      integer nb                    ! [input], no. of bases associated with the argument 'bas'. 
      integer ns                    ! [input], no. of spin channels in the system.              
      integer g_rd(ns)              ! [input], GA handle to an RDM matrix. If ns=1, then g_rd(1) is the alpha RDM.
      integer nv                    ! [input], no. of vectors in 'vin'.                          
      double precision vin(nb, nv)  ! [input], input vectors to be acted on by Vhf matrix.       
      integer sid(nv)               ! [input], spin id (i.e. '1' or '2') of the vectors in 'vin'.
      double precision chfx         ! [input], prefactor for exact exchange.
      double precision tol2e        ! [input], tolerance for elements of 2e integrals.
      logical oskel                 ! [input], toggle skeleton Fock matrix.
      double precision vout(nb, nv) ! [output], coefficients of vhf|vin>.

      character*9 hfx_name(2)
      logical has_a, has_b
      integer g_trd(1), g_j(1), g_k(ns)
      integer l_ve, k_ve, i_ve
      integer i, is, nve
      double precision zero, one, mone, jfac(1), kfac(1)

      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0)


      hfx_name = (/ 'hfx_alpha', 'hfx_beta' /)
      


c     ==== Determine which spin channels are in sid ====
      has_a = .false.
      has_b = .false.
      do i = 1, nv
         if (sid(i) == 1) then
            has_a = .true.
         else if (sid(i) == 2) then
            has_b = .true.
         else if (sid(i) == 2 .and. ns == 1) then
            call errquit('WFAT: An element of input array ''sid'' '//
     &           'of subroutine wfat_vhfvec_coef is equal to 2 while '//
     &           'the input ''ns'' is equal to 1.', 1, UERR)
         else
            call errquit('WFAT: Elements of input array ''sid'' of '//
     &           'subroutine wfat_vhfvec_coef must be either 1 or 2.',
     &           1, UERR)
         endif
      enddo

      
c     ==== Calculate total RDM ====
      if (.not. ga_create(mt_dbl, nb, nb, 'tot_rdm', 0, 0, g_trd(1)))
     &     call errquit('WFAT: Cannot create ''tot_rdm'' global array.',
     &     1, GA_ERR)
      call ga_zero(g_trd(1))
      if (ns == 1) then
         call ga_copy(g_rd(1), g_trd(1))
         call ga_scale(g_trd(1), 2.0d0)
      else if (ns == 2) then
         call ga_dadd(one, g_rd(1), one, g_rd(2), g_trd(1))
      endif
      

c     ==== Calculate Hartree (classical repulsion) matrix ====      
      if (.not. ga_create(mt_dbl, nb, nb, 'hartree_mat', 0, 0,
     &     g_j(1)))
     &     call errquit('WFAT: Cannot create ''hartree_mat'' global '//
     &     'array.', 1, GA_ERR)
      call ga_zero(g_j(1))
      jfac(1) = one
      kfac(1) = zero
      call fock_2e(geom, bas, 1, jfac, kfac, tol2e, oskel,
     &     g_trd, g_j, .false.)

      
c     ==== Calculate HF exchange matrix and V_ee matrix ====
      do is = 1, ns
         if (.not. ga_create(mt_dbl, nb, nb, hfx_name(is), 0, 0,
     &        g_k(is)))
     &        call errquit('WFAT: Cannot create '//hfx_name(is)//'.',
     &        1, GA_ERR)
         call ga_zero(g_k(is))
      enddo

      jfac(1) = zero
      kfac(1) = one
      if (has_a) then
         call fock_2e(geom, bas, 1, jfac, kfac, tol2e, oskel,
     &        g_rd(1), g_k(1), .false.)
         call ga_dadd(one, g_j, mone*chfx, g_k(1), g_k(1))
c         write(luout, '("HF J+K matrix, spin:", i2)') 1
c         call ga_print(g_k(1))
      endif
      if (has_b) then
         call fock_2e(geom, bas, 1, jfac, kfac, tol2e, oskel,
     &        g_rd(2), g_k(2), .false.)
         call ga_dadd(one, g_j, mone*chfx, g_k(2), g_k(2))
      endif


c     ===== Copy V_ee matrix to local buffer ====
      nve = nb * nb * ns
      if (.not. ma_alloc_get(mt_dbl, nve, 'HF_vee', l_ve, k_ve))
     &     call errquit('WFAT: Cannot allocate HF_vee.', 1, MA_ERR)
      do is = 1, ns
         i_ve = k_ve + (is-1)*nb**2
         call ga_get(g_k(is), 1, nb, 1, nb, dbl_mb(i_ve), nb)
c         write(luout, '("HF Vee matrix, spin:", i2)') is
c         call wfat_dprintmat(nb, nb, 'e', dbl_mb(i_ve))
      enddo
      

c     ==== Calculate vout = (J + K)|vin> ====
      do i = 1, nv
         is = sid(i)
         i_ve = k_ve + (is-1)*nb*nb
         call wfat_sopvec_coef(bas, nb, dbl_mb(i_ve), 1, vin(:,i),
     &        vout(:,i))
      enddo

      
c     ==== Deallocations ====
      if (.not. ga_destroy(g_trd))
     &     call errquit('WFAT: Cannot destroy tot_rdm.', 1, GA_ERR)
      if (.not. ga_destroy(g_j))
     &     call errquit('WFAT: Cannot destroy hartree_mat.', 1, GA_ERR)
      do is = 1, ns
         if (.not. ga_destroy(g_k(is)))
     &        call errquit('WFAT: Cannot destroy hfx.', 1, GA_ERR)
      enddo
      if (.not. ma_free_heap(l_ve))
     &     call errquit('WFAT: Cannot deallocate HF_vee.', 1, MA_ERR)
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_vxcvec_coef(rtdb, geom, bas, nb, ns, g_rd, nv,
     &           vin, sid, ixcfit, iolgc, numnel, exc, vout)
      

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"      
      integer rtdb                  ! [input], RTDB handle.
      integer geom                  ! [input], geometry handle.
      integer bas                   ! [input], basis set handle.
      integer nb                    ! [input], no. of bases associated with the argument 'bas'.
      integer ns                    ! [input], no. of spin channels in the system.              
      integer g_rd(2)               ! [input], GA handle to an RDM matrix. If ns=1, then g_rd(1) is the alpha RDM while g_rd(2) is unused.
      integer nv                    ! [input], no. of vectors in 'vin'.                          
      double precision vin(nb,nv)   ! [input], input vectors to be acted on by Vxc matrix.       
      integer sid(nv)               ! [input], spin id (i.e. '1' or '2') of the vectors in 'vin'.
      integer ixcfit                ! [input], if 1 then do density fitting to calculate exchange potential.
      logical iolgc                 ! [input], if .true. then does not use disk for exchange fitting, otherwise, use disk.
      double precision numnel       ! [output], numerically calculated no. of electrons.
      double precision exc(2)       ! [output], XC energies.                            
      double precision vout(nb,nv)  ! [output], coefficients of vxc|vin>.

      character*9 rdmname(2), vxcname(2)
      logical atmblk
      integer i, is
      integer nvxc
      integer l_vxc, k_vxc, i_vxc
      integer g_rd0(2)
      integer g_xcinv           ! This GA must not be created here because will be in a deeper routine.
      integer g_vxc(4)          ! This GA must be created here, i.e. not in deeper routines.
      integer g_rd_fde(2)       ! Is unused (should never be used for the present routine) as long as the 10th argument of subroutine xc_getvxc below is equal to zero. 

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      
      rdmname = (/ 'rdm_alpha', 'rdm_beta' /)
      vxcname = (/ 'vxc_alpha', 'vxc_beta' /)
      

c     ==== Determine which spin channels are in sid ====
      do i = 1, nv
         if (sid(i) == 2 .and. ns == 1) then
            call errquit('WFAT: An element of input array ''sid'' '//
     &           'of subroutine wfat_vxcvec_coef is equal to 2 while '//
     &           'the input ''ns'' is equal to 1.', 1, UERR)
         else if (sid(i) /= 1 .and. sid(i) /= 2) then
            call errquit('WFAT: Elements of input array ''sid'' of '//
     &           'subroutine wfat_vxcvec_coef must be either 1 or 2.',
     &           1, UERR)
         endif
      enddo

      
c     ==== Modify RDM based on ns ====
      do is = 1, ns
         if (.not. ga_duplicate(g_rd(is), g_rd0(is), rdmname(is)))
     &        call errquit('WFAT: Cannot duplicate '//rdmname(is)//'.',
     &        1, GA_ERR)
         call ga_copy(g_rd(is), g_rd0(is))
      enddo
      if (ns == 1) call ga_scale(g_rd0(1), 2.0d0)
         

c     ==== Obtain Vxc matrix ====
      if (.not. rtdb_get(rtdb, 'dft:atomblock', mt_log, 1, atmblk))
     &     atmblk = .true.
      do is = 1, ns
         if (atmblk) then
            g_vxc(is) = ga_create_atom_blocked(geom, bas, vxcname(is))
         else
            if (.not. ga_create(mt_dbl, nb, nb, vxcname(is),
     &           0, 0, g_vxc(is)))
     &           call errquit('WFAT: Error creating '//vxcname(is)//
     &           ' global array.', 0, GA_ERR)
         endif
         call ga_zero(g_vxc(is))
      enddo      
      call xc_getvxc(rtdb, exc, 2, ixcfit, g_xcinv, g_rd0, g_vxc,
     &     iolgc, numnel, 0, g_rd_fde)

      
c     ===== Copy V_ee matrix to local buffer ====
      nvxc = nb * nb * ns
      if (.not. ma_alloc_get(mt_dbl, nvxc, 'XC_vee', l_vxc, k_vxc))
     &     call errquit('WFAT: Cannot allocate XC_vee.', 1, MA_ERR)
      do is = 1, ns
         i_vxc = k_vxc + (is-1)*nb**2
         call ga_get(g_vxc(is), 1, nb, 1, nb, dbl_mb(i_vxc), nb)
c         write(luout, '("DFT XC matrix, spin:", i2)') is
c         call wfat_dprintmat(nb, nb, 'e', dbl_mb(i_vxc))
      enddo
      

c     ==== Calculate vout = Vxc|vin> ====
      do i = 1, nv
         is = sid(i)
         i_vxc = k_vxc + (is-1)*nb**2
         call wfat_sopvec_coef(bas, nb, dbl_mb(i_vxc), 1, vin(:,i),
     &        vout(:,i))
      enddo      


c     ==== Deallocations ====
      do is = 1, ns
         if (.not. ga_destroy(g_vxc(is)))
     &        call errquit('WFAT: Cannot destroy Vxc.', 1, GA_ERR)
         if (.not. ga_destroy(g_rd0(is)))
     &        call errquit('WFAT: Cannot destroy rdm0.', 1, GA_ERR)         
      enddo
      if (.not. ma_free_heap(l_vxc))
     &     call errquit('WFAT: Cannot deallocate XC_vee.', 1, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_vhydvec_coef(bas, nb, rc, cg, nv, vin, vout)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer bas
      integer nb
      double precision rc(3)
      double precision cg
      integer nv
      double precision vin(nb,nv)
      double precision vout(nb,nv)

      integer l_hy, k_hy
      integer nhy
      

      nhy = nb * nb
      if (.not. ma_alloc_get(MT_DBL, nhy, 'vhyd', l_hy, k_hy))
     &     call errquit('WFAT: Cannot allocate vhyd in '//
     &     'wfat_vhydvec_coef.', 1, MA_ERR)

      call int_1ehyd_all(bas, nb, rc, cg, dbl_mb(k_hy))
      call wfat_sopvec_coef(bas, nb, dbl_mb(k_hy), nv, vin, vout)

      if (.not. ma_free_heap(l_hy))
     &     call errquit('WFAT: Cannot deallocate vhyd in '//
     &     'wfat_vhydvec_coef.', 1, MA_ERR)
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_sopvec_coef(bas, nb, op, nv, vin, vout)
c     op = symmetric matrix
c     op|vin> = |vout>
c     [op] * [vin] = [S] * [vout]
      
      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"      
      integer bas                  ! [input]  
      integer nb                   ! [input]
      double precision op(nb,nb)   ! [input]
      integer nv                   ! [input]
      double precision vin(nb,nv)  ! [input]
      double precision vout(nb,nv) ! [output]

      integer l_ov, k_ov
      integer nov


c     ==== Do op * vin = vrhs ====
      call wfat_dsymm_wrap(nb, nv, 1.0d0, op, vin, 0.0d0, vout)

c     ==== Solve S * vout = vrhs for vout ====
      nov = nb * nb
      if (.not. ma_alloc_get(MT_DBL, nov, 'ovlap', l_ov, k_ov))
     &     call errquit('WFAT: Cannot allocate ovlap.', 1, MA_ERR)

      call int_1eov_all(bas, nb, dbl_mb(k_ov))
      call wfat_dsysv_wrap(nb, nv, dbl_mb(k_ov), vout)


      if (.not. ma_free_heap(l_ov))
     &     call errquit('WFAT: Cannot deallocate ovlap.', 1, MA_ERR)
      

      end
c=======================================================      
