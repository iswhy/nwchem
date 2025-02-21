c=======================================================
c=======================================================
      subroutine wfat_get_fock(rtdb, geom, bas, nb, ns, ne, g_mo,
     &           chg, mult, tol2e, oskel, ixcfit, iolgc, doxc, g_fk)

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
      logical doxc            ! [input], whether to include DFT XC matrix (.true.) or not (.false.).
      integer g_fk(ns)        ! [output]
      
      integer is


      do is = 1, ns
         call ga_zero(g_fk(is))
      enddo


c     ==== Calculate kinetic energy matrix ====
      call int_1e_ga(bas, bas, g_fk(1), 'kinetic', oskel)
      if (ns == 2) call ga_copy(g_fk(1), g_fk(2))


c     ==== Calculate potentials ====
      call wfat_get_dft_veffmat(rtdb, geom, bas, nb, ns, ne, g_mo, 
     &     chg, mult, tol2e, oskel,  ixcfit, iolgc, doxc, g_fk)
      
            
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_dft_veffmat(rtdb, geom, bas, nb, ns, ne, g_mo,
     &           chg, mult, tol2e, oskel, ixcfit, iolgc, doxc, g_veff)

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
      logical doxc            ! [input], whether to include DFT XC matrix (.true.) or not (.false.).
      integer g_veff(ns)      ! [in/output]
      
      character*9  vxcname(2)
      logical atmblk
      integer g_trd(2)        ! Total RDM, in practice only g_trd(1) is ever created and used.
      integer g_rd(2), g_k(2), g_j(1), g_vxc(4)
      integer g_xcinv, g_rd_fde(2)
      integer is, ndid
      double precision zero, one, mone, jfac(2), kfac(2), chfx
      double precision numnel, exc(2), normm

      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0)

      logical nwxc_has_cam
      integer ga_create_atom_blocked
      double precision nwxc_wght_hfx
      external nwxc_wght_hfx, nwxc_has_cam, ga_create_atom_blocked


      ndid = ga_nodeid()
      vxcname = (/ 'vxc_alpha', 'vxc_beta' /)


cUSELESS      if (.not. rtdb_put(rtdb, 'charge', mt_dbl, 1, chg))
cUSELESS     &     call errquit('wfat_get_dft_veffmat: Failed writing chg to '//
cUSELESS     &     'rtdb', 1, INPUT_ERR)
cUSELESS      if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, mult))
cUSELESS     &     call errquit('wfat_get_dft_veffmat: Failed writing mult to'//
cUSELESS     &     ' rtdb', 1, INPUT_ERR)
      

c     ==== Calculate nuclear potential matrix ====
      do is = 1, ns
         call int_1e_ga(bas, bas, g_veff(is), 'potential', oskel)
      enddo
      

c     ==== Calculate RDM and total RDM ====
      do is = 1, ns
         g_rd(is) = ga_create_atom_blocked(geom, bas, 'rdm')
         call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     &        g_mo(is), 1, nb, 1, ne(is),
     &        g_mo(is), 1, ne(is), 1, nb,
     &        g_rd(is), 1, nb, 1, nb)
      enddo
      g_trd(1) = ga_create_atom_blocked(geom, bas, 'tot_rdm')
      call ga_zero(g_trd(1))
      if (ns == 1) then
         call ga_copy(g_rd(1), g_trd(1))
         call ga_scale(g_trd(1), 2.0d0)
      else if (ns == 2) then
         call ga_dadd(one, g_rd(1), one, g_rd(2), g_trd(1))
      endif


c     ==== Calculate Hartree (classical repulsion) matrix ====
      g_j(1) = ga_create_atom_blocked(geom, bas, 'hartree_mat')
      call ga_zero(g_j(1))
      jfac(1) = one
      kfac(1) = zero
      call fock_2e(geom, bas, 1, jfac, kfac, tol2e, oskel,
     &     g_trd, g_j, .false.)
cfordebug      if (ndid==0) write(luout, '("HF Hartree matrix (REF)")')
cfordebug      call ga_print(g_j(1))
      

c     ==== Calculate HF exchange matrix ====      
      do is = 1, ns
         g_k(is) = ga_create_atom_blocked(geom, bas, 'hfexch_mat')
         call ga_zero(g_k(is))
      enddo
      jfac = zero
      kfac = one
      
      if (nwxc_has_cam()) call case_setflags (.true.)
      call fock_2e(geom, bas, ns, jfac, kfac, tol2e, oskel,
     &     g_rd, g_k, .false.)
cfordebug      if (ndid==0) write(luout, '("HF exchange matrix (REF), is=1")')
cfordebug      call ga_print(g_k(1))      
      if (nwxc_has_cam()) call case_setflags (.false.)


c     ==== Calculate HF 2e matrix (J+c*K) ====
      chfx = nwxc_wght_hfx()
      do is  = 1, ns
         call ga_dadd(one, g_j, mone*chfx, g_k(is), g_k(is))
cfordebug         if (ndid==0) write(luout, '("HF e-e matrix (REF), is=1")')
cfordebug         call ga_print(g_k(1))
      enddo
      
      
c     ==== Obtain Vxc matrix ====
      if (doxc) then
         if (.not. rtdb_get(rtdb, 'dft:atomblock', mt_log, 1, atmblk))
     &        atmblk = .true.
         do is = 1, ns
            if (atmblk) then
               g_vxc(is) = ga_create_atom_blocked(geom, bas,
     &                     vxcname(is))
            else
               if (.not. ga_create(mt_dbl, nb, nb, vxcname(is),
     &              0, 0, g_vxc(is)))
     &              call errquit('WFAT: Error creating '//vxcname(is)//
     &              ' global array.', 0, GA_ERR)
            endif
            call ga_zero(g_vxc(is))
         enddo
         
         if (ns == 1) then
            call xc_getvxc(rtdb, exc, 2, ixcfit, g_xcinv, g_trd, g_vxc,
     &           iolgc, numnel, 0, g_rd_fde)
cfordebug         if (ndid==0)
cfordebug     &        write(luout, '("XC potential matrix (REF), is=1")')
cfordebug         call ga_print(g_vxc(1))
         else if (ns == 2) then
            call xc_getvxc(rtdb, exc, 2, ixcfit, g_xcinv, g_rd, g_vxc,
     &           iolgc, numnel, 0, g_rd_fde)
         endif
      endif
      

c     ==== Calculate (total) Fock matrices ====
      do is = 1, ns
         call ga_dadd(one, g_veff(is), one, g_k(is), g_k(is))
         if (doxc) then
            call ga_dadd(one, g_k(is), one, g_vxc(is), g_veff(is))
         else
            call ga_copy(g_k(is), g_veff(is))
         endif
      enddo
      

c     ==== Deallocations ====
      
      do is = 1, ns
         if (.not. ga_destroy(g_k(is)))
     &        call errquit('wfat_test_hamiltonian: Cannot destroy '//
     &        'hfexch_mat', 41 , GA_ERR)
         if (.not. ga_destroy(g_rd(is)))
     &        call errquit('wfat_test_hamiltonian: Cannot destroy '//
     &        'rdm.', 1 , GA_ERR)
         if (doxc) then
            if (.not. ga_destroy(g_vxc(is)))
     &           call errquit('wfat_test_hamiltonian: Cannot destroy '//
     &           vxcname(is)//'.', 1 , GA_ERR)
         endif
      enddo      
      if (.not. ga_destroy(g_trd(1)))
     &     call errquit('wfat_test_hamiltonian: Cannot destroy '//
     &     'tot_rdm', 1 , GA_ERR)
      if (.not. ga_destroy(g_j(1)))
     &     call errquit('wfat_test_hamiltonian: Cannot destroy '//
     &     'hartree_mat', 1 , GA_ERR)
      

      end
c=======================================================
