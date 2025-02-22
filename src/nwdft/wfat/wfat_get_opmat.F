c=======================================================
c=======================================================      
      subroutine wfat_get_vhfmat(geom, bas, nb, ns, g_rd,
     &           chfx, tol2e, oskel, vhfmat)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer geom                  ! [input], geometry handle.                                 
      integer bas                   ! [input], basis set handle.                                
      integer nb                    ! [input], no. of bases associated with the argument 'bas'. 
      integer ns                    ! [input], no. of spin channels in the system.              
      integer g_rd(ns)              ! [input], GA handle to an RDM matrix. If ns=1, then g_rd(1) is the alpha RDM.
      double precision chfx         ! [input], prefactor for exact exchange.
      double precision tol2e        ! [input], tolerance for elements of 2e integrals.
      logical oskel                 ! [input], toggle skeleton Fock matrix.
      double precision vhfmat(nb,nb,ns) ! [output], coefficients of vhf|vin>.

      character*9 hfx_name(2)
      integer g_trd(1), g_j(1), g_k(ns)
      integer i, is
      double precision zero, one, mone, jfac(1), kfac(1)

      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0)


      hfx_name = (/ 'hfx_alpha', 'hfx_beta' /)
      


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
      jfac(1) = zero
      kfac(1) = one
      do is = 1, ns
         if (.not. ga_create(mt_dbl, nb, nb, hfx_name(is), 0, 0,
     &        g_k(is)))
     &        call errquit('WFAT: Cannot create '//hfx_name(is)//'.',
     &        1, GA_ERR)
         call ga_zero(g_k(is))
         call fock_2e(geom, bas, 1, jfac, kfac, tol2e, oskel,
     &        g_rd(is), g_k(is), .false.)
         call ga_dadd(one, g_j, mone*chfx, g_k(is), g_k(is))
      enddo


c     ===== Copy V_ee matrix to local buffer ====
      do is = 1, ns
         call ga_get(g_k(is), 1, nb, 1, nb, vhfmat(1,1,is), nb)
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


      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_vxcmat(rtdb, geom, bas, nb, ns, g_rd,
     &           ixcfit, iolgc, numnel, exc, vxcmat)
      

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
      integer ixcfit                ! [input], if 1 then do density fitting to calculate exchange potential.
      logical iolgc                 ! [input], if .true. then does not use disk for exchange fitting, otherwise, use disk.
      double precision numnel       ! [output], numerically calculated no. of electrons.
      double precision exc(2)       ! [output], XC energies.                            
      double precision vxcmat(nb,nb,ns)  ! [output], coefficients of vxc|vin>.

      character*9 rdmname(2), vxcname(2)
      logical atmblk
      integer i, is
      integer g_rd0(2)
      integer g_xcinv           ! This GA must not be created here because will be in a deeper routine.
      integer g_vxc(4)          ! This GA must be created here, i.e. not in deeper routines.
      integer g_rd_fde(2)       ! Is unused (should never be used for the present routine) as long as the 10th argument of subroutine xc_getvxc below is equal to zero. 

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      
      rdmname = (/ 'rdm_alpha', 'rdm_beta' /)
      vxcname = (/ 'vxc_alpha', 'vxc_beta' /)
      
      
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
      do is = 1, ns
         call ga_get(g_vxc(is), 1, nb, 1, nb, vxcmat(1,1,is), nb)
      enddo
      

c     ==== Deallocations ====
      do is = 1, ns
         if (.not. ga_destroy(g_vxc(is)))
     &        call errquit('WFAT: Cannot destroy Vxc.', 1, GA_ERR)
         if (.not. ga_destroy(g_rd0(is)))
     &        call errquit('WFAT: Cannot destroy rdm0.', 1, GA_ERR)         
      enddo
      
      end
c=======================================================      
