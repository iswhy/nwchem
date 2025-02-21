c=======================================================
c=======================================================
      subroutine wfat_test_dtrm(rtdb, geom, bas, nb, ns, ne, g_mo,
     &     chg, mult, moerg, tol2e, oskel, ixcfit, iolgc)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb                    ! [input]
      integer geom                    ! [input]
      integer bas                     ! [input]
      integer nb                      ! [input]
      integer ns                      ! [input]
      integer ne(ns)                  ! [input]
      integer g_mo(ns)                ! [input]
      double precision chg            ! [input]
      integer mult                    ! [input]
      double precision moerg(nb, ns)  ! [input]
      double precision tol2e          ! [input]
      logical oskel                   ! [input]
      integer ixcfit                  ! [input]
      logical iolgc                   ! [input]

      character*9 hfx_name(2)
      integer is, io, ndid, g_fk(ns), k_du, l_du, g_s, g_du
      double precision dets, detv0(2), detv(2)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      

      ndid = ga_nodeid()

c     ==== Calc. determinant by eigenvalues ====
      g_s = ga_create_atom_blocked(geom, bas, 'overlap')
      call ga_zero(g_s)
      call int_1e_ga(bas, bas, g_s, 'overlap', oskel)
      call wfat_et_det(g_s, nb, dets)      
      detv0 = 1.0d0
      if (ndid==0)
     &     write(luout, '(/, "Determinant by eigenvalues")')
      do is = 1, ns
         do io = 1, nb
            if (ndid==0) write(luout, '(2i4, es14.5)')
     &                   is, io, moerg(io,is)
            detv0(is) = detv0(is) * moerg(io, is)
         enddo
         if (ndid==0) write(luout, '(i1, ") det. = ", es16.8)') is,
     &        detv0(is) * dets
      enddo
      

c     ==== Calc. determinant explicitly ====
      hfx_name = (/ 'hfx_alpha', 'hfx_beta' /)
      do is = 1, ns
         g_fk(is) = ga_create_atom_blocked(geom, bas, hfx_name(is))
         call ga_zero(g_fk(is))
      enddo
      call wfat_get_fock(rtdb, geom, bas, nb, ns, ne, g_mo, chg, mult,
     &                   tol2e, oskel, ixcfit, iolgc, .true., g_fk)
      if (ndid==0)
     &     write(luout, '(/, "Determinant by Gauss elimination")')
      do is = 1, ns
         call wfat_et_det(g_fk(is), nb, detv(is))
         if (ndid==0) write(luout, '(i1, ") det. = ", es16.8)') is,
     &        detv(is)
      enddo


c     ==== Dummy matrix test ====
      if (.not. ma_alloc_get(mt_dbl, 4, 'dummy', l_du, k_du))
     &     call errquit('WFAT: Error here 1', 1, MA_ERR)
      dbl_mb(k_du) = 35.0d0
      dbl_mb(k_du+1) = 2.0d0
      dbl_mb(k_du+2) = 2.0d0
      dbl_mb(k_du+3) = -2.0d0

      if (.not. ga_create(mt_dbl, 2, 2, 'dummy', 0, 0, g_du))
     &     call errquit('WFAT: Error here 2', 1, GA_ERR)
      call ga_put(g_du, 1, 2, 1, 2, dbl_mb(k_du), 2)
      if (ndid==0) write(luout, '("Dummy matrix:")')
      call ga_print(g_du)
      call wfat_et_det(g_du, 2, detv0(1))
c      call et_det('t', g_du, 2, detv0(2))
      if (ndid==0) then
         write(luout, '("Determinant A:", es16.8)') detv0(1)
c         write(luout, '("Determinant A^T:", es14.6)') detv0(2)
      endif

      if (.not. ma_free_heap(l_du))
     &     call errquit('WFAT: Error here 3', 1, MA_ERR)
      if (.not. ga_destroy(g_du))
     &     call errquit('WFAT: Error here 4', 1, GA_ERR)
      if (.not. ga_destroy(g_s))
     &     call errquit('WFAT: Error here 4', 1, GA_ERR)
      do is = 1, ns
         if (.not. ga_destroy(g_fk(is)))
     &        call errquit('WFAT: Error here 7', 1, GA_ERR)
      enddo
      
      
      end
c=======================================================
      
