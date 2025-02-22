c==========================================================
c==========================================================
      subroutine wfat_2a1efd_target(geom, mch, sphi, orbid, pottype, 
     &     nel, potparam, ni, n, maxn, xq, wq, nr, nt, g_cf, erg,
     &     dipole)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "global.fh"
#include "errquit.fh"
      integer geom
      integer mch
      integer sphi
      integer orbid
      character*3 pottype
      integer nel                    ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)
      integer ni
      integer n(ni)
      integer maxn
      double precision xq(maxn,ni)
      double precision wq(maxn,ni)
      integer nr
      integer nt
      integer g_cf                  ! [input], GA of size nr*nt x nr*nt
      double precision erg(nr*nt)
      double precision dipole(3)

      character*5 sm
      integer l_dt, l_x, k_dt, k_x, l_cf, k_cf, ndid

      
      ndid = ga_nodeid()

      
c     ==== Solve the eigenproblem for GS ====
      call wfat_2a1efd_solve(mch, ni, n, maxn, xq, wq, nt, geom, 
     &     pottype, nel, potparam, g_cf, erg)


c     ==== Electronic dipole moment ====
      dipole = 0.0d0
      if (.not. ma_alloc_get(MT_DBL, (nt*nt)**2, 'dt', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt', 26, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad_', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad_', 86, MA_ERR)
      call wfat_fedvr_basnodes(ni, n, maxn, xq, dbl_mb(k_x))
      
      if (.not. ma_alloc_get(mt_dbl, nr*nt, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Could not allocate coef.', 56, MA_ERR)
      call ga_get(g_cf, 1, nr*nt, orbid, orbid, dbl_mb(k_cf), nr*nt)
      call wfat_2a1efd_dpo('z', mch, 0, nr, nt, nt, sphi, sphi,
     &     dbl_mb(k_cf), dbl_mb(k_cf), dbl_mb(k_x), .false., 
     &     dbl_mb(k_dt), .false., dipole(3))
      
      
c     ==== Analysis ====
      if (ndid == 0) then
         write(sm, '(i5)') mch
         write(luout, '(2/)')
         call util_print_centered(luout,
     &        'Eigenvalue spectrum of m = '//trim(adjustl(sm))//
     &        ' manifold', 60, .true.)
         call wfat_2a1e_print_eigvals(nr*nt, 20, erg)
      endif
      
      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt', 1, MA_ERR)
      if (.not. ma_free_heap(l_x))
     &     call errquit('WFAT: Cannot deallocate xquad_', 1, MA_ERR)
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Could not deallocate coef', 11, MA_ERR)

      end
c==========================================================
      
