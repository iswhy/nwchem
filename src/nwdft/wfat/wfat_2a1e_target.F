c===============================================================
c===============================================================
      subroutine wfat_2a1e_target(geom, m, sphi, target, pottype, nel, 
     &     potparam, nr, nt, init_ke, mult, maxiter, g_cf, erg, dipole,
     &     ke)

      implicit none
#include "global.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer m
      integer sphi
      integer target
      character*3 pottype
      integer nel                    ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)
      integer nr
      integer nt
      double precision init_ke
      double precision mult
      integer maxiter
      integer g_cf                   ! [input], GA of size nr*nt x nr*nt
      double precision erg(nr*nt)
      double precision dipole(3)
      double precision ke

      character*5 sm
      integer l_dr, k_dr, l_dt, k_dt, l_cf, k_cf, ndid

      
      ndid = ga_nodeid()

c     ==== Solve the eigenproblem for GS ====
      ke = init_ke
      call wfat_2a1e_solve(m, nr, nt, geom, ke, mult, .true., maxiter, 
     &     target, pottype, nel, potparam, g_cf, erg)


c     ==== Electronic dipole moment ====
      dipole = 0.0d0
      if (.not. ma_alloc_get(MT_DBL, (nr*nr)**2, 'dr', l_dr, k_dr))
     &     call errquit('WFAT: Cannot allocate dr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, (nt*nt)**2, 'dt', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt', 1, MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, nr*nt, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Could not allocate coef.', 34, MA_ERR)
      call ga_get(g_cf, 1, nr*nt, target, target, dbl_mb(k_cf), nr*nt)
      call wfat_2a1e_dpo('z', m, 0, nr, nt, nr, nt, sphi, sphi, ke,
     &     dbl_mb(k_cf), dbl_mb(k_cf), .false., dbl_mb(k_dr), .false.,
     &     dbl_mb(k_dt), .false., dipole(3))

      
c     ==== Analysis ====
      if (ndid == 0) then
         write(sm, '(i5)') m
         write(luout, '(2/)')
         call util_print_centered(luout,
     &        'Eigenvalue spectrum of m = '//trim(adjustl(sm))//
     &        ' manifold', 60, .true.)
         call wfat_2a1e_print_eigvals(nr*nt, 20, erg)
      endif
      
      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dr))
     &     call errquit('WFAT: Cannot deallocate dr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt', 1, MA_ERR)
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Could not deallocate coef', 7, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_print_eigvals(nb, nprint, erg)

      implicit none
#include "stdio.fh"
      integer nb
      integer nprint
      double precision erg(nb)

      character*10 ss
      character*19 sf
      integer i

      write(luout, '("")')
      write(ss, '(i10)') nprint
      write(luout, '(4x, a, " lowest energies : ")') trim(adjustl(ss))
      do i = 1, min(nb, nprint)
         if (mod(i,5) == 1) write(luout, '(6x)', advance='no')
         write(sf, '(sp, es19.8)') erg(i)
         write(luout, '(3x, i4, ") ", a)', advance='no')
     &        i, trim(adjustl(sf))
         if (mod(i,5) == 0) write(luout, '("")')
      enddo

      end
c===============================================================
