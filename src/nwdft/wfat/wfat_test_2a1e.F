c==========================================================
c==========================================================
      subroutine wfat_test_2a1e(m, nr1, nt1, nr2, nt2, nr3, nt3)

      implicit none
#include "geom.fh"
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr1
      integer nt1
      integer nr2
      integer nt2
      integer nr3
      integer nt3

      character*16 tag(2)
      logical numz
      integer l_cf1, l_eg1, l_cf2, l_eg2, l_cf3, l_eg3, l_dr, l_dt
      integer k_cf1, k_eg1, k_cf2, k_eg2, k_cf3, k_eg3, k_dr, k_dt
      integer geom, nr1_, nt1_, i, maxiter, orbid, i_cf2, sphi
      double precision ch(2), cd(3,2), potparam(2)
      double precision q1, q2, z1, z2, ke, polb(3), eps, dipole(3)


      ke = sqrt(2 * 0.5d0)
      ke = 1.0d0
      write(luout, '("Initial ke = ", f18.8, " /Bohr")') ke
      call wfat_test_dvr_l0(nr2, ke)
      call wfat_test_fbr_l0(nr2, ke)
      
      write(luout, '(a)') 'Testing 2a1e ...'

c     ==== H_2^+ ====
      q1 = 1.0d0
      q2 = 1.0d0
      z1 = 1.0d0
      z2 = -1.0d0
      orbid = 1
      eps = 0.09d0
      potparam = (/ eps, -999.0d0 /)
      sphi = 2

c     ==== HeH^{2+} ====
c      q1 = 1.0d0
c      q2 = 2.0d0
c      z1 = 3.112d0
c      z2 = -0.778d0
c      orbid = 2
c      eps = 0.1d0
c      sphi = 1


c     ==== Set up the geometry ====
      if (.not. geom_create(geom, 'wfat_test_geom'))
     &     call errquit('WFAT: Error in creating geometry object.', 50,
     &     GEOM_ERR)
      tag = (/ 'atom1', 'atom2' /)
      cd(:,1) = (/ 0.0d0, 0.0d0, z1 /)
      cd(:,2) = (/ 0.0d0, 0.0d0, z2 /)
      ch = (/ q1, q2 /)
      if (.not. geom_cart_set(geom, 2, tag, cd, ch))
     &     call errquit('WFAT: Error in initializing geometry object.',
     &     39, GEOM_ERR)

      
c     ==== Obtain eigenstates with m'=m ====
      maxiter = 10
      if (.not. ma_alloc_get(MT_DBL, (nr2*nt2)**2, 'hm2', l_cf2, k_cf2))
     &     call errquit('WFAT: Cannot allocate hm2', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr2*nt2, 'erg2', l_eg2, k_eg2))
     &     call errquit('WFAT: Cannot allocate erg2', 1, MA_ERR)
      call wfat_2a1e_solve(m, nr2, nt2, geom, ke, .true., maxiter,
     &     orbid, 'cou', 1, potparam, dbl_mb(k_cf2), dbl_mb(k_eg2))


c     ==== Print eigenvalues ====
      write(luout, '("Output ke = ", f18.8, " /Bohr")') ke
      write(luout, '("First 10 eigenvalues for m = ", i3, " and with '
     &     //'ke adaptation:")') m
      do i = 1, 10
         write(luout, '(i4, f18.8)') i, dbl_mb(k_eg2+i-1)
      enddo
      write(luout, '("Eigenvalue range = ", f18.8, " H")')
     &     dbl_mb(k_eg2+nr2*nt2-1) - dbl_mb(k_eg2)

      dipole = 0.0d0
      i_cf2 = k_cf2 + (orbid-1)*nr2*nt2
      if (.not. ma_alloc_get(MT_DBL, (nr2*nr2)**2, 'dr', l_dr, k_dr))
     &     call errquit('WFAT: Cannot allocate dr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, (nt2*nt2)**2, 'dt', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt', 1, MA_ERR)
      call wfat_2a1e_dpo('z', m, 0, nr2, nt2, nr2, nt2, sphi, sphi, ke,
     &     dbl_mb(i_cf2), dbl_mb(i_cf2), .False., dbl_mb(k_dr), .False.,
     &     dbl_mb(k_dt), .True., dipole(3))
      write(luout, '("Dipole moment = ", 3es17.8)') dipole

      
c     ==== Obtain eigenstates with m'=m-1 ====
      if (m == 0) then
         nr1_ = 1
         nt1_ = 1
      else if (m >= 1) then
         nr1_ = nr1
         nt1_ = nt1
      endif
      if (.not. ma_alloc_get(MT_DBL, (nr1_*nt1_)**2, 'hm1', l_cf1,
     &     k_cf1))
     &     call errquit('WFAT: Cannot allocate hm1', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr1_*nt1_, 'erg1', l_eg1, k_eg1))
     &     call errquit('WFAT: Cannot allocate erg1', 1, MA_ERR)
      if (m >= 1) then
         write(luout, '("Input ke = ", es17.8)') ke
         call wfat_2a1e_solve(m-1, nr1_, nt1_, geom, ke, .False., 1, 1,
     &        'cou', 1, potparam, dbl_mb(k_cf1), dbl_mb(k_eg1))
         write(luout, '("Output ke = ", es17.8)') ke
         write(luout, '("First 10 eigenvalues for m = ", i3, " :")') m-1
         do i = 1, 10
            write(luout, '(i4, f18.8)') i, dbl_mb(k_eg1+i-1)
         enddo
      endif

      
c     ==== Obtain eigenstates with m'=m+1 ====
      if (.not. ma_alloc_get(MT_DBL, (nr3*nt3)**2, 'hm3', l_cf3, k_cf3))
     &     call errquit('WFAT: Cannot allocate hm3', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr3*nt3, 'erg3', l_eg3, k_eg3))
     &     call errquit('WFAT: Cannot allocate erg3', 1, MA_ERR)
      write(luout, '("Input ke = ", es17.8)') ke
      call wfat_2a1e_solve(m+1, nr3, nt3, geom, ke, .False., 1, 1,
     &     'cou', 1, potparam, dbl_mb(k_cf3), dbl_mb(k_eg3))
      write(luout, '("Output ke = ", es17.8)') ke
      write(luout, '("First 10 eigenvalues for m = ", i3, " :")') m+1
      do i = 1, 10
         write(luout, '(i4, f18.8)') i, dbl_mb(k_eg3+i-1)
      enddo


c     ==== Calculate the static polarizability ====
      numz = .false.
cOLD      call wfat_2a1e_polb('c', m, 1, nr, nt, ke, dbl_mb(k_eg),
cOLD     &     dbl_mb(k_cf), numz, polb)
cOLD      call wfat_2a1e_polb(m, sphi, orbid, nr1_, nt1_, nr2, nt2, nr3, 
cOLD     &     nt3, ke, dbl_mb(k_eg1), dbl_mb(k_eg2), dbl_mb(k_eg3), 
cOLD     &     dbl_mb(k_cf1), dbl_mb(k_cf2), dbl_mb(k_cf3), numz, polb)
cOLD      write(luout, '("MF static polarizability = ", 3es17.8)') polb

      

c     ==== Don't adapt ke ====
      q1 = 1.0d0
      q2 = 1.0d0
      z1 = 1.0d0
      z2 = -1.0d0
      cd(:,1) = (/ 0.0d0, 0.0d0, z1 /)
      cd(:,2) = (/ 0.0d0, 0.0d0, z2 /)
      ch = (/ q1, q2 /)
      if (.not. geom_cart_set(geom, 2, tag, cd, ch))
     &     call errquit('WFAT: Error in initializing geometry object.',
     &     39, GEOM_ERR)
      
      maxiter = 10
      call wfat_2a1e_solve(0, nr2, nt2, geom, ke, .false., 0, 1, 'cou',
     &     1, potparam, dbl_mb(k_cf2), dbl_mb(k_eg2))
      write(luout, '("The lowest eigenvalue for m = 0 and without ke '//
     &     'adaptation = ", f18.8, " H")') dbl_mb(k_eg2)

      
c     ==== Deallocations ====
      if (.not. geom_destroy(geom))
     &     call errquit('WFAT: Error in destroying geometry object.',
     &     50, GEOM_ERR)
      
      if (m >= 1) then
         if (.not. ma_free_heap(l_cf1))
     &        call errquit('WFAT: Cannot deallocate hm1', 1, MA_ERR)
         if (.not. ma_free_heap(l_eg1))
     &        call errquit('WFAT: Cannot deallocate erg1', 1, MA_ERR)
      endif
      if (.not. ma_free_heap(l_cf2))
     &     call errquit('WFAT: Cannot deallocate hm2', 1, MA_ERR)
      if (.not. ma_free_heap(l_eg2))
     &     call errquit('WFAT: Cannot deallocate erg2', 1, MA_ERR)
      if (.not. ma_free_heap(l_cf3))
     &     call errquit('WFAT: Cannot deallocate hm3', 1, MA_ERR)
      if (.not. ma_free_heap(l_eg3))
     &     call errquit('WFAT: Cannot deallocate erg3', 1, MA_ERR)
      if (.not. ma_free_heap(l_dr))
     &     call errquit('WFAT: Cannot deallocate dr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt', 1, MA_ERR)
      
      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_test_dvr_l0(nr, ke)
c     DESCRIPTION:
c     This subroutine tests the computation of the eigenvalues of the hydrogen
c     Hamiltonian matrix in the DVR bases.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer nr
      double precision ke

      integer l_lr, l_ur, l_hm, l_x, l_rw, l_eg, l_po0, l_po
      integer k_lr, k_ur, k_hm, k_x, k_rw, k_eg, k_po0, k_po
      integer i, j, i_hm, i_po0
      double precision x, pot, normi, normj, v, laguerre_norm_i, mm
      external laguerre_norm_i


      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'lapr', l_lr, k_lr))
     &     call errquit('WFAT: Cannot allocate lapr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'u_r', l_ur, k_ur))
     &     call errquit('WFAT: Cannot allocate u_r', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'g_r', l_hm, k_hm))
     &     call errquit('WFAT: Cannot allocate g_r', 1, MA_ERR)


      call wfat_laplr_fbrr(nr, 0, dbl_mb(k_lr))
      call wfat_dvr_fbr_r(0, nr, dbl_mb(k_ur))
c     T_r^{DVR} = transpose(U_r) . T_r^{FBR} . U_r
      call wfat_represent_sym(nr, nr, dbl_mb(k_ur), dbl_mb(k_lr),
     &     dbl_mb(k_hm))
      dbl_mb(k_hm:k_hm+nr**2-1) = -0.5d0 * dbl_mb(k_hm:k_hm+nr**2-1)
     &     * ke**2
c      write(luout, '("Kinetic energy matrix:")')
c      call wfat_dprintmat(nr, nr, 'f', dbl_mb(k_hm))


      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'xweight', l_rw, k_rw))
     &     call errquit('WFAT: Cannot allocate xweight', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)
      

      call wfat_lf_quadrature_rule(nr, 0+2.0d0, dbl_mb(k_x),
     &     dbl_mb(k_rw))
      do i = 1, nr
         x = dbl_mb(k_x+i-1)
         pot = - ke / x
         i_hm = k_hm + (i-1)*nr + (i-1)
         dbl_mb(i_hm) = dbl_mb(i_hm) + pot
      enddo


c     ================================================================
c     ==== The block of code below computes the potential matrix  ====
c     ==== elements analytically instead of numerically, that is, ====
c     ==== it does not use the special property of the DVR bases  ====
c     ====  in which the matrix of a local operator is diagonal.  ====
c     ====  Remove the comment below but comment out the entire   ====
c     ====  loop above if you want to use the analytic potential  ====
c     ==== matrix. Using the analytic potential matrix will give  ====
c     ====  a more accurate result which should be identical to   ====
c     ==== the results obtained by diagonalizing in the FBR basis ====
c     ==== (which is done by the wfat_test_fbr_l0 routine) at the ====
c     ====          expense of a longer simulation time.          ====
c     ================================================================
cc      if (.not. ma_alloc_get(MT_DBL, nr**2, 'pot0', l_po0, k_po0))
cc     &     call errquit('WFAT: Cannot allocate pot0', 1, MA_ERR)
cc      if (.not. ma_alloc_get(MT_DBL, nr**2, 'pot', l_po, k_po))
cc     &     call errquit('WFAT: Cannot allocate pot', 1, MA_ERR)
cc      mm = dble(0)
cc      do j = 1, nr
cc         normj = sqrt(laguerre_norm_i(j-1,0+2))
cc         do i = 1, nr
cc            normi = sqrt(laguerre_norm_i(i-1,0+2))
cc            call wfat_laguerre_intg1(mm+2, mm+2, mm+1, i-1, j-1, v)
cc            i_po0 = k_po0 + (j-1)*nr + (i-1)
cc            dbl_mb(i_po0) = - v / normi / normj / ke**2
cc         enddo
cc      enddo
cc      call wfat_represent_sym(nr, nr, dbl_mb(k_ur), dbl_mb(k_po0),
cc     &     dbl_mb(k_po))
cc
cc      dbl_mb(k_hm:k_hm+nr**2-1) = dbl_mb(k_hm:k_hm+nr**2-1) +
cc     &                            dbl_mb(k_po:k_po+nr**2-1)
      
      
      call wfat_dsyev_wrap(nr, dbl_mb(k_hm), dbl_mb(k_eg))
      write(luout, '("Hydrogen eigenvalues for m = 0 using DVR :")')
      do i = 1, nr
         if (dbl_mb(k_eg+i-1) < 0.0d0)
     &        write(luout, '(i4, f18.8)') i, dbl_mb(k_eg+i-1)
      enddo


      if (.not. ma_free_heap(l_lr))
     &     call errquit('WFAT: Cannot deallocate lapr', 1, MA_ERR)
      if (.not. ma_free_heap(l_ur))
     &     call errquit('WFAT: Cannot deallocate u_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_hm))
     &     call errquit('WFAT: Cannot deallocate hm', 1, MA_ERR)
      if (.not. ma_free_heap(l_x))
     &     call errquit('WFAT: Cannot deallocate xquad', 1, MA_ERR)
      if (.not. ma_free_heap(l_rw))
     &     call errquit('WFAT: Cannot deallocate xweight', 1, MA_ERR)
      if (.not. ma_free_heap(l_eg))
     &     call errquit('WFAT: Cannot deallocate erg', 1, MA_ERR)

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_test_fbr_l0(nr, ke)
c     DESCRIPTION:
c     This subroutine tests the computation of the eigenvalues of the hydrogen
c     Hamiltonian matrix in the FBR bases.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer nr
      double precision ke

      integer l_hm, k_hm, l_eg, k_eg, i, j, i_hm
      double precision mm, normi, normj, pot, v, laguerre_norm_i
      external laguerre_norm_i


      ke = sqrt(2 * 0.5d0)
      mm = dble(0)
      
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'hm', l_hm, k_hm))
     &     call errquit('WFAT: Cannot allocate hm', 1, MA_ERR)

      
      call wfat_laplr_fbrr(nr, 0, dbl_mb(k_hm))
      dbl_mb(k_hm:k_hm+nr**2-1) = -0.5d0 * dbl_mb(k_hm:k_hm+nr**2-1)
     &     / ke

      do j = 1, nr
         normj = sqrt(laguerre_norm_i(j-1,0+2))
         do i = 1, nr
            normi = sqrt(laguerre_norm_i(i-1,0+2))
            call wfat_laguerre_intg1(mm+2, mm+2, mm+1, i-1, j-1, v)
            pot = - v / normi / normj / ke**2
            i_hm = k_hm + (j-1)*nr + (i-1)
            dbl_mb(i_hm) = dbl_mb(i_hm) + pot
         enddo
      enddo


      if (.not. ma_alloc_get(MT_DBL, nr, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)
      call wfat_dsyev_wrap(nr, dbl_mb(k_hm), dbl_mb(k_eg))
      write(luout, '("Hydrogen eigenvalues with m = 0 using FBR :")')
      do i = 1, nr
         if (dbl_mb(k_eg+i-1) < 0.0d0)
     &        write(luout, '(i4, f18.8)') i, dbl_mb(k_eg+i-1)
      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_hm))
     &     call errquit('WFAT: Cannot deallocate hm', 32, MA_ERR)
      
      end
c==========================================================
      
      
