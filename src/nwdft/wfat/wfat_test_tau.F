c=======================================================
c=======================================================
      subroutine wfat_test_tau(n, m, erg, zq, lmax, chmmt)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      integer lmax
      double precision chmmt(3)

      character*255 fname
      integer npt, i, i_r, i_th, i_ph, i_t0, i_t0l, i_t1, i_t1l
      integer l_r, l_th, l_ph, l_t0, l_t0l, l_t1, l_t1l
      integer k_r, k_th, k_ph, k_t0, k_t0l, k_t1, k_t1l
      double precision kappa, eta, xi, phi, bt, gm
      parameter(npt=6, bt=0.0d0, gm=10.0d0)


      write(luout, '("TESTING FOR TAU")')
      kappa = sqrt(2 * abs(erg))
      
      
c     ==== Allocate r, theta, and phi ====
      if (.not. MA_alloc_get(MT_DBL, npt, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'theta', l_th, k_th))
     &     call errquit('WFAT: theta allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'phi', l_ph, k_ph))
     &     call errquit('WFAT: phi allocation failed', 84, MA_ERR)

      dbl_mb(k_r) =    8.0d0
      dbl_mb(k_th) =   1.5d0
      dbl_mb(k_ph) =   3.2d0

      dbl_mb(k_r+1) =    8.0d0
      dbl_mb(k_th+1) =   0.2d0
      dbl_mb(k_ph+1) =   3.2d0
      
      dbl_mb(k_r+2) =  2.5d0
      dbl_mb(k_th+2) = 0.4d0
      dbl_mb(k_ph+2) = 1.2d0
      
      dbl_mb(k_r+3) =  1.8d1
      dbl_mb(k_th+3) = 2.1d0
      dbl_mb(k_ph+3) = 0.7d0
      
      dbl_mb(k_r+4) =  17.44d0
      dbl_mb(k_th+4) = 1.34d0
      dbl_mb(k_ph+4) = 3.93d0

      dbl_mb(k_r+5) =  1.0d1
      dbl_mb(k_th+5) = 2.9d0
      dbl_mb(k_ph+5) = 4.1d0

      
c     ==== Allocate tau's ====
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau0', l_t0, k_t0))
     &     call errquit('WFAT: tau0 allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau0_l', l_t0l, k_t0l))
     &     call errquit('WFAT: tau0_l allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau1', l_t1, k_t1))
     &     call errquit('WFAT: tau1 allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau1_l', l_t1l, k_t1l))
     &     call errquit('WFAT: tau1_l allocation failed', 84, MA_ERR)


c     ==== Compare tau_\nu^{(0)} ====
      call wfat_calc_tau0(n, m, npt, dbl_mb(k_r), dbl_mb(k_th),
     &     dbl_mb(k_ph), erg, zq, dcpl_mb(k_t0))
      call wfat_calc_tau0_lexp(n, m, npt, dbl_mb(k_r), dbl_mb(k_th),
     &     dbl_mb(k_ph), erg, zq, lmax, dcpl_mb(k_t0l))


c     ==== Compare tau_\nu^{(1)} ====
      call wfat_calc_tau1(n, m, npt, dbl_mb(k_r), dbl_mb(k_th),
     &     dbl_mb(k_ph), erg, zq, bt, gm, chmmt, dcpl_mb(k_t1))
      call wfat_calc_tau1_lexp(n, m, npt, dbl_mb(k_r), dbl_mb(k_th),
     &     dbl_mb(k_ph), erg, zq, bt, gm, chmmt, lmax, dcpl_mb(k_t1l))


c     ==== Open the test result file and print the header ====
      call util_file_name('test.wtau', .false., .false., fname)
      open(WF_TAUTEST_UNIT, file=fname, action='write')
      
      write(WF_TAUTEST_UNIT, '("n = ", i5)') n
      write(WF_TAUTEST_UNIT, '("m = ", i5)') m
      write(WF_TAUTEST_UNIT, '("energy = ", es15.6)') erg
      write(WF_TAUTEST_UNIT, '("ion charge = ", es15.6)') zq
      write(WF_TAUTEST_UNIT, '("lmax = ", i5)') lmax
      write(WF_TAUTEST_UNIT, '("orbital moment (x y z) = ", 3es15.6)')
     &     chmmt
      write(WF_TAUTEST_UNIT, '("beta, gamma (deg) = ", 2es15.6)') bt, gm

      
      call wfat_test_tau_header()      
      do i = 1, npt
c        == Print r, theta, phi ==
         i_r = k_r + i - 1
         i_th = k_th + i - 1
         i_ph = k_ph + i - 1
         write(WF_TAUTEST_UNIT, '(3f15.8)', advance='no')
     &        dbl_mb(i_r), dbl_mb(i_th), dbl_mb(i_ph)

c        == Print eta, xi, phi ==
         eta = dbl_mb(i_r) * (1 - cos( dbl_mb(i_th) ))
         xi = dbl_mb(i_r) * (1 + cos( dbl_mb(i_th) ))
         phi = dbl_mb(i_ph)
         write(WF_TAUTEST_UNIT, '(3f15.8)', advance='no') eta, xi, phi

c        == Print kappa*r, kappa*eta ==
         write(WF_TAUTEST_UNIT, '(2f15.8)', advance='no')
     &        kappa*dbl_mb(i_r), kappa*eta

c        == Print tau0 and tau1 ==
         i_t0 = k_t0 + i - 1
         i_t0l = k_t0l + i - 1
         i_t1 = k_t1 + i - 1
         i_t1l = k_t1l + i - 1
         write(WF_TAUTEST_UNIT, '(2(3x, 2es21.12))', advance='no')
     &        dcpl_mb(i_t0), dcpl_mb(i_t0l)
         write(WF_TAUTEST_UNIT, '(6x, " ")', advance='no')
         write(WF_TAUTEST_UNIT, '(2(3x, 2es21.12))', advance='no')
     &        dcpl_mb(i_t1), dcpl_mb(i_t1l)

         write(WF_TAUTEST_UNIT, '(" ")')
      enddo
      close(WF_TAUTEST_UNIT)

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_th))
     &     call errquit('WFAT: theta deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_ph))
     &     call errquit('WFAT: phi deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t0))
     &     call errquit('WFAT: tau0 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t0l))
     &     call errquit('WFAT: tau0_l deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t1))
     &     call errquit('WFAT: tau1 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t1l))
     &     call errquit('WFAT: tau1_l deallocation failed', 46, MA_ERR)
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_test_tau_header()

      implicit none
#include "wfat.fh"
#include "errquit.fh"

      logical isopen


c     ==== Check if the unit is opened ====
      inquire(unit=WF_TAUTEST_UNIT, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: WF_TAUTEST_UNIT must first be opened '//
     &     'before calling wfat_test_tau_header', 1, UERR)
      

c     ==== Column number ====
      write(WF_TAUTEST_UNIT, '("#", i14, 2i15)', advance='no') 1, 2, 3     ! r, theta, phi
      write(WF_TAUTEST_UNIT, '(3i15)', advance='no') 4, 5, 6               ! eta, xi, phi
      write(WF_TAUTEST_UNIT, '(2i15)', advance='no') 7, 8                  ! kappa*r, kappa*eta
      write(WF_TAUTEST_UNIT, '(2(3x, 2i21))', advance='no')
     &     9, 10, 11, 12                                                   ! tau0, tau0_pw
      write(WF_TAUTEST_UNIT, '(6x, " ")', advance='no')
      write(WF_TAUTEST_UNIT, '(2(3x, 2i21))', advance='no')
     &     13, 14, 15, 16                                                  ! tau1, tau1_pw
      write(WF_TAUTEST_UNIT, '(" ")')
      

c     ==== Column title ====
      write(WF_TAUTEST_UNIT, '("#", a14, 2a15)', advance='no')
     &     'r', 'theta', 'phi'
      write(WF_TAUTEST_UNIT, '(3a15)', advance='no')
     &     'eta', 'xi', 'phi'
      write(WF_TAUTEST_UNIT, '(2a15)', advance='no')
     &     'kappa*r', 'kappa*eta'
      write(WF_TAUTEST_UNIT, '(2(3x, 2a21))', advance='no')
     &     'real tau0', 'imag tau0', 'real tau0_pw', 'imag tau0_pw'
      write(WF_TAUTEST_UNIT, '(6x, " ")', advance='no')
      write(WF_TAUTEST_UNIT, '(2(3x, 2a21))', advance='no')
     &     'real tau1', 'imag tau1', 'real tau1_pw', 'imag tau1_pw'
      write(WF_TAUTEST_UNIT, '(" ")')

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_test2_tau(n, m, erg, zq, lmax, chmmt)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      integer lmax
      double precision chmmt(3)

      character*255 fname
      integer npt, nz, ny, l, i, j, i_r, i_t, i_p, i_t0, i_t0l, i_t1
      integer i_t1l
      integer l_r, l_t, l_p, l_t0, l_t0l, l_t1, l_t1l
      integer k_r, k_t, k_p, k_t0, k_t0l, k_t1, k_t1l
      double precision kappa, y0, yn, dy, z0, zn, dz, x, y, z, r, w
      double precision eta, bt, gm
      parameter(ny=1, nz=181, bt=30.0d0, gm=10.0d0)


      write(luout, '("TESTING FOR TAU v2")')
      kappa = sqrt(2 * abs(erg))

      x = 0.0d0
      
      y0 =   0.0d0
      yn =  24.0d0
      if (ny > 1) then
         dy = (yn-y0)/(ny-1)
      else
         dy = 0.0d0
      endif

      z0 = -50.0d0
      zn =  10.0d0
      if (nz > 1) then
         dz = (zn-z0)/(nz-1)
      else
         dz = 0.0d0
      endif
      
      npt = ny * nz

      
c     ==== Allocate r, theta, and phi ====
      if (.not. MA_alloc_get(MT_DBL, npt, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 24, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'theta', l_t, k_t))
     &     call errquit('WFAT: theta allocation failed', 24, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'phi', l_p, k_p))
     &     call errquit('WFAT: phi allocation failed', 24, MA_ERR)

      
c     ==== Allocate tau's ====
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau0', l_t0, k_t0))
     &     call errquit('WFAT: tau0 allocation failed', 24, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau0_l', l_t0l, k_t0l))
     &     call errquit('WFAT: tau0_l allocation failed', 24, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau1', l_t1, k_t1))
     &     call errquit('WFAT: tau1 allocation failed', 24, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau1_l', l_t1l, k_t1l))
     &     call errquit('WFAT: tau1_l allocation failed', 24, MA_ERR)
      

c     ==== Open the test result file and print the header ====
      call util_file_name('test.wtau2', .false., .false., fname)
      open(WF_TAUTEST2_UNIT, file=fname, action='write')      

      write(WF_TAUTEST2_UNIT, '("# ny, nz = ", 2i5)') ny, nz
      write(WF_TAUTEST2_UNIT, '("# n = ", i5)') n
      write(WF_TAUTEST2_UNIT, '("# m = ", i5)') m
      write(WF_TAUTEST2_UNIT, '("# energy = ", es15.6)') erg
      write(WF_TAUTEST2_UNIT, '("# ion charge = ", es15.6)') zq
      write(WF_TAUTEST2_UNIT, '("# lmax = ", i5)') lmax
      write(WF_TAUTEST2_UNIT, '("# orb. moment (x y z) = ", 3es15.6)')
     &     chmmt
      write(WF_TAUTEST2_UNIT, '("# beta, gamma (deg) = ", 2es15.6)')
     &     bt, gm
      

      call wfat_test2_tau_header()
      l = 0
      do i = 1, ny
         y = y0 + (i-1)*dy
         do j = 1, nz
            z = z0 + (j-1)*dz
            l = l + 1

            i_r = k_r + l - 1
            dbl_mb(i_r) = sqrt(x**2 + y**2 + z**2)
            i_t = k_t + l - 1
            dbl_mb(i_t) = acos(z/dbl_mb(i_r))
            i_p = k_p + l - 1
            dbl_mb(i_p) = atan(y, x)
         enddo
      enddo


c     ==== Compare tau_\nu^{(0)} ====
      call wfat_calc_tau0(n, m, npt, dbl_mb(k_r), dbl_mb(k_t),
     &     dbl_mb(k_p), erg, zq, dcpl_mb(k_t0))
      call wfat_calc_tau0_lexp(n, m, npt, dbl_mb(k_r), dbl_mb(k_t),
     &     dbl_mb(k_p), erg, zq, lmax, dcpl_mb(k_t0l))


c     ==== Compare tau_\nu^{(1)} ====
      call wfat_calc_tau1(n, m, npt, dbl_mb(k_r), dbl_mb(k_t),
     &     dbl_mb(k_p), erg, zq, bt, gm, chmmt, dcpl_mb(k_t1))
      call wfat_calc_tau1_lexp(n, m, npt, dbl_mb(k_r), dbl_mb(k_t),
     &     dbl_mb(k_p), erg, zq, bt, gm, chmmt, lmax, dcpl_mb(k_t1l))


      l = 0
      do i = 1, ny
         y = y0 + (i-1)*dy
         do j = 1, nz
            z = z0 + (j-1)*dz
            write(WF_TAUTEST2_UNIT, '(i10, 3f15.8)', advance='no')
     &           l, x, y, z
            r = sqrt(x**2 + y**2 + z**2)
            w = exp(-kappa*r)
            eta = r - z

c           == Print kappa*r, kappa*eta ==
            write(WF_TAUTEST2_UNIT, '(2f15.8)', advance='no')
     &           kappa*r, kappa*eta
            
            l = l + 1
c           == Print tau0 and tau1 ==
            i_t0  = k_t0  + l - 1
            i_t0l = k_t0l + l - 1
            i_t1  = k_t1  + l - 1
            i_t1l = k_t1l + l - 1
            write(WF_TAUTEST2_UNIT, '(4(3x, 2es21.12))', advance='no')
     &           dcpl_mb(i_t0), dcpl_mb(i_t0l),
     &           w*dcpl_mb(i_t0), w*dcpl_mb(i_t0l)
            write(WF_TAUTEST2_UNIT, '(6x, " ")', advance='no')
            write(WF_TAUTEST2_UNIT, '(4(3x, 2es21.12))', advance='no')
     &           dcpl_mb(i_t1), dcpl_mb(i_t1l),
     &           w*dcpl_mb(i_t1), w*dcpl_mb(i_t1l)
            
            write(WF_TAUTEST2_UNIT, '(" ")')
         enddo
         write(WF_TAUTEST2_UNIT, '(" ")')
      enddo

      close(WF_TAUTEST2_UNIT)

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 73, MA_ERR)
      if (.not. MA_free_heap(l_t))
     &     call errquit('WFAT: theta deallocation failed', 73, MA_ERR)
      if (.not. MA_free_heap(l_p))
     &     call errquit('WFAT: phi deallocation failed', 73, MA_ERR)
      if (.not. MA_free_heap(l_t0))
     &     call errquit('WFAT: tau0 deallocation failed', 73, MA_ERR)
      if (.not. MA_free_heap(l_t0l))
     &     call errquit('WFAT: tau0_l deallocation failed', 73, MA_ERR)
      if (.not. MA_free_heap(l_t1))
     &     call errquit('WFAT: tau1 deallocation failed', 73, MA_ERR)
      if (.not. MA_free_heap(l_t1l))
     &     call errquit('WFAT: tau1_l deallocation failed', 73, MA_ERR)

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_test2_tau_header()

      implicit none
#include "wfat.fh"
#include "errquit.fh"

      logical isopen


c     ==== Check if the unit is opened ====
      inquire(unit=WF_TAUTEST2_UNIT, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: WF_TAUTEST2_UNIT must first be opened '//
     &     'before calling wfat_test_tau_header', 1, UERR)
      

c     ==== Column number ====
      write(WF_TAUTEST2_UNIT, '("#", i9, 3i15)', advance='no')
     &     1, 2, 3, 4                                                       ! r, theta, phi
      write(WF_TAUTEST2_UNIT, '(2i15)', advance='no') 5, 6                  ! kappa*r, kappa*eta
      write(WF_TAUTEST2_UNIT, '(4(3x, 2i21))', advance='no')
     &     7, 8, 9, 10, 11, 12, 13, 14                                      ! tau0, tau0_pw
      write(WF_TAUTEST2_UNIT, '(6x, " ")', advance='no')
      write(WF_TAUTEST2_UNIT, '(4(3x, 2i21))', advance='no')
     &     15, 16, 17, 18, 19, 20, 21, 22                                   ! tau1, tau1_pw
      write(WF_TAUTEST2_UNIT, '(" ")')
      

c     ==== Column title ====
      write(WF_TAUTEST2_UNIT, '("#", a9, 3a15)', advance='no')
     &     'No.', 'x', 'y', 'z'
      write(WF_TAUTEST2_UNIT, '(2a15)', advance='no')
     &     'kappa*r', 'kappa*eta'
      write(WF_TAUTEST2_UNIT, '(4(3x, 2a21))', advance='no')
     &     'Re tau0', 'Im tau0', 'Re tau0_pw', 'Im tau0_pw',
     &     'Re w*tau0', 'Im w*tau0', 'Re w*tau0_pw', 'Im w*tau0_pw'
      write(WF_TAUTEST2_UNIT, '(6x, " ")', advance='no')
      write(WF_TAUTEST2_UNIT, '(4(3x, 2a21))', advance='no')
     &     'Re tau1', 'Im tau1', 'Re tau1_pw', 'Im tau1_pw',
     &     'Re w*tau1', 'Im w*tau1', 'Re w*tau1_pw', 'Im w*tau1_pw'
      write(WF_TAUTEST2_UNIT, '(" ")')

      end
c=======================================================
