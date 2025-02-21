c=======================================================
c=======================================================
      subroutine wfat_test_q(n, m, erg, zq, lmax, chmmt)

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
      integer nr, nlgl, i, nl, l
      integer l_r, l_qp1, l_qp2, k_r, k_qp1, k_qp2
      integer i_q1_2h, i_q1_1h, i_q1_0h, i_q1_1h_, i_q1_2h_
      integer i_q2_2h, i_q2_1h, i_q2_0h, i_q2_1h_, i_q2_2h_
      double precision r_first, r_last, kappa, lmb, a, b, dr, r, epsi
      double precision f2h, f1h, f0h, f1h_, f2h_, dif1, dif2, bt, gm
      double precision k_m, k_p, c_m, c_p, c, r0_m, r0_p, r0, lhs, rhs
      double precision protmat(3,3), mu(3), chgm(abs(m):lmax+1)
      integer nrr, k
      parameter (nrr=4)
      integer ii(nrr)
      double precision rr(nrr)
      double precision q_p1_(abs(m):lmax,nrr), q_p2_(abs(m):lmax,nrr) 
      

      nr = 10001
      nlgl = 10

      r_first = 0.0d0
      r_last = 100.0d0

      bt = 0.0d0
      gm = 0.0d0

      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0

      write(luout, '("TESTING FOR Q_l")')

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      mu = matmul(transpose(protmat), -chmmt)      
      

c     ==== Generate the radial points ====
      if (.not. MA_alloc_get(MT_DBL, nr, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 52, MA_ERR)
      dr = (r_last - r_first) / (nr-1)
      do i = 1, nr
         dbl_mb(k_r+i-1) = r_first + (i-1)*dr
      enddo

      nl = lmax - abs(m) + 1
      if (.not. MA_alloc_get(MT_DBL, nl*nr, 'q_p1', l_qp1, k_qp1))
     &     call errquit('WFAT: q_p1 allocation failed', 52, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl*nr, 'q_p2', l_qp2, k_qp2))
     &     call errquit('WFAT: q_p2 allocation failed', 52, MA_ERR)
            

c     ==== Calculate Q_l ====
      call wfat_q_particular(nr, dbl_mb(k_r), n, m, lmax, erg, zq,
     &     nlgl, dbl_mb(k_qp1), dbl_mb(k_qp2))
      

c     ==== Open the test result file and print the header ====
      call util_file_name('test.wq', .false., .false., fname)
      open(WF_QTEST_UNIT, file=fname, action='write')
      write(WF_QTEST_UNIT, '("# nr = ", i4)') nr
      write(WF_QTEST_UNIT, '("# nlgl = ", i4)') nlgl
      write(WF_QTEST_UNIT, '("# n = ", i5)') n
      write(WF_QTEST_UNIT, '("# m = ", i5)') m
      write(WF_QTEST_UNIT, '("# energy = ", es15.6)') erg
      write(WF_QTEST_UNIT, '("# kappa = ", es15.6)') kappa
      write(WF_QTEST_UNIT, '("# ion charge = ", es15.6)') zq
      write(WF_QTEST_UNIT, '("# lmax = ", i5)') lmax
      write(WF_QTEST_UNIT, '("# orbital moment (x y z) = ", 3es15.6)')
     &     chmmt
      write(WF_QTEST_UNIT, '("# beta, gamma (deg) = ", 2es15.6)') bt, gm

      call wfat_test_q_header(m, lmax)      
      

c     ==== Loop over radial points ====
      do i = 1, nr-2
         r = dbl_mb(k_r+i-1)
         write(WF_QTEST_UNIT, '(1x, i5, es15.6)', advance='no') i, lmb*r
         call wfat_chgm_1pt(r, abs(m), lmax+1, a, b, lmb, chgm)

         do l = abs(m), lmax
c           ==== Calculate the derivatives ====
            i_q1_2h = k_qp1 + ((i+2)-1)*nl + l - abs(m)
            i_q1_1h = k_qp1 + ((i+1)-1)*nl + l - abs(m)
            i_q1_0h = k_qp1 + ((i+0)-1)*nl + l - abs(m)
            if (i > 1) i_q1_1h_ = k_qp1 + ((i-1)-1)*nl + l - abs(m)
            if (i > 2) i_q1_2h_ = k_qp1 + ((i-2)-1)*nl + l - abs(m)

            i_q2_2h = k_qp2 + ((i+2)-1)*nl + l - abs(m)
            i_q2_1h = k_qp2 + ((i+1)-1)*nl + l - abs(m)
            i_q2_0h = k_qp2 + ((i+0)-1)*nl + l - abs(m)
            if (i > 1) i_q2_1h_ = k_qp2 + ((i-1)-1)*nl + l - abs(m)
            if (i > 2) i_q2_2h_ = k_qp2 + ((i-2)-1)*nl + l - abs(m)

            f2h = dbl_mb(i_q1_2h) + mu(3)*dbl_mb(i_q2_2h)
            f1h = dbl_mb(i_q1_1h) + mu(3)*dbl_mb(i_q2_1h)
            f0h = dbl_mb(i_q1_0h) + mu(3)*dbl_mb(i_q2_0h)
            if (i > 1) f1h_ = dbl_mb(i_q1_1h_) + mu(3)*dbl_mb(i_q2_1h_)
            if (i > 2) f2h_ = dbl_mb(i_q1_2h_) + mu(3)*dbl_mb(i_q2_2h_)

            if (i == 1) then
               dif2 = 0.0d0
               dif1 = 0.0d0
            else if (i == 2) then
               dif2 = (f1h - 2*f0h + f1h_) / (dr**2)
               dif1 = (f1h - f1h_) / (2*dr)
            else
               dif2 = (-f2h + 16*f1h - 30*f0h + 16*f1h_ - f2h_)
     &                / (12*dr**2)
               dif1 = (-f2h + 8*f1h - 8*f1h_ + f2h_) / (12*dr)
            endif


c           ==== Calculate the terms on the RHS ====
            if (l == abs(m)) then
               r0_m = 0.0d0
               k_m = 0.0d0
            else
               call wfat_calc_radnorm(l-1, n, m, int(zq), kappa, c_m)
               k_m = sqrt(dble(l**2-m**2) / (4*l**2-1))
               r0_m = c_m * (kappa*r)**(l-1)*exp(-kappa*r) * chgm(l-1)
            endif

            call wfat_calc_radnorm(l+1, n, m, int(zq), kappa, c_p)
            k_p = sqrt(dble((l+1)**2-m**2) / (4*(l+1)**2-1))
            r0_p = c_p * (kappa*r)**(l+1)*exp(-kappa*r) * chgm(l+1)

            call wfat_calc_radnorm(l, n, m, int(zq), kappa, c)
            r0 = c * (kappa*r)**l*exp(-kappa*r) * chgm(l)


c           ==== Calculate the LHS and RHS ====
            lhs = r**2*dif2 + 2*r*dif1 +
     &            (-l*(l+1) + 2*zq*r + 2*erg*r**2) * f0h
            rhs = 2*k_m*r**3*r0_m + 2*k_p*r**3*r0_p + 2*mu(3)*r**2*r0
            epsi = lhs - rhs

            write(WF_QTEST_UNIT, '(6x, 5es21.12)', advance='no')
     &           lhs, rhs, epsi, f0h, r0
         enddo

         write(WF_QTEST_UNIT, '(" ")')
      enddo



      ii(1) = int(1.0d0/4 * nr)
      rr(1) = dbl_mb(k_r + ii(1) - 1)
      ii(2) = int(2.0d0/4 * nr)
      rr(2) = dbl_mb(k_r + ii(2) - 1)
      ii(3) = int(3.0d0/4 * nr)
      rr(3) = dbl_mb(k_r + ii(3) - 1)
      ii(4) = int(3.8d0/4 * nr)
      rr(4) = dbl_mb(k_r + ii(4) - 1)
      call wfat_q_particular(nrr, rr, n, m, lmax, erg, zq,
     &     nlgl, q_p1_, q_p2_)
      do i = 1, nrr
         write(WF_QTEST_UNIT, '("#", 2x, i7, es15.6)', advance='no')
     &        ii(i), lmb*rr(i)
         do l = abs(m), lmax
            i_q1_0h = k_qp1 + ((ii(i)+0)-1)*nl + l - abs(m)
            i_q2_0h = k_qp2 + ((ii(i)+0)-1)*nl + l - abs(m)
            write(WF_QTEST_UNIT, '(8x, 2es15.6)', advance='no')
     &           q_p1_(l,i) + mu(3)*q_p2_(l,i), 
     &           dbl_mb(i_q1_0h) + mu(3)*dbl_mb(i_q2_0h)
         enddo
         write(WF_QTEST_UNIT, '(" ")')
      enddo
      
      close(WF_QTEST_UNIT)

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_r))
     &     call errquit('WFAT: Cannot deallocate r', 49, MA_ERR)
      if (.not. ma_free_heap(l_qp1))
     &     call errquit('WFAT: Cannot deallocate q_p1', 49, MA_ERR)
      if (.not. ma_free_heap(l_qp2))
     &     call errquit('WFAT: Cannot deallocate q_p2', 49, MA_ERR)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_test_q_header(m, lmax)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer m
      integer lmax

      logical isopen
      integer i, j, l


c     ==== Check if the unit is opened ====
      inquire(unit=WF_QTEST_UNIT, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: WF_QTEST_UNIT must first be opened '//
     &     'before calling wfat_test_q_header', 1, UERR)
      

c     ==== Column number ====
      write(WF_QTEST_UNIT, '("#", i5, i15)', advance='no') 1, 2
      do i = 1, lmax-abs(m)+1
         j = 5*i - 2
         write(WF_QTEST_UNIT, '(6x, 5i21)', advance='no')
     &        j, j+1, j+2, j+3, j+4
      enddo
      write(WF_QTEST_UNIT, '(" ")')


c     ==== Angular momenta ====
      write(WF_QTEST_UNIT, '("#", a20)', advance='no') " "
      do l = abs(m), lmax
         write(WF_QTEST_UNIT, '(6x, a100, i5)', advance='no')
     &        "l = ", l
      enddo
      write(WF_QTEST_UNIT, '(" ")')
      

c     ==== Column title ====
      write(WF_QTEST_UNIT, '("#", a5, a15)', advance='no') 'No.', 'y'
      do l = abs(m), lmax
         write(WF_QTEST_UNIT, '(6x, 5a21)', advance='no')
     &        'LHS', 'RHS', 'epsilon', 'Q_l', 'R_l'
      enddo
      write(WF_QTEST_UNIT, '(" ")')


      end
c=======================================================
      
      
