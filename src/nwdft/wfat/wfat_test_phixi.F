c=======================================================
c=======================================================
      subroutine wfat_test_phi1xi(n, m, erg, zq, chmmt)

      implicit none
#include "mafdecls.fh"
#include "stdio.fh"      
#include "errquit.fh"
#include "wfat.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      double precision chmmt(3)

      character*255 fname
      integer i, nxi, l_xi, k_xi, l_p1, k_p1
      double precision protmat(3,3), mu(3), kappa, beta0, beta1
      double precision dxi, xi_first, xi_last, bt, gm, xi, phi1, phi0
      double precision f2h, f1h, f0h, f1h_, f2h_, dif2, dif1, lhs, rhs
      double precision epsi
      parameter (bt=0.0d0, gm=0.0d0)
      

c     ==== Get the passive rotation matrix ===
      call wfat_passrot(0.0d0, 0.0d0, protmat)
      mu = matmul(transpose(protmat), -chmmt)
      
      
      kappa = sqrt(2 * abs(erg))
      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
      beta1 = -mu(3)/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      

c     ==== Calculate xi ====
      nxi = 101
      if (.not. MA_alloc_get(MT_DBL, nxi, 'xi', l_xi, k_xi))
     &     call errquit('WFAT: xi allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nxi, 'phi1', l_p1, k_p1))
     &     call errquit('WFAT: phi1 allocation failed', 86, MA_ERR)
      xi_first = 0.0d0
      xi_last = 10.0d0
      dxi = (xi_last-xi_first) / (nxi-1)
      do i = 1, nxi
         dbl_mb(k_xi+i-1) = xi_first + (i-1)*dxi
      enddo


      write(luout, '("TESTING FOR phi^{(1)}")')
      call util_file_name('test.wp1', .false., .false., fname)
      open(WF_PHI1TEST_UNIT, file=fname, action='write')
      write(WF_PHI1TEST_UNIT, '("nxi = ", i4)') nxi
      write(WF_PHI1TEST_UNIT, '("(n,m) = ", 2i4)') n, m
      write(WF_PHI1TEST_UNIT, '("kappa = ", es15.6)') kappa
      write(WF_PHI1TEST_UNIT, '("beta0 = ", es15.6)') beta0
      write(WF_PHI1TEST_UNIT, '("beta1 = ", es15.6)') beta1
      write(WF_PHI1TEST_UNIT, '("beta, gamma (deg) = ", 2es15.6)')
     &      bt, gm
      

c     ==== Calculate phi^{(1)} ====
      do i = 1, nxi
         xi = dbl_mb(k_xi + i - 1)
         call wfat_calc_phi1xi(xi, n, m, erg, bt, gm, chmmt,
     &        dbl_mb(k_p1+i-1))
      enddo


c     ==== Check differential equation ====
      write(WF_PHI1TEST_UNIT, '("#", a5, a15)', advance='no')
     &     'No.', 'kappa*xi'
      write(WF_PHI1TEST_UNIT, '(5a21)')
     &     'LHS', 'RHS', 'Epsilon', 'phi_\nu^(1)', 'phi_\nu^(0)'
      do i = 1, nxi-2
         f2h = dbl_mb(k_p1+(i+2)-1)
         f1h = dbl_mb(k_p1+(i+1)-1)
         f0h = dbl_mb(k_p1+(i+0)-1)
         f1h_ = dbl_mb(k_p1+(i-1)-1)
         f2h_ = dbl_mb(k_p1+(i-2)-1)
         dif2 = (-f2h + 16*f1h - 30*f0h + 16*f1h_ - f2h_)
     &          / (12*dxi**2)
         dif1 = (-f2h + 8*f1h - 8*f1h_ + f2h_) / (12*dxi)

         
         xi = dbl_mb(k_xi + i - 1)
         phi1 = f0h
         call wfat_calc_phi0xi(xi, n, m, erg, phi0)
         lhs = xi*dif2 + dif1 + (-m**2/(4*xi) + zq + erg*xi/2 - beta0) *
     &         phi1
         rhs = (mu(3)*xi/2 + xi**2/4 + beta1) * phi0
         epsi = lhs - rhs


         write(WF_PHI1TEST_UNIT, '(1x, i5, es15.6)', advance='no')
     &        i, kappa*xi
         write(WF_PHI1TEST_UNIT, '(5es21.12)')
     &        lhs, rhs, epsi, phi1, phi0
      enddo
      close(WF_PHI1TEST_UNIT)


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_xi))
     &     call errquit('WFAT: Cannot deallocate xi', 49, MA_ERR)
      if (.not. ma_free_heap(l_p1))
     &     call errquit('WFAT: Cannot deallocate phi1', 49, MA_ERR)
      
      end
c=======================================================
