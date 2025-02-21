c=======================================================
c=======================================================
      subroutine wfat_test_r1(n, m, erg, zq, chmmt)

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
      integer l_et, k_et, l_r1, k_r1
      integer nlgl, neta, i, j
      double precision kappa, a, b, et_first, et_last, de, et, z
      double precision bt, gm, chgm, mu(3), f, o0_1, r0, r1
      double precision protmat(3,3), beta0, beta1, d, g
      double precision cr_1_ref, f2h, f1h, f0h, f1h_, f2h_, dif2, lhs
      double precision rhs, epsi
      double complex zz, chgm_
      double precision factorial
      external factorial
      parameter(bt=0.0d0, gm=0.0d0)


c     ==== Get the passive rotation matrix ===
      call wfat_passrot(0.0d0, 0.0d0, protmat)
      mu = matmul(transpose(protmat), -chmmt)

      nlgl = 10
      
      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
      beta1 = -mu(3)/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))

      neta = 301
      if (.not. MA_alloc_get(MT_DBL, neta, 'eta', l_et, k_et))
     &     call errquit('WFAT: eta allocation failed', 84, MA_ERR)
      et_first = 0.0d0
      et_last = 15.0d0
      de = (et_last-et_first) / (neta-1)
      do i = 1, neta
         dbl_mb(k_et+i-1) = et_first + (i-1)*de
      enddo

      
      write(luout, '("TESTING FOR R^{(1)}")')
      call util_file_name('test.wr1', .false., .false., fname)
      open(WF_R1TEST_UNIT, file=fname, action='write')
      write(WF_R1TEST_UNIT, '("npt = ", i4)') neta
      write(WF_R1TEST_UNIT, '("(n,m) = ", 2i4)') n, m
      write(WF_R1TEST_UNIT, '("a, b = ", 2es15.6)') a, b
      write(WF_R1TEST_UNIT, '("kappa = ", es15.6)') kappa
      write(WF_R1TEST_UNIT, '("beta0 = ", es15.6)') beta0
      write(WF_R1TEST_UNIT, '("beta1 = ", es15.6)') beta1
      write(WF_R1TEST_UNIT, '("beta, gamma (deg) = ", 2es15.6)') bt, gm


      if (.not. MA_alloc_get(MT_DBL, neta, 'r1', l_r1, k_r1))
     &     call errquit('WFAT: r1 allocation failed', 86, MA_ERR)      
      call wfat_r1(neta, dbl_mb(k_et), n, m, erg, zq, bt, gm, chmmt,
     &     nlgl, dbl_mb(k_r1))
            
cfordebug      cr_1_ref = -g/kappa/factorial(abs(m)+1) * beta1


c     ==== Check differential equation ====
      write(WF_R1TEST_UNIT, '("#", a5, a15)', advance='no')
     &     'No.', 'kappa*eta'
      write(WF_R1TEST_UNIT, '(5a21)')
     &     'LHS', 'RHS', 'Epsilon', 'R_\nu^(1)', 'R_\nu^(0)'
      do i = 1, neta-2
         if (i <= 2) then
            dif2 = 0.0d0
         else
            f2h = dbl_mb(k_r1+(i+2)-1)
            f1h = dbl_mb(k_r1+(i+1)-1)
            f0h = dbl_mb(k_r1+(i+0)-1)
            f1h_ = dbl_mb(k_r1+(i-1)-1)
            f2h_ = dbl_mb(k_r1+(i-2)-1)
            dif2 = (-f2h + 16*f1h - 30*f0h + 16*f1h_ - f2h_)
     &             / (12*de**2)
         endif

         r1 = dbl_mb(k_r1+(i+0)-1)

         et = dbl_mb(k_et+i-1)
         z = kappa * et
         zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
         call wfat_cchg(a, b, zz, chgm_)
         chgm = real(chgm_)
         r0 = d * z**(b/2) * exp(-z/2) * chgm
               
         lhs = dif2 +
     &        ((1-m**2)/(4*et**2) + beta0/et + erg/2) * r1
         rhs = -(beta1/et - mu(3)/2 + et/4) * r0
         epsi = lhs - rhs

         write(WF_R1TEST_UNIT, '(1x, i5, es15.6)', advance='no') i, z
         write(WF_R1TEST_UNIT, '(5es21.12)') lhs, rhs, epsi, r1, r0
      enddo
      close(WF_R1TEST_UNIT)
      
               
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_et))
     &     call errquit('WFAT: Cannot deallocate eta', 49, MA_ERR)
      if (.not. ma_free_heap(l_r1))
     &     call errquit('WFAT: Cannot deallocate r1', 49, MA_ERR)

      end
c=======================================================
