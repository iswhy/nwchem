c=======================================================
c=======================================================
      subroutine wfat_test_o(n, m, erg, zq)

      implicit none
#include "stdio.fh"
#include "wfat.fh"
      integer n
      integer m
      double precision erg
      double precision zq

      character*255 fname
      integer i, neta
      double precision kp, a, b, d, g, o(2), eta0, etan, de, eta, ee(2)
      double precision r0o0, factorial
      external factorial
      

      kp = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kp
      b = dble(abs(m) + 1)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))
      
      write(luout, '("TESTING FOR o")')
      call util_file_name('test.wo', .false., .false., fname)
      open(WF_OTEST_UNIT, file=fname, action='write')
      write(WF_OTEST_UNIT, '("# n = ", i5)') n
      write(WF_OTEST_UNIT, '("# m = ", i5)') m
      write(WF_OTEST_UNIT, '("# energy = ", es15.6)') erg
      write(WF_OTEST_UNIT, '("# ion charge = ", es15.6)') zq
      
      call wfat_calc_o(n, m, erg, zq, o)
      write(WF_OTEST_UNIT, '("# o_\nu = ", 2es21.12)') o


      eta0 = 1.0d-12
      etan = 1.0d4
      neta = 10000
      de = (etan-eta0) / (neta-1)
      
      write(WF_OTEST_UNIT, '("#", a6, a15, 2a21, 5x, a21)') 
     &     'No.', 'eta', 'Integrand(1)', 'Integrand(2)',
     &     'R_\nu*O_\nu'
      do i = 1, neta
         eta = eta0 + (i-1)*de
         call wfat_calc_o_intd(eta, n, m, erg, zq, ee)
         call wfat_calc_r0o0(eta, n, m, erg, zq, r0o0)

         write(WF_OTEST_UNIT, '(i7, f15.6, 2es21.12, 5x, es21.12)') 
     &        i, eta, ee, r0o0
      enddo
      close(WF_OTEST_UNIT)
      
      end
c=======================================================
