c=======================================================
c=======================================================
      subroutine wfat_test_erfc(c, a, t, xx, n)

      implicit none
#include "wfat.fh"
#include "stdio.fh"
      double precision c
      double precision a
      double precision t
      double precision xx
      integer n

      character*255 fname
      integer i
      double precision dx, x, f

      dx = xx/(n-1)

      write(luout, '("TESTING FOR ERFC")')
      write(luout, '("  c = ", es14.6)') c
      write(luout, '("  a = ", es14.6)') a
      write(luout, '("  t = ", es14.6)') t
      write(luout, '("  xx = ", es14.6)') xx
      write(luout, '("  n = ", i6)') n
      call util_file_name('test.werfc', .false., .false., fname)
      open(WF_ERFCTEST_UNIT, file=fname, action='write')

      do i = 1, n
         x = 0.0d0 + (i-1)*dx
         f = c * erfc(a*(x-t))
         write(WF_ERFCTEST_UNIT, '(i6, 2es25.12E3)') i, x, f
      enddo

      close (WF_ERFCTEST_UNIT)

      end
c=======================================================
