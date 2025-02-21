c=======================================================
c=======================================================
      subroutine wfat_test_chg(lmin_, lmax_, zq, erg)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer lmin_, lmax_
      double precision zq, erg

      integer l_m, k_m, l_u, k_u
      integer i, l, nl, nr
      parameter (nr=5)
      double precision kappa, lmb, a, aa, b, bb, r(nr)


      write(luout, '("TESTING FOR M(a,b,z) and U(a,b,z)")')

      nl = lmax_ - lmin_ + 1
      if (.not. MA_alloc_get(MT_DBL, nl, 'chgf', l_m, k_m))
     &     call errquit('WFAT: chgf allocation failed', 30, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl, 'chgu', l_u, k_u))
     &     call errquit('WFAT: chgu allocation failed', 31, MA_ERR)

      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0
      write(luout, '(1x)')
      write(luout, '(" erg = ", es14.6)') erg
      write(luout, '(" kappa = ", es14.6)') kappa
      r = (/ 0.005d0, 1.0d0, 4.8d0, 8.3d0, 240.7d0 /)

      do i = 1, nr
         write(luout, '(1x)')
         write(luout, '("------------------------------------")')
         write(luout, '("z = ", es14.6)') lmb * r(i)
         
         call wfat_chgm_1pt(r(i), lmin_, lmax_, a, b, lmb, dbl_mb(k_m))
         call wfat_chgu_1pt(r(i), lmin_, lmax_, a, b, lmb, dbl_mb(k_u))
         
         write(luout, '(1x)')
         write(luout, '(a14, a4, 2a14)') "a", "b", "M(a,b,z)",
     &        "U(a,b,z)"
         do l = lmin_, lmax_
            aa = l + a
            bb = 2*l + b
            write(luout, '(es14.6, i4, 2es14.6)') aa, int(bb),
     &           dbl_mb(k_m+l-lmin_), dbl_mb(k_u+l-lmin_)
         enddo

         write(luout, '("------------------------------------")')
      enddo


      if (.not. MA_free_heap(l_m))
     &     call errquit('WFAT: Cannot deallocate chgf', 27, MA_ERR)
      if (.not. MA_free_heap(l_u))
     &     call errquit('WFAT: Cannot deallocate chgu', 28, MA_ERR)


      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_test_chgu(a, b)

      implicit none
#include "mafdecls.fh"
#include "stdio.fh"      
#include "errquit.fh"
#include "wfat.fh"
      double precision a
      double precision b

      character*255 fname
      integer nz, i, md, l_u, k_u, l_uo, k_uo, l_m, k_m
      double precision z0, zn, dz, z, f2h, f1h, f0h, f1h_, f2h_
      double precision dif1, dif2, u, lhs, u_o, lhs_o


      call wfat_reset_lfqlib()
      if (a-1 >= 0.0d0) call wfat_update_lfqlib(a-1)
      

c     ==== Calculate U(a,b,z) ====
      nz = 1000
      if (.not. MA_alloc_get(MT_INT, nz, 'md', l_m, k_m))
     &     call errquit('WFAT: md allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nz, 'chgu', l_u, k_u))
     &     call errquit('WFAT: chgu allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nz, 'chgu_o', l_uo, k_uo))
     &     call errquit('WFAT: chgu_o allocation failed', 86, MA_ERR)
      z0 = 1.0d-3
      zn = 50.0d0
      dz = (zn - z0) / (nz - 1)
      do i = 1, nz
         z = z0 + (i-1)*dz
         call wfat_chgu(a, b, z, dbl_mb(k_u+i-1), int_mb(k_m+i-1))
         call wfat_chgu_OLD(a, b, z, dbl_mb(k_uo+i-1), md)
      enddo
      

      write(luout, '("TESTING FOR U(a,b,z)")')
      call util_file_name('test.whu', .false., .false., fname)
      open(WF_UTEST_UNIT, file=fname, action='write')
      write(WF_UTEST_UNIT, '("# nz = ", i4)') nz
      write(WF_UTEST_UNIT, '("# a = ", es15.6)') a
      write(WF_UTEST_UNIT, '("# b = ", es15.6)') b

      write(WF_UTEST_UNIT, '("#", a5, a21)', advance='no')
     &     'No.', 'z'
      write(WF_UTEST_UNIT, '(2a21, a3, 2a21)')
     &     'U', 'LHS', 'md', 'U_old', 'LHS_old'
      do i = 1, nz-2
         z = z0 + (i-1)*dz

c        ==== New code ====
         f2h = dbl_mb(k_u+(i+2)-1)
         f1h = dbl_mb(k_u+(i+1)-1)
         f0h = dbl_mb(k_u+(i+0)-1)
         f1h_ = dbl_mb(k_u+(i-1)-1)
         f2h_ = dbl_mb(k_u+(i-2)-1)
         dif2 = (-f2h + 16*f1h - 30*f0h + 16*f1h_ - f2h_)
     &          / (12*dz**2)
         dif1 = (-f2h + 8*f1h - 8*f1h_ + f2h_) / (12*dz)

         u = f0h
         lhs = z*dif2 + (b-z)*dif1 - a*u

         
c        ==== Old code ====
         f2h = dbl_mb(k_uo+(i+2)-1)
         f1h = dbl_mb(k_uo+(i+1)-1)
         f0h = dbl_mb(k_uo+(i+0)-1)
         f1h_ = dbl_mb(k_uo+(i-1)-1)
         f2h_ = dbl_mb(k_uo+(i-2)-1)
         dif2 = (-f2h + 16*f1h - 30*f0h + 16*f1h_ - f2h_)
     &          / (12*dz**2)
         dif1 = (-f2h + 8*f1h - 8*f1h_ + f2h_) / (12*dz)

         u_o = f0h
         lhs_o = z*dif2 + (b-z)*dif1 - a*u_o

         
c        ==== Print results ====
         write(WF_UTEST_UNIT, '(1x, i5, es21.12)', advance='no')
     &        i, z
         write(WF_UTEST_UNIT, '(2es21.12, i3, 2es21.12)')
     &        u, lhs, int_mb(k_m+i-1), u_o, lhs_o
      enddo
      close(WF_UTEST_UNIT)

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_m))
     &     call errquit('WFAT: Cannot deallocate md', 49, MA_ERR)
      if (.not. ma_free_heap(l_u))
     &     call errquit('WFAT: Cannot deallocate chgu', 49, MA_ERR)
      if (.not. ma_free_heap(l_uo))
     &     call errquit('WFAT: Cannot deallocate chgu_o', 49, MA_ERR)

      end
c=======================================================
      
      
