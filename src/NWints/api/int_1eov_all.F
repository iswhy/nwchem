c=======================================================
c=======================================================
c======== Imam =======      
      subroutine int_1eov_all(ibas, nbas, ovmat)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"      
      integer ibas                 ! [input]
      integer nbas                 ! [input]
      double precision ovmat(nbas, nbas)   ! [output]

      integer l_s, l_ov
      integer k_s, k_ov
      integer nbas0, maxs, maxov, ncont, nib, njb
      integer i, j, ib1, ib2, jb1, jb2      


c     ==== Initial checks ====
      if (.not. bas_numbf(ibas, nbas0)) then
         call errquit('int_1eov_all: Getting nbas0 failed',
     &        53, BASIS_ERR)
      else
         if (nbas0 /= nbas) then
            call errquit('int_1eov_all: Argument nbas is '//
     &           'inconsistent with handle ibas', 53,
     &           BASIS_ERR)
         endif
      endif
      
c     ==== Some allocations ====
      call int_mem_1e(maxov, maxs)
      if (.not. MA_alloc_get(MT_DBL, maxs, 'scr00', l_s, k_s))
     &  call errquit('int_1eov_all: 1e scratch allocation failed',
     &     31, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, maxov, '1eov_bl', l_ov, k_ov))
     &  call errquit('int_1eov_all: ovmat_bl allocation failed',
     &     32, MA_ERR)

      
c     ==== Begin filling ovmat ====
      ovmat = 0.0d0
      if (.not. bas_numcont(ibas, ncont))
     &     call errquit('int_1eov_all: Getting ncont failed',
     &          54, BASIS_ERR)

c     ==== Loop over column contractions ====
      do j = 1, ncont
         if (.not. bas_cn2bfr(ibas, j, jb1, jb2))
     &        call errquit('int_1eov_all: a cn2bfr call failed',
     &        57, BASIS_ERR)
         njb = jb2 - jb1 + 1

c        ==== Loop over row contractions ====         
         do i = 1, j
            if (.not. bas_cn2bfr(ibas, i, ib1, ib2))
     &           call errquit('int_1eov_all: a cn2bfr call failed',
     &           57, BASIS_ERR)
            nib = ib2 - ib1 + 1

c           ==== Calculate a block of ovmat ====
            call int_1eov(ibas, j, ibas, i, maxs, dbl_mb(k_s),     ! 1)
     &           maxov, dbl_mb(k_ov))                           
c            write(luout, '("block : ", 2i3)') i, j
c            call wfat_dprintmat(1, nib*njb, 'e', dbl_mb(k_ov))

c           ==== Fill the corresponding block of ovmat ====
            call int_1eov_fill(ib1, ib2, jb1, jb2, nbas,
     &           dbl_mb(k_ov), ovmat)
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_s))
     &     call errquit('int_1eov_all: scratch deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_ov))
     &     call errquit('int_1eov_all: 1eov_bl deallocation failed',
     &     4, MA_ERR)
      
      end
c     Notes
c     1) There doesn't seem to be a place where the ordering convention of the argument 'Ov' of the
c        subroutine int_1eov can be found. Such convention may be there but is not well accessible.
c        But experience of the author of the present subroutine indicates that the fastest moving 
c        elements of 'Ov' correspond to the basis functions within the 'jsh' contraction index (the 4th
c        argument of int_1eov).
c=======================================================


c=======================================================
c=======================================================
      subroutine int_1eov_fill(ib1, ib2, jb1, jb2, nbas,
     &     ovbl, ovmat)

      implicit none
      integer ib1, ib2                    ! [input], ranges of row of ovmat
      integer jb1, jb2                    ! [input], ranges of column of ovmat
      integer nbas                        ! [input], number of rows as well as columns of ovmat
      double precision ovbl(1:*)          ! [input], the block of ovmat corresponding to row ib1->ib2 and column jb1->jb2      1)
      double precision ovmat(nbas, nbas)  ! [output], array to be filled.

      integer jb, ib, k


      k = 0
      do jb = jb1, jb2         
         do ib = ib1, ib2
            k = k + 1
            
            ovmat(ib, jb) = ovbl(k)
            ovmat(jb, ib) = ovmat(ib, jb)
         enddo
      enddo
         
      end
c     Notes:
c     1) ovbl is a 1D array and it is arranged such that rows of ovmat for a
c        fixed column are adjacent in ovbl. This arrangement follows the apparent
c        ordering convention of the elements of the 1D output array 'Ov' returned
c        by the subroutine int_1eov.
c=======================================================
