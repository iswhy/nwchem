c=======================================================
c=======================================================
c======== Imam =======      
      subroutine int_1epe_all(ibas, nbas, pemat)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"      
      integer ibas                 ! [input]
      integer nbas                 ! [input]
      double precision pemat(nbas, nbas)   ! [output]

      integer l_s, l_pe
      integer k_s, k_pe
      integer nbas0, maxs, maxpe, ncont, nib, njb
      integer i, j, ib1, ib2, jb1, jb2      


c     ==== Initial checks ====
      if (.not. bas_numbf(ibas, nbas0)) then
         call errquit('int_1epe_all: Getting nbas0 failed',
     &        53, BASIS_ERR)
      else
         if (nbas0 /= nbas) then
            call errquit('int_1epe_all: Argument nbas is '//
     &           'inconsistent with handle ibas', 53,
     &           BASIS_ERR)
         endif
      endif
      
c     ==== Some allocations ====
      call int_mem_1e(maxpe, maxs)
      if (.not. MA_alloc_get(MT_DBL, maxs, 'scr00', l_s, k_s))
     &  call errquit('int_1epe_all: 1e scratch allocation failed',
     &     31, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, maxpe, '1epe_bl', l_pe, k_pe))
     &  call errquit('int_1epe_all: pemat_bl allocation failed',
     &     32, MA_ERR)

      
c     ==== Begin filling pemat ====
      pemat = 0.0d0
      if (.not. bas_numcont(ibas, ncont))
     &     call errquit('int_1epe_all: Getting ncont failed',
     &          54, BASIS_ERR)

c     ==== Loop over column contractions ====
      do j = 1, ncont
         if (.not. bas_cn2bfr(ibas, j, jb1, jb2))
     &        call errquit('int_1epe_all: a cn2bfr call failed',
     &        57, BASIS_ERR)
         njb = jb2 - jb1 + 1

c        ==== Loop over row contractions ====         
         do i = 1, j
            if (.not. bas_cn2bfr(ibas, i, ib1, ib2))
     &           call errquit('int_1epe_all: a cn2bfr call failed',
     &           57, BASIS_ERR)
            nib = ib2 - ib1 + 1

c           ==== Calculate a block of pemat ====
            call int_1epe(ibas, j, ibas, i, maxs, dbl_mb(k_s),     ! 1)
     &           maxpe, dbl_mb(k_pe))                           
c            write(luout, '("block : ", 2i3)') i, j
c            call wfat_dprintmat(1, nib*njb, 'e', dbl_mb(k_pe))

c           ==== Fill the corresponding block of pemat ====
            call int_1epe_fill(ib1, ib2, jb1, jb2, nbas,
     &           dbl_mb(k_pe), pemat)
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_s))
     &     call errquit('int_1epe_all: scratch deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_pe))
     &     call errquit('int_1epe_all: 1epe_bl deallocation failed',
     &     4, MA_ERR)
      
      end
c     Notes
c     1) There doesn't seem to be a place where the ordering convention of the argument 'Pe' of the
c        subroutine int_1epe can be found. Such convention may be there but is not well accessible.
c        But experience of the author of the present subroutine indicates that the fastest moving 
c        elements of 'Pe' correspond to the basis functions within the 'jsh' contraction index (the 4th
c        argument of int_1epe).
c=======================================================


c=======================================================
c=======================================================
      subroutine int_1epe_fill(ib1, ib2, jb1, jb2, nbas,
     &     pebl, pemat)

      implicit none
      integer ib1, ib2                    ! [input], ranges of row of pemat
      integer jb1, jb2                    ! [input], ranges of column of pemat
      integer nbas                        ! [input], number of rows as well as columns of pemat
      double precision pebl(1:*)          ! [input], the block of pemat corresponding to row ib1->ib2 and column jb1->jb2      1)
      double precision pemat(nbas, nbas)  ! [output], array to be filled.

      integer jb, ib, k


      k = 0
      do jb = jb1, jb2         
         do ib = ib1, ib2
            k = k + 1
            
            pemat(ib, jb) = pebl(k)
            pemat(jb, ib) = pemat(ib, jb)
         enddo
      enddo
         
      end
c     Notes:
c     1) pebl is a 1D array and it is arranged such that rows of pemat for a
c        fixed column are adjacent in pebl. This arrangement follows the apparent
c        ordering convention of the elements of the 1D output array 'Pe' returned
c        by the subroutine int_1epe.
c=======================================================
