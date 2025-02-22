c=======================================================
c=======================================================
c======== Imam =======      
      subroutine int_mpole_all(ibas, nbas, lmax, center, mpmat)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"      
      integer ibas                 ! [input]
      integer nbas                 ! [input]
      integer lmax
      double precision center(1:3) ! [Input] coordinate of the origin for the multipole.
      double precision mpmat(nbas, nbas, 1:*)   ! [output]

      integer l_s, l_mp
      integer k_s, k_mp
      integer nbas0, maxs, maxmp, maxq, ncont, nib, njb
      integer i, j, iq, ib1, ib2, jb1, jb2      


c     ==== Initial checks ====
      if (.not. bas_numbf(ibas, nbas0)) then
         call errquit('int_mpole_all: Getting nbas0 failed',
     &        53, BASIS_ERR)
      else
         if (nbas0 /= nbas) then
            call errquit('int_mpole_all: Argument nbas is '//
     &           'inconsistent with handle ibas', 53,
     &           BASIS_ERR)
         endif
      endif
      
c     ==== Some allocations ====
      call int_mem_dipole(maxmp, maxs, ibas, ibas, lmax)
      if (.not. MA_alloc_get(MT_DBL, maxs, 'scr00', l_s, k_s))
     &  call errquit('int_mpole_all: 1e scratch allocation failed',
     &     31, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, maxmp, 'mpole_bl', l_mp, k_mp))
     &  call errquit('int_mpole_all: mpmat_bl allocation failed',
     &     32, MA_ERR)

      
c     ==== Begin filling mpmat ====
      if (.not. bas_numcont(ibas, ncont))
     &     call errquit('int_mpole_all: Getting ncont failed',
     &          54, BASIS_ERR)
      maxq = ((lmax+1) * (lmax+2) * (lmax+3)) / 6

      do j = 1, ncont
         if (.not. bas_cn2bfr(ibas, j, jb1, jb2))
     &        call errquit('int_mpole_all: a cn2bfr call failed',
     &        57, BASIS_ERR)
         njb = jb2 - jb1 + 1
         
         do i = 1, j
            if (.not. bas_cn2bfr(ibas, i, ib1, ib2))
     &           call errquit('int_mpole_all: a cn2bfr call failed',
     &           57, BASIS_ERR)
            nib = ib2 - ib1 + 1

c           ==== Calculate a block of mpmat ====
            call int_mpole(ibas, j, ibas, i, lmax, center, maxs,
     &           dbl_mb(k_s), maxmp, dbl_mb(k_mp))
c            write(luout, '("block : ", 2i3)') i, j
c            call wfat_dprintmat(1, nib*njb*maxq, 'e', dbl_mb(k_mp))

c           ==== Fill the corresponding block of mpmat ====
            call int_mpole_fill(ib1, ib2, jb1, jb2, nbas, lmax,
     &           dbl_mb(k_mp), mpmat)            
         enddo
      enddo

c     ==== Apply the Hermitian property ====
      do iq = 1, maxq
         do j = 1, nbas
            do i = 1, j
               mpmat(j,i,iq) = mpmat(i,j,iq)
            enddo
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_s))
     &     call errquit('int_mpole_all: scratch deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_mp))
     &     call errquit('int_mpole_all: mpole_bl deallocation failed',
     &     4, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine int_mpole_fill(ib1, ib2, jb1, jb2, nbas, lmax,
     &           mpbl, mpmat)

      implicit none

      integer ib1, ib2
      integer jb1, jb2
      integer nbas
      integer lmax
      double precision mpbl(1:*)
      double precision mpmat(nbas, nbas, 1:*)

      integer maxq, nq, nq_
      integer iq, iiq, jb, jb_, ib, ib_, k, l

      
      k = 0
      do l = 0, lmax          !! 1)         
         nq = ((l+1)*(l+2))/2
         nq_ = ((l+0)*(l+1)*(l+2)) / 6

         do jb = jb1, jb2
            jb_ = jb - jb1 + 1
            do iq = 1, nq
               do ib = ib1, ib2
                  ib_ = ib - ib1 + 1
                  k = k + 1
                  
                  iiq = nq_ + iq
                  mpmat(ib, jb, iiq) = mpbl(k)
               enddo
            enddo
         enddo
         
      enddo

      end
c  Notes:
c    1) The indexing of the one-dimensional array mpbl is such that for a given
c       l, the fastest moving index is ib_, then iq, and finally jb_.
c=======================================================
