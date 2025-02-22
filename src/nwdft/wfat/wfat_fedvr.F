c==========================================================
c==========================================================
      subroutine wfat_fedvr_rule(ni, n, maxn, qnodes, xq, wq)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
      integer ni                     ! [input].
      integer n(ni)                  ! [input].
      integer maxn                   ! [input].
      double precision qnodes(ni)    ! [input].
      double precision xq(maxn,ni)   ! [output].
      double precision wq(maxn,ni)   ! [output].

      integer i
      double precision x1, x2

      if (minval(n) < 2) then
         write(luout, '("min. val of n = ", i3)') minval(n)
         call errquit('The smallest element of the array argument n '//
     &        'of wfat_fedvr_rule must be bigger than or equal to 2.',
     &        73, UERR)
      endif
      
      do i = 1, ni
c        ==== Obtain the grids and weights ====
         if (i == 1) then
            x1 = 0.0d0
         else
            x1 = qnodes(i-1)
         endif
         x2 = qnodes(i)
         call wfat_lgl(x1, x2, n(i)-1, xq(:,i), wq(:,i))
      enddo

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_fedvr_basnodes(ni, n, maxn, xq, x)

      implicit none
#include "stdio.fh"
      integer ni                     ! [input].
      integer n(ni)                  ! [input].
      integer maxn                   ! [input].
      double precision xq(maxn,ni)   ! [input].
      double precision x(1:*)        ! [output].

      logical incl
      integer j, i, m

c     ==== Loop over ket ====
      j = 0
      do i = 1, ni
         do m = 1, n(i)

c           == Exclude (1,1) (the 1st basis), (ni,n(ni)) (the last ==
c           == basis), and (1,:) (have been absorbed into bridges) ==
            incl = .true.
            if ((i==1 .and. m==1) .or. (i==ni .and. m==n(ni))
     &           .or. m==1) incl = .false.
            if (incl) j = j + 1

            if (incl) x(j) = xq(m,i)

c            write(luout, '(i5)') j
         enddo
      enddo

      end
c==========================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_get_fedvrnr(nfe, ndvr, nr)

      implicit none
      integer nfe
      integer ndvr(nfe)
      integer nr

      nr = (sum(ndvr) - nfe*2) + (nfe - 1)

      end
c=======================================================
      
