!==================================================
!==================================================
      subroutine wfat_gauntp_lscan(l1, m1, l2, m2, gn)
c     DESCRIPTION:
c     Calculates Gaunt's integral defined as
c     $$
c       I(l1,m2,l2,m2,l3) = \int_{-1}^1 dx P_{l1,m1}(x) P_{l2,m2}(x) P_{l3,m1+m2}(x)
c     $$
c     for $l3=[l_min, l1+l2]$, where $l_min=max(|l1-l2|,|m1+m2|)$ and $P_{lm}(x)$ 
c     is a NORMALIZED associated Legendre function, that is,
c     $$
c       \int_{-1}^1 dx (P_{lm}(x))^2 = 1 .
c     $$
      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer l1                ! [input], l1
      integer m1                ! [input], m1
      integer l2                ! [input], l2
      integer m2                ! [input], m2
      double precision gn(1:*)  ! [output], Gaunt's integrals, the minimum size should be l1+l2 - max(abs(l1-l2),abs(m1+m2)) + 1.

      integer l_min, n, n0, l_3j, k_3j, l_3j0, k_3j0, i, l3, i0
      double precision f


      if (l1<0 .or. l2<0 .or. l3<0)
     &     call errquit('WFAT: The argument l1, l2, and l3 of '//
     &     'wfat_gauntp_lscan and wfat_gaunty_lscan must be zero '//
     &     'or positive integers.', 63, UERR)

      
c     / l1  l2      l3   \
c     \ m1  m2  -(m1+m2) /
c     where l3=[l_min,l1+l2].
      l_min = max(abs(l1-l2), abs(m1+m2))
      n = (l1+l2) - l_min + 1
      if (.not. MA_alloc_get(MT_DBL, n, '3j', l_3j, k_3j))
     &     call errquit('WFAT: 3j allocation failed', 59, MA_ERR)
      call wfat_3j_jscan(2*l1, 2*m1, 2*l2, 2*m2, dbl_mb(k_3j))

      
c     / l1  l2  l3 \
c     \  0   0   0 /
c     where l3 = [|l1-l2|,l1+l2].
      n0 = 2*min(l1,l2) + 1
      if (.not. MA_alloc_get(MT_DBL, n0, '3j_000', l_3j0, k_3j0))
     &     call errquit('WFAT: 3j_000 allocation failed', 59, MA_ERR)
      call wfat_3j_jscan000(l1, l2, dbl_mb(k_3j0))

      
      gn(1:n) = 0.0d0
      do i = 1, n
         l3 = l_min + i - 1
         i0 = i + l_min - abs(l1-l2)
         if (mod(l1+l2+l3,2) == 0) then
            f = (2*l1+1) * (2*l2+1) * (2*l3+1) / 2.0d0
            gn(i) = dbl_mb(k_3j0+i0-1) * dbl_mb(k_3j+i-1)
            gn(i) = (-1)**(m1+m2) * sqrt(f) * gn(i)
         endif
      enddo


      if (.not. MA_free_heap(l_3j))
     &     call errquit('WFAT: 3j deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_3j0))
     &     call errquit('WFAT: 3j000 deallocation failed', 61, MA_ERR)
      end
!==================================================


!==================================================
!==================================================
      subroutine wfat_gaunty_lscan(l1, m1, l2, m2, gn)
c     DESCRIPTION:
c     Calculates Gaunt's integral defined as
c     $$
c       I(l1,m2,l2,m2,l3,m3) = \int_{-1}^1 dx Y_{l1,m1}(x) Y_{l2,m2}(x) Y_{l3,m3}(x)
c     $$
c     for $l3=[l_min, l1+l2]$, where $l_min=max(|l1-l2|,|m1+m2|)$.
      implicit none
      integer l1                ! [input], l1
      integer m1                ! [input], m1
      integer l2                ! [input], l2
      integer m2                ! [input], m2
      double precision gn(1:*)  ! [output], Gaunt's integrals, the minimum size should be l1+l2 - max(abs(l1-l2),abs(m1+m2)) + 1.

      integer l_min, n, i
      double precision pi_const
      parameter (pi_const=3.14159265359d0)

      call wfat_gauntp_lscan(l1, m1, l2, m2, gn)

      l_min = max(abs(l1-l2), abs(m1+m2))
      n = (l1+l2) - l_min + 1
      do i = 1, n
         gn(i) = gn(i)/(-1)**(m1+m2)/sqrt(2*pi_const)
      enddo

      end
!==================================================
      
