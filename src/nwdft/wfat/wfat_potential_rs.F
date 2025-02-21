c=======================================================
c=======================================================
      subroutine wfat_1ehyd_erf_all(bas, nb, mu, rcent, crg, vhyd)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"      
      integer bas                     ! [input]
      integer nb                      ! [input]
      double precision mu             ! [input]
      double precision rcent(1:3)     ! [Input] coordinate of potential center
      double precision crg            ! [Input] charge of potential center
      double precision vhyd(nb, nb)   ! [output]

      integer l_s, l_hb
      integer k_s, k_hb
      integer nb0, leri, lscr, ncont, nib, njb, nblock
      integer i, j, ib1, ib2, jb1, jb2



c     ==== Initial checks ====
      if (.not. bas_numbf(bas, nb0)) then
         call errquit('int_1ehyd_all: Getting nb0 failed',
     &        53, BASIS_ERR)
      else
         if (nb0 /= nb) then
            call errquit('int_1ehyd_all: Argument nb is '//
     &           'inconsistent with handle bas', 53,
     &           BASIS_ERR)
         endif
      endif
      
c     ==== Some allocations ====
      call int_mem_2e3c(leri, lscr)      
cOLD      call int_mem_1e(max1e, mem1)
      if (.not. MA_alloc_get(MT_DBL, lscr, 'scr00', l_s, k_s))
     &     call errquit('int_1ehyd_all: 1e scratch allocation failed',
     &     31, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, leri, 'vhyd_bl', l_hb, k_hb))
     &     call errquit('int_1ehyd_all: vhyd_bl allocation failed',
     &     32, MA_ERR)

c     ==== Begin filling vhyd ====
      if (.not. bas_numcont(bas, ncont))
     &     call errquit('int_1ehyd_all: Getting ncont failed',
     &          54, BASIS_ERR)

      vhyd = 0.0d0
      do j = 1, ncont
         if (.not. bas_cn2bfr(bas, j, jb1, jb2))
     &        call errquit('int_1ehyd_all: a cn2bfr call failed',
     &        57, BASIS_ERR)
         njb = jb2 - jb1 + 1
         
         do i = 1, j
            if (.not. bas_cn2bfr(bas, i, ib1, ib2))
     &           call errquit('int_1ehyd_all: a cn2bfr call failed',
     &           57, BASIS_ERR)
            nib = ib2 - ib1 + 1

c           ==== Calculate a block of vhyd ====
            dbl_mb(k_s : k_s+lscr-1) = 0.0d0            
            call wfat_1ehyd_erf(bas, j, i, mu, rcent, crg, 
     &           lscr, dbl_mb(k_s), leri, dbl_mb(k_hb))

c           ==== Fill the corresponding block of vhyd ====
            nblock = nib * njb
            call wfat_1ehyd_erf_fill(nib, njb, 
     &           dbl_mb(k_hb:k_hb+nblock-1), vhyd(ib1:ib2, jb1:jb2))
         enddo
      enddo

c     ==== Apply the Hermitian property ====
      do j = 1, nb
         do i = 1, j
            vhyd(j,i) = vhyd(i,j)
         enddo
      enddo

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_s))
     &     call errquit('int_1ehyd_all: scratch deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_hb))
     &     call errquit('int_1ehyd_all: vhyd_bl deallocation failed',
     &     4, MA_ERR)

      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_1ehyd_erf(bas, ish, jsh, mu, vpt, crg, 
     &           lscr, scr, leri, eri)

      implicit none
#include "errquit.fh"      
#include "bas.fh"
#include "geom.fh"
#include "wfat_rs.fh"
      integer bas
      integer ish
      integer jsh
      double precision mu
      double precision vpt(3)
      double precision crg
      integer lscr
      double precision scr(1:*)
      integer leri
      double precision eri(1:*)

      double precision nrm, pi_const
      parameter (pi_const=3.14159265359d0)
      
      logical int_normalize
      external int_normalize


      if (.not. bas_check_handle(rsbas, 'rsbas basis must be created.'))
     &     call errquit('WFAT: Error rsbas basis handle has not yet '//
     &     'been created upon calling wfat_1ehyd_erf.', 1, BASIS_ERR)
      

      if (.not. geom_cent_set(rsgeom, 1, rstag, vpt, rschg))
     &     call errquit('WFAT: Cannot set RS center', 1, GEOM_ERR)
      if (.not. bas_set_exponent(rsbas, 1, mu**2, 1))
     &     call errquit('WFAT: Cannot set RS exponent', 1, BASIS_ERR)
      nrm = crg * mu**3 / sqrt(pi_const) / pi_const
      if (.not. bas_set_coeff(rsbas, 1, nrm, 1))
     &     call errquit('WFAT: Cannot set RS coeff', 1, BASIS_ERR)
      

      call wfat_int_2e3c(rsbas, 1, bas, ish, jsh, lscr, scr, leri,
     &     .false., eri)
cnope      call int_2e3c(rsbas, 1, bas, jsh, ish, lscr, scr, leri, eri)      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_1ehyd_erf_fill(m, n, vhbl, vh)

      implicit none
      integer m                   ! [input]
      integer n                   ! [input]
      double precision vhbl(m,n)  ! [input] 
      double precision vh(m,n)    ! [output]

      
      vh = vhbl
      
      end
c=======================================================
      
      
