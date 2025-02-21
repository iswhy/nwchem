c=======================================================
c=======================================================
      subroutine wfat_eval_vector(geom, bas, vpt, cf, val)

      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "errquit.fh"
      integer geom
      integer bas
      double precision vpt(3)
      double precision cf(1:*)
      double precision val

      integer nb, l_tm, k_tm, i


      if (.not. bas_numbf(bas, nb))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)

      if (.not. MA_alloc_get(mt_dbl, 20*nb, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: cannot allocate tmp.',0, MA_ERR)
      
      call wfat_eval_basis_1pt(geom, bas, 0, vpt, dbl_mb(k_tm),
     &     dbl_mb(k_tm+nb), dbl_mb(k_tm+3*nb), dbl_mb(k_tm+10*nb))

      val = 0.0d0
      do i = 1, nb
         val = val + cf(i) * dbl_mb(k_tm+i-1)
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 90, MA_ERR)
      
      end
c=======================================================
      
