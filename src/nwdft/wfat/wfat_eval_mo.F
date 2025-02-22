c=======================================================
c=======================================================
      subroutine wfat_eval_mo_1pt(ghan, bhan, spin, maxder, vpt, scr,
     &           moval, d1moval, d2moval, d3moval)

      implicit none
#include "stdio.fh"
#include "bas.fh"
#include "global.fh"
#include "errquit.fh"
c     #include "cscf.fh"
#include "wfat_qm_param.fh"      
      integer ghan
      integer bhan
      integer spin
      integer maxder
      double precision vpt(3)
      double precision scr(*)    ! minimum no. of elements is nbas * (20 + nbas)
      double precision moval(*)
      double precision d1moval(3,*)
      double precision d2moval(6,*)
      double precision d3moval(10,*)

      integer nb0, nb1, nb2, nb3, ncf
      integer pb0, pb1, pb2, pb3, pcf
      integer nbas, i, j, id


      if (maxder > 3)
     &     call errquit('WFAT: Argument maxder of wfat_eval_mo_1pt '//
     &                  'cannot be larger than 3', 1, UERR)

      if (.not. bas_numbf(bhan, nbas))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)

      nb0 =    1 * nbas
      nb1 =    3 * nbas
      nb2 =    6 * nbas
      nb3 =   10 * nbas
      ncf = nbas * nbas

      pb0 = 1             !! basis functions
      pb1 = pb0 + nb0     !! 1st derivative of basis functions
      pb2 = pb1 + nb1     !! 2nd derivative of basis functions
      pb3 = pb2 + nb2     !! 3rd derivative of basis functions
      pcf = pb3 + nb3
      
      call wfat_eval_basis_1pt(ghan, bhan, maxder, vpt,
     &     scr(pb0), scr(pb1), scr(pb2), scr(pb3))


      do i = 1, nbas
         j = pcf + (i-1)*nbas
         call ga_get(g_mov(spin), 1, nbas, i, i, scr(j), nbas)
         if (maxder >= 0)
     &        moval(i) = sum(scr(pb0:pb0+nb0-1) * scr(j:j+nbas-1))
         if (maxder >= 1) then
            do id = 1, 3
               d1moval(id,i) = sum(scr(pb1:3:pb1+nb1-1) *
     &                             scr(j:j+nbas-1))
            enddo
         endif
         if (maxder >= 2) then
            do id = 1, 6
               d2moval(id,i) = sum(scr(pb2:6:pb2+nb2-1) *
     &                             scr(j:j+nbas-1))
            enddo
         endif
         if (maxder >= 3) then
            do id = 1, 10
               d3moval(id,i) = sum(scr(pb3:10:pb3+nb3-1) *
     &                             scr(j:j+nbas-1))
            enddo
         endif
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval_mo0_1pt(geom, bas, vpt, nv, g_v, val)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"      
      integer geom
      integer bas
      double precision vpt(3)
      integer nv
      integer g_v
      double precision val(nv)

      integer nb, l_bv, k_bv, l_mo, k_mo, i
      double precision dummy1, dummy2, dummy3


      if (.not. bas_numbf(bas, nb))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)

      if (.not. ma_alloc_get(mt_dbl, nb, 'basval', l_bv, k_bv))
     &     call errquit('WFAT: Cannot allocate basval', 1, MA_ERR)
      call wfat_eval_basis_npt(geom, bas, 0, 1, vpt,
     &     dbl_mb(k_bv), dummy1, dummy2, dummy3)

      if (.not. ma_alloc_get(mt_dbl, nb, 'MO', l_mo, k_mo))
     &     call errquit('WFAT: Cannot allocate MO', 1, MA_ERR)

      do i = 1, nv
         call ga_get(g_v, 1, nb, i, i, dbl_mb(k_mo), nb)
         val(i) = sum(dbl_mb(k_bv:k_bv+nb-1) * dbl_mb(k_mo:k_mo+nb-1))
      enddo

      if (.not. ma_free_heap(l_bv))
     &     call errquit('WFAT: Cannot deallocate basval', 1, MA_ERR)
      if (.not. ma_free_heap(l_mo))
     &     call errquit('WFAT: Cannot deallocate MO', 1, MA_ERR)

      end
c=======================================================      
