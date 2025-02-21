c=================================================
c=================================================
      subroutine wfat_print_orb1(rtdb, nch, orbid, spinid)

      implicit none
#include "global.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "wfat_qm_param.fh"
#include "wfat.fh"
      integer rtdb
      integer nch
      integer orbid(nch)
      integer spinid(nch)

      character*256 name0
      character*11 idstring
      integer np(3), ncf, l_cf, k_cf, i_cf0, i_cf1, i, j, k1, k2, is, id
      double precision wfn(nbas0), erg(nch), pb(3,3,nch), ef, xr(2)
      double precision yr(2), zr(2), scr(2), bt, gm, protmat(3,3), cu(3)


      if (ga_nodeid() == 0)
     &     write(luout, '("Printing first-order corrected orbitals '//
     &     'to cube files")') 
      

      if (.not. rtdb_get(rtdb, 'wfat:corb_orient', mt_dbl, 2, scr))
     &     call errquit('WFAT: Failed reading erg from RTDB.',
     &     1, RTDB_ERR)
      bt = scr(1)
      gm = scr(2)
      if (.not. rtdb_get(rtdb, 'wfat:corb_field', mt_dbl, 1, ef))
     &     call errquit('WFAT: Failed reading corb_field from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:corb_np', mt_int, 3, np))
     &     call errquit('WFAT: Failed reading corb_np from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:corb_xr', mt_dbl, 2, xr))
     &     call errquit('WFAT: Failed reading corb_xr from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:corb_yr', mt_dbl, 2, yr))
     &     call errquit('WFAT: Failed reading corb_yr from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:corb_zr', mt_dbl, 2, zr))
     &     call errquit('WFAT: Failed reading corb_zr from RTDB.',
     &     1, RTDB_ERR)
      

      call wfat_passrot(bt, gm, protmat)
      cu = protmat(:,3)

      
      ncf = nbas0 * (nch + 3*nch)
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vch', l_cf, k_cf))
     &     call errquit('WFAT: vch allocation failed', 28, MA_ERR)
      i_cf0 = k_cf
      call wfat_get_coef(nch, orbid, spinid, dbl_mb(i_cf0), erg)
      i_cf1 = k_cf + nbas0*nch
      call wfat_orb1_property(nch, spinid, orbid, dbl_mb(i_cf1), pb)


      do i = 1, nch
         wfn = 0.0d0
         do j = 1, 3
            k1 = i_cf1 + nbas0*3*(i-1) + nbas0*(j-1)
            k2 = k1 + nbas0 - 1
            wfn = wfn + cu(j)*dbl_mb(k1:k2)
         enddo

         k1 = i_cf0 + nbas0*(i-1)
         k2 = k1 + nbas0 - 1
         wfn = dbl_mb(k1:k2) + ef*wfn

         is = spinid(i)
         id = orbid(i) - nel(is)
         write(idstring, '(sp, i10)') id
         if (is == 1) then
            name0 = 'corb.a.'//trim(adjustl(idstring))
         else if (is == 2) then
            name0 = 'corb.b.'//trim(adjustl(idstring))
         endif
         call wfat_print_cube(name0, WF_CORB_UNIT, np, 
     &        xr, yr, zr, geom_hd, bas_hd, wfn, .True.)
      enddo
      

      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: vch deallocation failed', 28, MA_ERR)

      end
c=================================================
