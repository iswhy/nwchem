c=======================================================
c=======================================================
      subroutine wfat_overlap_numinteg_xy(geom, bas, nx, xlim, ny, ylim,
     &     z, ovl)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "bas.fh"
#include "stdio.fh"
#include "global.fh"      
      integer geom
      integer bas
      integer nx
      double precision xlim(2)
      integer ny
      double precision ylim(2)
      double precision z
      double precision ovl(*)

      integer nb, nmat, npt
      integer ix, iy, i_bv
      integer l_bv, l_ov0, k_bv, k_ov0
      double precision x, y, dx, dy, vpt(3), dummy1, dummy2, dummy3

     
      if (.not. bas_numbf(bas, nb))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)
      nmat = nb * nb


      if (.not. ma_alloc_get(mt_dbl, nb, 'basval', l_bv, k_bv))
     &     call errquit('WFAT: Cannot allocate basval', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nmat, 'ovl0', l_ov0, k_ov0))
     &     call errquit('WFAT: Cannot allocate ovl0', 1, MA_ERR)
      
      
c     ==== xy-integrate overlap matrix at z ====
      dx = (xlim(2) - xlim(1)) / (nx - 1)
      dy = (ylim(2) - ylim(1)) / (ny - 1)
      ovl(1:nmat) = 0.0d0      
      do iy = 1, ny
         y = ylim(1) + dy*(iy-1)
         do ix = 1, nx
            x = xlim(1) + dx*(ix-1)
            
            vpt = (/ x, y, z /)
            call wfat_eval_basis_npt(geom, bas, 0, 1, vpt,
     &           dbl_mb(k_bv), dummy1, dummy2, dummy3)            
            call wfat_project_mat(nb, 1, dbl_mb(k_bv), dbl_mb(k_ov0))

            ovl(1:nmat) = ovl(1:nmat) + dbl_mb(k_ov0:k_ov0+nmat-1)
         enddo
      enddo
      ovl(1:nmat) = ovl(1:nmat) * dx*dy


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_bv))
     &     call errquit('WFAT: Cannot deallocate basval', 1, MA_ERR)
      if (.not. ma_free_heap(l_ov0))
     &     call errquit('WFAT: Cannot deallocate ovl0', 1, MA_ERR)
      
      end
c=======================================================
      

cOLD c=======================================================
cOLD c=======================================================
cOLD       subroutine wfat_overlap_numinteg_xy(geom, bas, nx, xlim, ny, ylim,
cOLD      &     z, ovl)
cOLD 
cOLD       implicit none
cOLD #include "mafdecls.fh"
cOLD #include "errquit.fh"
cOLD #include "bas.fh"
cOLD #include "stdio.fh"
cOLD #include "global.fh"      
cOLD       integer geom
cOLD       integer bas
cOLD       integer nx
cOLD       double precision xlim(2)
cOLD       integer ny
cOLD       double precision ylim(2)
cOLD       double precision z
cOLD       double precision ovl(*)
cOLD 
cOLD       integer nb, nmat, npt
cOLD       integer i, ix, iy, i_bv, i_pt1, i_pt2
cOLD       integer l_pt, l_bv, l_ov0, k_pt, k_bv, k_ov0
cOLD       double precision x, y, dx, dy, dummy1, dummy2, dummy3
cOLD 
cOLD      
cOLD       dx = (xlim(2) - xlim(1)) / (nx - 1)
cOLD       dy = (ylim(2) - ylim(1)) / (ny - 1)
cOLD       if (.not. bas_numbf(bas, nb))
cOLD      &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)
cOLD       nmat = nb * nb
cOLD 
cOLD 
cOLD       ovl(1:nmat) = 0.0d0
cOLD 
cOLD       npt = nx * ny
cOLD       if (.not. ma_alloc_get(mt_dbl, 3*npt, 'points', l_pt, k_pt))
cOLD      &     call errquit('WFAT: Cannot allocate points', 1, MA_ERR)
cOLD 
cOLD       
cOLD c     ==== Construct point array containing (x,y,z) triplet ====
cOLD       i = 0
cOLD       do iy = 1, ny
cOLD          y = ylim(1) + dy*(iy-1)
cOLD          do ix = 1, nx
cOLD             x = xlim(1) + dx*(ix-1)
cOLD 
cOLD             i = i + 1                  ! count the point
cOLD             i_pt1 = k_pt + (i-1)*3     ! x
cOLD             i_pt2 = i_pt1 + 2          ! z
cOLD             dbl_mb(i_pt1:i_pt2) = (/ x, y, z /)
cOLD          enddo
cOLD       enddo
cOLD 
cOLD 
cOLD c     ==== Calculate basis functions at the above points ====
cOLD       if (.not. ma_alloc_get(mt_dbl, npt*nb, 'basval', l_bv, k_bv))
cOLD      &     call errquit('WFAT: Cannot allocate basval', 1, MA_ERR)
cOLD       call wfat_eval_basis_npt(geom, bas, 0, npt, dbl_mb(k_pt),
cOLD      &     dbl_mb(k_bv), dummy1, dummy2, dummy3)
cOLD c      write(luout, '("here b1")')
cOLD 
cOLD       
cOLD c     ==== Transpose basval ====
cOLD       call wfat_dtranspose(npt, nb, dbl_mb(k_bv), dbl_mb(k_bv))   ! Now, basval is conceptually an [nb] by [npt] array. 
cOLD c      write(luout, '("here b2")')      
cOLD       
cOLD c     ==== xy-integrate overlap matrix at z ====
cOLD       if (.not. ma_alloc_get(mt_dbl, nmat, 'ovl0', l_ov0, k_ov0))
cOLD      &     call errquit('WFAT: Cannot allocate ovl0', 1, MA_ERR)
cOLD       i = 0
cOLD       do iy = 1, ny
cOLD          do ix = 1, nx
cOLD             i = i + 1
cOLD             i_bv = k_bv + (i-1)*nb
cOLD             call wfat_project_mat(nb, 1, dbl_mb(i_bv), dbl_mb(k_ov0))
cOLD 
cOLD             ovl(1:nmat) = ovl(1:nmat) + dbl_mb(k_ov0:k_ov0+nmat-1)
cOLD          enddo
cOLD       enddo
cOLD c      write(luout, '("here b3")')      
cOLD       ovl(1:nmat) = ovl(1:nmat) * dx*dy
cOLD c      write(luout, '("here b4")')
cOLD 
cOLD c     ==== Deallocations ====
cOLD       if (.not. ma_free_heap(l_pt))
cOLD      &     call errquit('WFAT: Cannot deallocate points', 1, MA_ERR)
cOLD       if (.not. ma_free_heap(l_bv))
cOLD      &     call errquit('WFAT: Cannot deallocate basval', 1, MA_ERR)
cOLD       if (.not. ma_free_heap(l_ov0))
cOLD      &     call errquit('WFAT: Cannot deallocate ovl0', 1, MA_ERR)
cOLD       
cOLD       end
cOLD c=======================================================
      
      
