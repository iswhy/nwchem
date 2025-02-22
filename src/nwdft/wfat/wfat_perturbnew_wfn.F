c=======================================================
c=======================================================
      subroutine wfat_orb1_property(nch, spinid, orbid, vch1, pb)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"      
#include "wfat_qm_param.fh"
#include "global.fh"
#include "stdio.fh"
      integer nch
      integer spinid(nch)
      integer orbid(nch)
      double precision vch1(nbas0, 3, nch)   ! [output], 
      double precision pb(3,3,nch)           ! [output], 

      integer i, is, imo, io, k, j, m, id0, idm, ieg0, iegm, i_dpo
      integer is1, is2, l_mo, k_mo, l_dpo, k_dpo, ndid
      double precision derg, ctr(3), vpt(3), val

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      ndid = ga_nodeid()

c     ==== Calculate dipole matrix in the MO basis ====      
      if (.not. MA_alloc_get(MT_DBL, nbas0**2*nspin0*3, 'dpmo',
     &     l_dpo, k_dpo)) 
     &     call errquit('WFAT: dpmo allocation failed', 19, MA_ERR)
      ctr = (/ 0.0d0, 0.0d0, 0.0d0 /)                !! 1)
      call wfat_calc_dpmo(ctr, dbl_mb(k_dpo))

      
c     ==== Compute the static polarizabilities in molecular frame ====
      call wfat_calc_polb(nch, orbid, spinid, dbl_mb(k_dpo), pb)
      if (ndid == 0) then
         call wfat_print_polb_tensors(nch, pb)
      endif
      call ga_sync()
      

c     ==== MO coefficients ====
      if (.not. MA_alloc_get(MT_DBL, nbas0**2*nspin0, 'mo',
     &     l_mo, k_mo)) 
     &     call errquit('WFAT: MO allocation failed', 19, MA_ERR)
      do is = 1, nspin0
         imo = k_mo + (is-1)*nbas0*nbas0
         call ga_get(g_mov(is), 1, nbas0, 1, nbas0, dbl_mb(imo), nbas0)
      enddo
      
         
c     ==== Calculate the GTO coefficients of ====
c     ==== the first order correction terms  ====
      do i = 1, nch             ! Loop over ionized channels
         io = orbid(i)
         is = spinid(i)
         do k = 1, 3            ! Loop over dipole components

            do j = 1, nbas0     ! Loop over GTO bases
               vch1(j,k,i) = 0.0d0
               do m = 1, nbas0  ! Sum over MO's
                  id0 = int_mb(k_modgid + (is-1)*nbas0 + (io-1))
                  idm = int_mb(k_modgid + (is-1)*nbas0 + (m-1))
                  if (id0 /= idm) then
                     ieg0 = k_moerg + (is-1)*nbas0 + (io-1)
                     iegm = k_moerg + (is-1)*nbas0 + (m-1)
                     derg = dbl_mb(ieg0) - dbl_mb(iegm)
                  
                     imo = k_mo + (is-1)*nbas0*nbas0 + (m-1)*nbas0 +
     &                     (j-1)
                     i_dpo = k_dpo + nbas0**2*nspin0*(k-1) +
     &                       nbas0**2*(is-1) + nbas0*(io-1) + (m-1)
                     vch1(j,k,i) = vch1(j,k,i) +
     &                             dbl_mb(imo) * dbl_mb(i_dpo) / derg
                  endif
               enddo
            enddo
            
         enddo
      enddo


cdebug      if (ndid==0) then
cdebug         write(luout, '("vch1 = ")')
cdebug         call wfat_dprintmat(nbas0, 3, 'e', vch1(:,:,1))
cdebug      endif
      
      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_mo))
     &     call errquit('WFAT: mo deallocation failed', 20, MA_ERR)
      if (.not. MA_free_heap(l_dpo))
     &     call errquit('WFAT: dpmo deallocation failed', 20, MA_ERR)

      end
c     NOTES:
c     1) The center of the dipole moments does not matter for its use in computing
c        the first order perturbation correction because the center only affects the
c        diagonal elements of the dipole matrix in MO basis while for the aforementioned
c        purpose, we only need the non-diagonal elements.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_wfat1_sppos(wfat1, nch, spinid, chpos, 
     &     spinid1, chpos1)

      implicit none
      logical wfat1
      integer nch
      integer spinid(nch)
      double precision chpos(3, nch)
      integer spinid1(1:*)
      double precision chpos1(3, 1:*)

      integer i, j1, j2

      spinid1(1:nch) = spinid
      chpos1(:, 1:nch) = chpos

      if (wfat1) then
         do i = 1, nch
            j1 = nch + 3*(i-1) + 1
            j2 = j1 + 2
            spinid1(j1:j2) = spinid(i)
            chpos1(1,j1:j2) = chpos(1,i)
            chpos1(2,j1:j2) = chpos(2,i)
            chpos1(3,j1:j2) = chpos(3,i)
         enddo
      endif

      end
c=======================================================
