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


c=======================================================
c=======================================================
      subroutine wfat_hint_diag(bt, gm, nsub, maxn, ndgs, dthr, dps, 
     &     isdiag, hi, ev)

      implicit none
      double precision bt
      double precision gm
      integer nsub                                ! [input], the number of degenerate subspaces, a nondegenerate state is considered a degenerate state in a subspace of dimension 1.
      integer maxn                                ! [input], the largest dimension of all degenerate subspaces, i.e. maxn=maxval(ndgs).
      integer ndgs(nsub)                          ! [input], the dimension of each degenerate subspace.
      double precision dthr
      double precision dps(maxn, maxn, 3, nsub)   ! [input], electronic dipole moment matrix in MF in each degenerate subspace. The 3 indices in the 3rd dimension corresponds to x, y, z.
      logical isdiag(nsub)                        ! [output], isdiag(i) is true if the H_int of the i-th subspace is diagonal, otherwise, false.
      double precision hi(maxn**2,nsub)           ! [output], mixing coeffs. of the H_int diagonalizing orbitals.
      double precision ev(maxn,nsub)              ! [output], the corresponding eigenvalues of the orbitals in hi.
c     ^^ On exit, the first ndgs(i)**2 elements of the i-th column of hi contain the flattened form of the mixing coeffs.
      
      integer i, nn, k
      double precision protmat(3,3), cu(3)
      

c     ==== Get a passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      cu = protmat(:,3)
      
      do i = 1, nsub
         nn = ndgs(i)
         isdiag(i) = .true.

c        ==== Calculate static polarizability tensor ====
         if (nn > 1) then
c           ==== Construct H_int within the degenerate MO's ====

            hi(1:nn**2,i) = 0.0d0
            do k = 1, 3
               call daxpy(nn**2, cu(k), dps(1:nn,1:nn,k,i), 1, hi(1,i), 
     &              1)
            enddo
            call wfat_isdiag(nn, hi(1,i), dthr, isdiag(i))

            if (isdiag(i)) then
c              == Identity matrix ==
               hi(:,i) = 0.0d0
               call dcopy(nn, hi(1,i), nn, ev(1,i), 1)    ! The eigvals are simply the diag. elements of hi(1,i).
               call wfat_get_identity(nn, hi(1,i))
            else
c              == Diagonalize H_int in this manifold ==
               call wfat_dsyev_wrap(nn, hi(1,i), ev(1,i)) ! 2)
            endif

         endif
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_diag_update_w0(bt, gm, zi, nch, chmmt, nsub, 
     &     maxn, dgs, ndgs, dthr, dps, cf, isdiag, scr, nb, nptp, psi, 
     &     vfpsi, pos)

      implicit none
#include "stdio.fh"
#include "global.fh"
      double precision bt
      double precision gm
      double precision zi
      integer nch
      double precision chmmt(3,nch)
      integer nsub                                ! [input], the number of degenerate subspaces, a nondegenerate state is considered a degenerate state in a subspace of dimension 1.
      integer maxn                                ! [input], the largest dimension of all degenerate subspaces, i.e. maxn=maxval(ndgs).
      integer dgs(maxn, nsub)                     ! [input], the channel ID arranged into degenerate subspaces (columns). Each subspace (column) must only have states of the same spin.
      integer ndgs(nsub)                          ! [input], the dimension of each degenerate subspace.
      double precision dthr
      double precision dps(maxn, maxn, 3, nsub) ! [input], electronic dipole moment matrix in MF in each degenerate subspace. The 3 indices in the 3rd dimension corresponds to x, y, z.
      double precision cf(maxn**2,nsub)
      logical isdiag(nsub)
      double precision scr(1:*)                   ! [input], scratch array, its minimum size is maxn**2 + maxn + max(nb,nptp)*maxn. 
      integer nb                                  ! [input], number of basis representing psi.
      integer nptp                                ! [input], number of points where vfpsi is evaluated.
      double precision psi(nb,nch)                ! [in/output], some vectors to be mixed.
      double precision vfpsi(nptp,nch)            ! [in/output], vectors resulting from the action of V_eff on the vectors in psi.
      double precision pos(3,nch)                 ! [in/output], channel positions.

      integer p_hi, p_ev, p_lc, i, nn, k, j
      double precision protmat(3,3), cu(3), pos0(3,nch), d_tot(3)
      

      if (ga_nodeid()==0) then
         do j = 1, nch
            write(luout, '("before", i4)') j
            write(luout, '("  chmmt = ", 3es12.5)') chmmt(:,j)
            write(luout, '("  pos = ", 3es12.5)') pos(:,j)
         enddo
      endif
      
c     ==== Double precision scratch pointer ====
c      p_hi = 1
c      p_ev = p_hi + maxn**2
c      p_lc = p_ev + maxn
      p_lc = 1
      
c     ==== Get a passive rotation matrix ===
c      call wfat_passrot(bt, gm, protmat)
c      cu = protmat(:,3)

c      d_tot = zi*pos(:,1) - chmmt(:,1) ! See eq. (B1) in Ref. [1]. Any orbital works, and here, the one with ID=1 is chosen.
      d_tot = -chmmt(:,1)    ! See eq. (B1) in Ref. [1]. Any orbital works, and here, the one with ID=1 is chosen.
      if (ga_nodeid()==0) write(luout, '("neutral D = "3es12.5)') d_tot
      
      do i = 1, nsub
         nn = ndgs(i)
c         isdiag = .true.

c        ==== Calculate static polarizability tensor ====
         if (nn > 1) then
c           ==== Construct H_int within the degenerate MO's ====
c            scr(p_hi:p_hi+nn**2-1) = 0.0d0
c            do k = 1, 3
c               call wfat_matrix_sumabsorb(nn, nn, scr(p_hi),
c     &              dps(1:nn,1:nn,k,i), 1.0d0, cu(k))
c            enddo
c            call wfat_isdiag(nn, scr(p_hi), dthr, isdiag)

            if (.not. isdiag(i)) then
c              == Diagonalize H_int in this manifold ==
c               call wfat_dsyev_wrap(nn, scr(p_hi), scr(p_ev))    ! 2)

c              == Update vectors ==
               call wfat_lincomb_vec(nb, nch, nn, cf(1,i), dgs(:,i),
     &              scr(p_lc), psi)
               call wfat_lincomb_vec(nptp, nch, nn, cf(1,i), dgs(:,i),
     &              scr(p_lc), vfpsi)

c              == Update channel positions ==
               call wfat_chpos_cf(nn, cf(1,i), dps(:nn,:nn,:,i), d_tot,
     &              zi, pos0)
               do j = 1, nn    ! Loop over H_int-diagonalizing orbitals
                  pos(:,dgs(j,i)) = pos0(:,j)
               enddo
            endif
         endif
      enddo

      if (ga_nodeid()==0) then
         do j = 1, nch
            write(luout, '("after", i4)') j
            write(luout, '("  pos = ", 3es12.5)') pos(:,j)
         enddo
      endif
      
      end
c     References:
c     [1] doi.org/10.1103/4xqd-yhdt
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_diag_update_w1(bt, gm, zi, nch, nsub, maxn, dgs, 
     &     ndgs, dthr, dps, pbs, cf, isdiag, scr, nb, nptp, psi, vfpsi,
     &     pb)

      implicit none
#include "global.fh"
#include "stdio.fh"
      double precision bt
      double precision gm
      double precision zi
      integer nch
      integer nsub                                ! [input], the number of degenerate subspaces, a nondegenerate state is considered a degenerate state in a subspace of dimension 1.
      integer maxn                                ! [input], the largest dimension of all degenerate subspaces, i.e. maxn=maxval(ndgs).
      integer dgs(maxn, nsub)                     ! [input], the channel ID arranged into degenerate subspaces (columns). Each subspace (column) must only have states of the same spin.
      integer ndgs(nsub)                          ! [input], the dimension of each degenerate subspace.
      double precision dthr
      double precision dps(maxn, maxn, 3, nsub)   ! [input], electronic dipole moment matrix in MF in each degenerate subspace. The 3 indices in the 3rd dimension corresponds to x, y, z.
      double precision pbs(maxn, maxn, 6, nsub) ! [input], static polarizability matrix in MF in each degenerate subspace. The 6 indices in the 3rd dimension corresponds to xx, yy, zz, xy, xz, yz.
      double precision cf(maxn**2,nsub)
      logical isdiag(nsub)
      double precision scr(1:*)                   ! [input], scratch array, its minimum size is maxn**2 + maxn + max(nb,nptp)*maxn.
      integer nb                                  ! [input], number of bases representing psi.
      integer nptp                                ! [input], number of points where vfpsi is evaluated.
      double precision psi(nb,3*nch)              ! [in/output], some vectors to be mixed, it can be orbitals or potential acting on orbitals.
      double precision vfpsi(nptp,3*nch)          ! [in/output], some vectors to be mixed, it can be orbitals or potential acting on orbitals.
      double precision pb(3,3,nch)                ! [output], MF polarizability tensor for each ionized channels.
      
      integer p_hi, p_ev, p_lc, i_hi, i, nn, k, j, jj
      double precision protmat(3,3), cu(3)


c     ==== Double precision scratch pointer ====
c      p_hi = 1
c      p_ev = p_hi + maxn**2
c      p_lc = p_ev + maxn
      p_lc = 1

c     ==== Get a passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      cu = protmat(:,3)


      do i = 1, nsub
         nn = ndgs(i)
c         isdiag = .true.

c        ==== Calculate static polarizability tensor ====
         if (nn > 1) then
c           ==== Construct H_int within the degenerate MO's ====
c            scr(p_hi:p_hi+nn**2-1) = 0.0d0
c            do k = 1, 3
c               call wfat_matrix_sumabsorb(nn, nn, scr(p_hi),
c     &              dps(1:nn,1:nn,k,i), 1.0d0, cu(k))
c            enddo
c            call wfat_isdiag(nn, scr(p_hi), dthr, isdiag)

c           ==== If the above H_int is not diagonal ... ====
            if (.not. isdiag(i)) then
               if (ga_nodeid()==0)
     &              write(luout, '("not diag. ", 2f8.4)') bt, gm
               
c              == Diagonalize H_int in this manifold ==
c               call wfat_dsyev_wrap(nn, scr(p_hi), scr(p_ev))     ! 2)

c              == Update x, y, and z correction vectors ==

               if (ga_nodeid()==0)
     &              write(luout, '("before sum x ", es13.6)')
     &              sum(psi(:,1::3))
               call wfat_lincomb_vec(nb, nch, nn, cf(1,i), dgs(:,i),
     &              scr(p_lc), psi(:,1::3))
               if (ga_nodeid()==0)
     &              write(luout, '("after sum x ", es13.6)')
     &              sum(psi(:,1::3))
               call wfat_lincomb_vec(nb, nch, nn, cf(1,i), dgs(:,i),
     &              scr(p_lc), psi(:,2::3))
               call wfat_lincomb_vec(nb, nch, nn, cf(1,i), dgs(:,i),
     &              scr(p_lc), psi(:,3::3))
               
c              == Update V_eff acting on the x, y, and z correction vectors ==
               call wfat_lincomb_vec(nptp, nch, nn, cf(1,i), dgs(:,i),
     &              scr(p_lc), vfpsi(:,1::3))
               call wfat_lincomb_vec(nptp, nch, nn, cf(1,i), dgs(:,i),
     &              scr(p_lc), vfpsi(:,2::3))
               call wfat_lincomb_vec(nptp, nch, nn, cf(1,i), dgs(:,i),
     &              scr(p_lc), vfpsi(:,3::3))

c              == Update the polarizabilities using ==
c              == the H_int-diagonalizing channels  ==           ! 1)
               do j = 1, nn    ! Loop over the H_int-diagonalizing channels.
                  i_hi = 1 + (j-1)*nn
                  jj = dgs(j,i)    ! Where in pb the j-th eigenstate should go.
                  call wfat_construct_polb_cf(nn, cf(i_hi,i),
     &                 pbs(1:nn,1:nn,:,i), pb(:,:,jj))
               enddo
c           ==== If the above H_int is diagonal ... ====
            else
               if (ga_nodeid()==0)
     &              write(luout, '("diag. ", 2f8.4)') bt, gm
               
               do j = 1, nn    ! Loop over the H_int-diagonalizing channels.
                  jj = dgs(j,i)      ! Where in pb the j-th eigenstate should go.
                  call wfat_construct_polb(pbs(j,j,:,i), pb(:,:,jj))
               enddo
            endif

c        ==== If non-degenerate ====
         else
            jj = dgs(1,i)      ! Where in pb the j-th eigenstate should go.
            call wfat_construct_polb(pbs(1,1,:,i), pb(:,:,jj))
         endif
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chpos_cf(n, cf, dps, d, zi, pos)
c     DESCRIPTION:
c     Calculates channel positions of orbitals given the expansion coefficients
c     of these orbitals in terms of another set of orbitals within which the
c     dipole matrices are known.
      implicit none
      integer n
      double precision cf(n,n)     ! [input], expansion coefficients in terms of the old orbitals.
      double precision dps(n,n,3)  ! [input], dipole matrices in the rep. of the old orbitals.
      double precision d(3)        ! [input], total dipole of the neutral.
      double precision zi
      double precision pos(3,n)  ! [output], channel positions of the desired/new orbitals.

      integer i, k
      double precision mu_new(3)

      do i = 1, n
         do k = 1, 3
            call wfat_1e_expval_1eopr(n, cf(:,i), dps(:,:,k), mu_new(k))
         enddo
         pos(:,i) = (d + mu_new) / zi
      enddo

      end
c=======================================================
      
      
