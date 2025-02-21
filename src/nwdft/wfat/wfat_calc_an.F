c============================================================
c============================================================
      subroutine wfat_get_srfc1(bt, gm, nat, nch, nqn, lmax, mid, erg, 
     &     chmmt, intl, intl_j, intl_k1, intl_k2, nsub, maxn, dgs, ndgs,
     &     dthr, dps, pbs, nopw, scr, scrx, pb, gn, jn, kn)

      implicit none
#include "stdio.fh"      
#include "errquit.fh"
      double precision bt
      double precision gm
      integer nat
      integer nch
      integer nqn
      integer lmax
      integer mid(nqn)                                      ! [input],
      double precision erg(nch)
      double precision chmmt(3,nch)
      double complex intl((1+lmax)**2, nat, nqn, nch)       ! [input],
      double complex intl_j((1+lmax)**2, nat, nqn, 3, nch)  ! [input],
      double complex intl_k1((1+lmax)**2, nat, nqn, nch)    ! [input],
      double complex intl_k2((1+lmax)**2, nat, nqn, nch)    ! [input],
      integer nsub                                ! [input], the number of degenerate subspaces, a nondegenerate state is considered a degenerate state in a subspace of dimension 1.
      integer maxn                                ! [input], the largest dimension of all degenerate subspaces, i.e. maxn=maxval(ndgs).
      integer dgs(maxn, nsub)                     ! [input], the channel ID arranged into degenerate subspaces (columns). Each subspace (column) must only have states of the same spin.
      integer ndgs(nsub)                          ! [input], the dimension of each degenerate subspace.
      double precision dthr
      double precision dps(maxn, maxn, 3, nsub)   ! [input], The electronic dipole moment matrix in each degenerate subspace. The 3 indices in the 3rd dimension corresponds to x, y, z.
      double precision pbs(maxn, maxn, 6, nsub)   ! [input], The static polarizability matrix in each degenerate subspace. The 6 indices in the 3rd dimension corresponds to xx, yy, zz, xy, xz, yz.
      logical nopw                                ! [input],
      double precision scr(1:*)                   ! [output], scratch array, its minimum size is maxn**2 + maxn.
      double complex scrx(1:*)                    ! [output], scratch array, its minimum size is nat*nqn*nch.
      double precision pb(3,3,nch)                ! [output], the molecular frame polarizability.
      double complex gn(nat, nqn, nch)            ! [output], G_n = \exp{-\kappa*\mu_z} * g_n
      double complex jn(nat, nqn, nch)            ! [output], J_n = \exp{-\kappa*\mu_z} * j_n
      double complex kn(nat, nqn, nch)            ! [output], K_n = \exp{-\kappa*\mu_z} * k_n

      logical isdiag
      integer p_hi, p_ev, p_lc, i, nn, k, j, i_hi, jj, iu
      double precision protmat(3,3), cu(3)


      if (nopw .and. lmax /= 0)
     &     call errquit('WFAT: If nopw is .true., then lmax must be '//
     &     'equal to zero (wfat_get_srfc1).', 43, UERR)

      
c     ==== Double precision scratch pointer ====
      p_hi = 1
      p_ev = p_hi + maxn**2

c     ==== Double complex scratch pointer ====
      p_lc = 1


c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      cu = protmat(:,3)
      
      
      do i = 1, nsub
         nn = ndgs(i)
         isdiag = .true.

c        ==== Calculate static polarizability tensor ====
         if (nn > 1) then
c           ==== Construct H_int within the degenerate MO's ====
            scr(p_hi:p_hi+nn**2-1) = 0.0d0
            do k = 1, 3
               call wfat_matrix_sumabsorb(nn, nn, scr(p_hi),
     &              dps(1:nn,1:nn,k,i), 1.0d0, cu(k))
            enddo
            call wfat_isdiag(nn, scr(p_hi), dthr, isdiag)

c           ==== If the above H_int is not diagonal ... ====
            if (.not. isdiag) then
c              == Diagonalize H_int in this manifold ==
               call wfat_dsyev_wrap(nn, scr(p_hi), scr(p_ev))

c              == Update the polarizabilities using ==
c              == the H_int-diagonalizing channels  ==
               do j = 1, nn    ! Loop over the H_int-diagonalizing channels.
                  i_hi = p_hi + (j-1)*nn
                  jj = dgs(j,i)      ! Where in pb the j-th eigenstate should go.
                  call wfat_construct_polb_cf(nn, scr(i_hi), 
     &                 pbs(1:nn,1:nn,:,i), pb(:,:,jj))             ! 1)
                  write(luout, '("nondiag: ", 2f10.5, i6, es21.12)',
     &                 advance='no') bt, gm, j, scr(p_ev+j-1)
                  do iu = 1, nn
                     write(luout, '(" u", es15.6)', advance='no')
     &                    scr(i_hi+iu-1)
                  enddo
                  write(luout, '("")')
               enddo
c           ==== If the above H_int is diagonal ... ====
            else
               do j = 1, nn    ! Loop over the H_int-diagonalizing channels.
                  jj = dgs(j,i)      ! Where in pb the j-th eigenstate should go.
                  call wfat_construct_polb(pbs(j,j,:,i), pb(:,:,jj))
               enddo
            endif

c        ==== If non-degenerate ====
         else
            jj = dgs(1,i)      ! Where in pb the j-th eigenstate should go.
            call wfat_construct_polb(pbs(1,1,:,i), pb(:,:,jj))
c            call wfat_dprintmat(1, 6, 'f', pbs(1,1,:,i))
c            call wfat_dprintmat(3, 3, 'f', pb(:,:,jj))
         endif

         
c        ==== Calculate g_\nu and a_\nu of the non-diagonalizing ====
c        ====   channels within the current degenerate subspace  ====
         do j = 1, nn
            jj = dgs(j,i)
            call wfat_calc_gn(nat, nqn, mid, lmax, bt, gm, 
     &           intl(:,:,:,jj), nopw, gn(:,:,jj))
            call wfat_calc_an(nat, nqn, mid, lmax, bt, gm, chmmt(:,jj), 
     &           intl_j(:,:,:,:,jj), intl_k1(:,:,:,jj), 
     &           intl_k2(:,:,:,jj), nopw, jn(:,:,jj), kn(:,:,jj))
         enddo

         
c        ==== Calculate g_\nu and a_\nu of the diagonalizing  ====
c        ==== channels within the current degenerate subspace ====
         if (nn > 1) then
            if (.not. isdiag) then
               call wfat_lincomb_gn(nat, nch, nqn, nn, scr(p_hi), 
     &              dgs(:,i), scrx(p_lc), gn)
               call wfat_lincomb_gn(nat, nch, nqn, nn, scr(p_hi), 
     &              dgs(:,i), scrx(p_lc), jn)
               call wfat_lincomb_gn(nat, nch, nqn, nn, scr(p_hi), 
     &              dgs(:,i), scrx(p_lc), kn)
            endif
         endif
         
      enddo


c     ==== Calculate structure factors ====
      do i = 1, nch
         call wfat_calc_srfc0(nat, nqn, bt, gm, erg(i), chmmt(:,i), 
     &        gn(:,:,i))
         call wfat_calc_srfc1(nat, nqn, bt, gm, erg(i), chmmt(:,i),
     &        jn(:,:,i), kn(:,:,i))
      enddo
      
      end
c     NOTE:
c     1) After computing the new degenerate channels above, only the 
c        polarizabilities are updated. The dipole moments are not because
c        it is assumed that the original degenerate channels have the same
c        dipole moment vector. This assumption is needed to avoid
c        recomputation of the WFAT(1) integrals as these integrals depend
c        on the dipole moment vector of the ionizing channel through the
c        choice of origin. Unfortunately, this means that degenerate channels
c        that do not have the same dipole moment vectors cannot be handled
c        by this code. According to my experience, however, I have never
c        encountered degenerate eigenstates of a molecular Hamiltonian that
c        do not share the same dipole moment vector.
c============================================================


c============================================================
c============================================================
      subroutine wfat_calc_an(nat, nqn, mid, lmax, bt, gm, mmt,
     &     intl_j, intl_k1, intl_k2, nopw, jn, kn)
c     DESCRIPTION:
c     Performs the summation of the partial wave expansion of the WFAT(1) integrals.
      implicit none
#include "errquit.fh"
      integer nat
      integer nqn
      integer mid(nqn)
      integer lmax
      double precision bt
      double precision gm
      double precision mmt(3)
      double complex intl_j((1+lmax)**2, nat, nqn, 3)  ! [input],
      double complex intl_k1((1+lmax)**2, nat, nqn)    ! [input],
      double complex intl_k2((1+lmax)**2, nat, nqn)    ! [input],
      logical nopw
      double complex jn(nat, nqn)              ! [output],
      double complex kn(nat, nqn)              ! [output],

      integer iqn, ic, l, m, ish, i, l0
      double precision protmat(3,3), cu(3), btr, gmr, mu(3), deg2rad
      double precision wgf
      double complex aa, bb, czero, imnum
      parameter ( deg2rad = 0.01745329251994d0 )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )


      if (nopw .and. lmax /= 0)
     &     call errquit('WFAT: If nopw is .true., then lmax must be '//
     &     'equal to zero (wfat_calc_an).', 43, UERR)

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      cu = protmat(:,3)      
      btr = bt * deg2rad
      gmr = gm * deg2rad
      

c     IMPORTANT:
c     At this point, if ich corresponds to a degenerate state, mmt (dipole moment)
c     must be associated with the interaction-Hamiltonian-diagonalizing eigenstates,
c     which are a linear combination of the original degenerate states. This requirement,
c     however, is always satisfied as long as the original degenerate states have the
c     same dipole moment vector. This is the assumption which must be satisfied,
c     otherwise the ionization rate cannot be calculated by this code.
      mu = matmul(transpose(protmat), -mmt)
      do iqn = 1, nqn

         do ic = 1, nat
            jn(ic, iqn) = czero
            kn(ic, iqn) = czero
            if (nopw) then
               l0 = 0
            else
               l0 = abs(mid(iqn))
            endif
            
            do l = l0, lmax
               do m = -l, l
                  ish = l**2 + l + m + 1
                  
                  aa = intl_k1(ish,ic,iqn) +
     &                 mu(3)*intl_k2(ish,ic,iqn)

                  bb = czero
                  do i = 1, 3
                     bb = bb + cu(i)*intl_j(ish,ic,iqn,i)
                  enddo

                  if (nopw) then
                     wgf = 1.0d0
                  else
                     call wfat_wignerf(l, mid(iqn), m, btr, wgf)
                  endif
                  
                  jn(ic, iqn) = jn(ic, iqn) +
     &                          bb * wgf * exp(-imnum * m * gmr)
                  kn(ic, iqn) = kn(ic, iqn) +
     &                          aa * wgf * exp(-imnum * m * gmr)
               enddo
            enddo
            
         enddo
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_lincomb_gn(nat, nch, nqn, n, cf, cid, scr, s)

      implicit none
      integer nat
      integer nch
      integer nqn
      integer n
      double precision cf(n,n)
      integer cid(n)                   ! [input], 1 <= cid(i) <= nch
      double complex scr(nat, nqn, n)  ! [output],
      double complex s(nat, nqn, nch)  ! [in/output],

      integer j, j_, i, i_
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )

      
c     ==== Construct the linear combinations ====
      do j = 1, n
         scr(:,:,j) = czero
         do i = 1, n
            i_ = cid(i)
            scr(:,:,j) = scr(:,:,j) + cf(i,j)*s(:,:,i_)
         enddo
      enddo

c     ==== Replace the old structure factors with the new ones ====
      do j = 1, n
         j_ = cid(j)
         s(:,:,j_) = scr(:,:,j)
      enddo

      end
c============================================================
