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
      double precision dps(maxn, maxn, 3, nsub)   ! [input], The electronic dipole moment matrix in MF in each degenerate subspace. The 3 indices in the 3rd dimension corresponds to x, y, z.
      double precision pbs(maxn, maxn, 6, nsub)   ! [input], The static polarizability matrix in MF in each degenerate subspace. The 6 indices in the 3rd dimension corresponds to xx, yy, zz, xy, xz, yz.
      logical nopw                                ! [input],
      double precision scr(1:*)                   ! [output], scratch array, its minimum size is maxn**2 + maxn.
      double complex scrx(1:*)                    ! [output], scratch array, its minimum size is nat*nqn*nch.
      double precision pb(3,3,nch)                ! [output], MF polarizability tensor for each ionized channels.
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
               call wfat_dsyev_wrap(nn, scr(p_hi), scr(p_ev)) ! 2)
cBAD               call wfat_setorder_dip(bt, gm, nn, dps(1:nn,1:nn,:,i),
cBAD     &              scr(p_ev), scr(p_hi))
cBAD               call wfat_setorder_polb2(bt, gm, nn, pbs(1:nn,1:nn,:,i),
cBAD     &              scr(p_ev), scr(p_hi))

c              == Update the polarizabilities using ==
c              == the H_int-diagonalizing channels  ==          ! 1)
               do j = 1, nn    ! Loop over the H_int-diagonalizing channels.
                  i_hi = p_hi + (j-1)*nn
                  jj = dgs(j,i)      ! Where in pb the j-th eigenstate should go.
                  call wfat_construct_polb_cf(nn, scr(i_hi), 
     &                 pbs(1:nn,1:nn,:,i), pb(:,:,jj))
c                  write(luout, '("nondiag: ", 2f10.5, i6, 3es21.12)',
c     &                 advance='no') bt, gm, j, scr(p_ev+j-1), dum(j),
c     &                 dum(3)
               enddo
c           ==== If the above H_int is diagonal ... ====
            else
               do j = 1, nn    ! Loop over the H_int-diagonalizing channels.
                  jj = dgs(j,i)      ! Where in pb the j-th eigenstate should go.
                  call wfat_construct_polb(pbs(j,j,:,i), pb(:,:,jj))
               enddo
            endif
cdebug            call wfat_print_degdip2(bt, gm, isdiag, nn, 
cdebug     &           dps(1:nn,1:nn,:,i), scr(p_hi))

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
c        the shifted dipole moment in WFAT is equal to the total dipole
c        vector of the neutral, hence independent of unitary transform.
c        The part where the dipole vector is affected by unitary
c        transform is the choice of origin. In this code, the new origin
c        caused by the unitary transform is assumed very close to the old
c        one. This requires that the off-diagonal elements of the
c        interaction matrix are much smaller than the diagonal elements.
c        Because of this assumption, the dipole moment is not updated.
c     2) I turns out that the eigenvalues of the interaction operator
c        provide a consistent way of ordering the diagonalizing degenerate
c        orbitals. Therefore, no additional ordering is needed beyond
c        that imposed by the diagonalizing subroutine which automatically
c        orders the these new orbitals based on their associated
c        eigenvalues. A consistent ordering can be checked by looking
c        at the 2D plot of the absolute value of structure factors, and
c        checking that no abrupt changes are present.
c============================================================


c============================================================
c============================================================
      subroutine wfat_calc_an(nat, nqn, mid, lmax, bt, gm, mmt,
     &     intl_j, intl_k1, intl_k2, nopw, jn, kn)
c     DESCRIPTION:
c     Performs the summation of the partial wave expansion of the WFAT(1) integrals.
      implicit none
#include "stdio.fh"
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
c     which are a linear combination of the original degenerate states.
      mu = matmul(transpose(protmat), -mmt)

cdebug      write(luout, '("check mu ", 2f12.6, es18.10)') bt, gm, mu(3)    ! 1)
c     1) This line is used to show that mu(3) is the same for all channels
c        since it is rotated from the total dipole moment of the neutral
c        molecule. This is further a consequence of requiring the origin
c        to be so that the total dipole of the cation is zero.
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


c============================================================
c============================================================
      subroutine wfat_print_degdip1(bt, gm, isdiag, ndeg, dip_old, u)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      double precision bt
      double precision gm
      logical isdiag
      integer ndeg
      double precision dip_old(ndeg,ndeg,3)    ! [input], MF dipole matrix in the degenerate subspace.
      double precision u(ndeg,ndeg)            ! [input], coefficients of the dipole-diagonalizing channels in terms of the dipole-non-diagonalizing channels.

      character*7 s
      integer i, k, iu, i_dn, k_dn, l_dn
      double precision dst, dvec(3), t

cdebug      double precision mu_old(3), mu_new(3), dip_new(3), protmat(3,3)
      

      if (.not.ma_push_get(mt_dbl, ndeg*ndeg*3, 'dnew', l_dn, k_dn))
     &     call errquit('WFAT: cannot allocate dnew', 56, MA_ERR)
      do k = 1, 3
         i_dn = k_dn + (k-1)*ndeg**2
         call wfat_represent_sym(ndeg, ndeg, u, dip_old(:,:,k),
     &        dbl_mb(i_dn))
      enddo

      if (isdiag) then
         s = 'diag'
      else
         s = 'ndiag'
      endif
      write(luout, '(a7, 2f10.5)', advance='no') s, bt, gm
      do i = 1, ndeg
         write(luout, '("        ")', advance='no')
         if (isdiag) then
            write(luout, '(es19.10)', advance='no') 0.0d0
            do iu = 1, ndeg
               t = 0.0d0
               if (iu == i) t = 1.0d0
               write(luout, '(es15.6)', advance='no') t
            enddo
         else
            do k = 1, 3
               i_dn = k_dn + (k-1)*ndeg**2 + (i-1)*ndeg + (i-1)
               dvec(k) = dip_old(i,i,k) - dbl_mb(i_dn)
            enddo
            dst = sqrt(sum(dvec**2))
            write(luout, '(es19.10)', advance='no') dst
            
cdebug            call wfat_passrot(bt, gm, protmat)
cdebug            mu_old = matmul(transpose(protmat), dip_old(i,i,:))
cdebug            do k = 1, 3
cdebug               i_dn = k_dn + (k-1)*ndeg**2 + (i-1)*ndeg + (i-1)
cdebug               dip_new(k) = dbl_mb(i_dn)
cdebug            enddo
cdebug            mu_new = matmul(transpose(protmat), dip_new)
cdebug            write(luout, '(2es19.10)', advance='no')
cdebug     &           mu_old(3), mu_new(3)
            
            do iu = 1, ndeg
               write(luout, '(es15.6)', advance='no') u(iu,i)
            enddo
         endif
      enddo
      write(luout, '("")')

      if (.not. MA_pop_stack(l_dn))
     &     call errquit('WFAT: cannot deallocate dnew', 56, MA_ERR)

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_print_degdip2(bt, gm, isdiag, ndeg, dip_old, u)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      double precision bt
      double precision gm
      logical isdiag
      integer ndeg
      double precision dip_old(ndeg,ndeg,3)    ! [input], MF dipole matrix in the degenerate subspace.
      double precision u(ndeg,ndeg)            ! [input], coefficients of the dipole-diagonalizing channels in terms of the dipole-non-diagonalizing channels.

      character*7 s
      integer i, iu, l_md, k_md
      double precision t


      if (.not.ma_push_get(mt_dbl, ndeg, 'mindst', l_md, k_md))
     &     call errquit('WFAT: cannot allocate mindst', 39, MA_ERR)

      
      if (isdiag) then
         s = 'diag'
      else
         s = 'ndiag'
         call wfat_get_mindst(ndeg, dip_old, u, dbl_mb(k_md))
      endif
      write(luout, '(a7, 2f10.5)', advance='no') s, bt, gm

      
      do i = 1, ndeg
         write(luout, '("        ")', advance='no')
         if (isdiag) then
            write(luout, '(es19.10)', advance='no') 0.0d0
            do iu = 1, ndeg
               t = 0.0d0
               if (iu == i) t = 1.0d0
               write(luout, '(es15.6)', advance='no') t
            enddo
         else
            write(luout, '(es19.10)', advance='no') dbl_mb(k_md+i-1)
            do iu = 1, ndeg
               write(luout, '(es15.6)', advance='no') u(iu,i)
            enddo
         endif
      enddo
      write(luout, '("")')

      
      if (.not. MA_pop_stack(l_md))
     &     call errquit('WFAT: cannot deallocate dnew', 56, MA_ERR)

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_get_mindst(ndeg, dip_old, u, mindst)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
c      double precision bt
c      double precision gm
      integer ndeg
      double precision dip_old(ndeg,ndeg,3)    ! [input], MF dipole matrix in the degenerate subspace.
      double precision u(ndeg,ndeg)            ! [input], coefficients of the dipole-diagonalizing channels in terms of the dipole-non-diagonalizing channels.
      double precision mindst(ndeg)
      
      integer l_vo, k_vo, l_vn, k_vn, l_d, k_d, l_p, k_p
      integer i, i_vo, k, i_vn, j, i_vn1, i_vn2, i_vo1, i_vo2, i_p
      integer mincol, minrow
      double precision dvec(3)

      
      if (.not.ma_push_get(mt_dbl, 3*ndeg, 'veco', l_vo, k_vo))
     &     call errquit('WFAT: cannot allocate veco', 56, MA_ERR)
      if (.not.ma_push_get(mt_dbl, 3*ndeg, 'vecn', l_vn, k_vn))
     &     call errquit('WFAT: cannot allocate vecn', 56, MA_ERR)
      
      do i = 1, ndeg
         i_vo = k_vo + (i-1)*3
         call dcopy(3, dip_old(i,i,1), ndeg**2, dbl_mb(i_vo), 1)
         do k = 1, 3
            i_vn = k_vn + (i-1)*3 + (k-1)
            call wfat_1e_expval_1eopr(ndeg, u(:,i), dip_old(:,:,k),
     &           dbl_mb(i_vn))
         enddo
      enddo
cdebug      call wfat_dprintmat(3, ndeg, 'e', dbl_mb(k_vo))
cdebug      call wfat_dprintmat(3, ndeg, 'e', dbl_mb(k_vn))
      

      if (.not.ma_push_get(mt_dbl, ndeg**2, 'dist', l_d, k_d))
     &     call errquit('WFAT: cannot allocate dist', 56, MA_ERR)

      do j = 1, ndeg   ! loop over columns (new orbs)
         i_vn1 = k_vn + (j-1)*3
         i_vn2 = i_vn1 + 2
         do i = 1, ndeg   ! loop over rows (old orbs)
            i_vo1 = k_vo + (i-1)*3
            i_vo2 = i_vo1 + 2
                     
            dvec = dbl_mb(i_vn1:i_vn2) - dbl_mb(i_vo1:i_vo2)
            dbl_mb(k_d+(j-1)*ndeg+(i-1)) = sqrt(sum(dvec**2))
         enddo
      enddo


      if (.not.ma_push_get(mt_int, 2*ndeg, 'pairs', l_p, k_p))
     &     call errquit('WFAT: cannot allocate pairs', 12, MA_ERR)
      call wfat_get_minpairs(ndeg, dbl_mb(k_d), int_mb(k_p))

      do i = 1, ndeg
         i_p = k_p + (i-1)*2
         mincol = int_mb(i_p)
         minrow = int_mb(i_p+1)
         if (i /= mincol) call errquit('WFAT: A bug in wfat_get_mindst',
     &        67, UERR)
         mindst(i) = dbl_mb(k_d+(mincol-1)*ndeg+(minrow-1))
      enddo

      if (.not. MA_pop_stack(l_p))
     &     call errquit('WFAT: cannot deallocate pairs', 56, MA_ERR)
      if (.not. MA_pop_stack(l_d))
     &     call errquit('WFAT: cannot deallocate dist', 56, MA_ERR)
      if (.not. MA_pop_stack(l_vn))
     &     call errquit('WFAT: cannot deallocate vecn', 56, MA_ERR)
      if (.not. MA_pop_stack(l_vo))
     &     call errquit('WFAT: cannot deallocate veco', 56, MA_ERR)
      
      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_get_minpairs(ndeg, dst, pair)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
      integer ndeg
      double precision dst(ndeg,ndeg)
      integer pair(2,ndeg)

      logical coltodo(ndeg), rowtodo(ndeg)
      integer k, j, i, mincol, minrow
      double precision minnow, minprev, maxvalue


cdebug      call wfat_dprintmat(ndeg, ndeg, 'e', dst)
      maxvalue = maxval(dst) + 10
            
      coltodo(:) = .true.
      rowtodo(:) = .true.
      do k = 1, ndeg

         minprev = maxvalue
         do j = 1, ndeg
            if (coltodo(j)) then
               minnow = minval(dst(:,j))
               if (minnow <= minprev) then
                  minprev = minnow
                  mincol = j
               endif
            endif
         enddo
cdebug         write(luout, '("mc: ", 2i4, es14.6)') k, mincol, minnow
         if (coltodo(mincol)) then
            coltodo(mincol) = .false.
         else
            call errquit('WFAT: Repeated column.', 28, UERR)
         endif

         minprev = maxvalue
         do i = 1, ndeg
            if (rowtodo(i)) then
               minnow = dst(i,mincol)
               if (minnow <= minprev) then
                  minprev = minnow
                  minrow = i
               endif
            endif
         enddo
         if (rowtodo(minrow)) then
            rowtodo(minrow) = .false.
         else
            call errquit('WFAT: Repeated row.', 28, UERR)
         endif

         pair(:,mincol) = (/ mincol, minrow /)
      enddo

      end
c============================================================
      
