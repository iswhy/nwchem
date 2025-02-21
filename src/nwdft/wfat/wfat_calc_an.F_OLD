c============================================================
c============================================================
      subroutine wfat_get_srfc1(bt, gm, nat, nb, ns, nch, nqn, 
     &     lmax, orbid, spinid, dgid0, dgn0, mid, dpmo, erg, chmmt, pb,
     &     chpos, mthr, dthr, intl, intl_j, intl_k1, intl_k2, scr, scrx,
     &     chmmtd, pbd, gn, jn, kn)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      double precision bt       ! [input], beta in degrees
      double precision gm     ! [input], gamma in degrees
      integer nat
      integer nb
      integer ns
      integer nch
      integer nqn
      integer lmax
      integer orbid(nch)
      integer spinid(nch)
      integer dgid0(nch)                                    ! [input],
      integer dgn0(nch)                                     ! [input],
      integer mid(nqn)                                      ! [input],
      double precision dpmo(nb,nb,ns,3)                     ! [input],
      double precision erg(nch)                             ! [input],
      double precision chmmt(3,nch)                         ! [input],
      double precision pb(3,3,nch)                          ! [input],
      double precision chpos(3,nch)                         ! [input],
      double precision mthr                                 ! [input], UNUSED. Threshold for determining if the dipole moments in the degenerate channels are identical.
      double precision dthr                                 ! [input], threhsold for determining if the interaction Hamiltonian in the degenerate subspace is diagonal.
      double complex intl((1+lmax)**2, nat, nqn, nch)       ! [input],
      double complex intl_j((1+lmax)**2, nat, nqn, 3, nch)  ! [input],
      double complex intl_k1((1+lmax)**2, nat, nqn, nch)    ! [input],
      double complex intl_k2((1+lmax)**2, nat, nqn, nch)    ! [input],
      double precision scr(1:*)                             ! [output], minimum size is 2*nb**2 + nb
      double complex scrx(1:*)                              ! [output], minimum size is nat*nqn*nch
      double precision chmmtd(3,nch)                        ! [output],
      double precision pbd(3,3,nch)                         ! [output],
      double complex gn(nat, nqn, nch)                      ! [output], G_n = \exp{-\kappa*\mu_z} * g_n
      double complex jn(nat, nqn, nch)                      ! [output], J_n = \exp{-\kappa*\mu_z} * j_n
      double complex kn(nat, nqn, nch)                      ! [output], K_n = \exp{-\kappa*\mu_z} * k_n

      logical iscalc(nch), isdiag, dgtest, ismeq
      integer nht, ndp, nev, nlc, p_ht0, p_htn, p_dp0, p_dpn, p_ev0
      integer p_evn, p_lc0, p_lcn, i, is, nn, k, j, j_, mfid(nch)
      integer orbid0(nch), ii, jj, ii_
      double precision protmat(3,3), cu(3)

cOLD      parameter(dthr=1.0d-12)

      
c     ==== Double precision scratch ====
      nht = nb**2
      ndp = nb**2
      nev = nb
      
      p_ht0 = 1                 ! p_ht0 = pointer to the (dynamic) submatrix of H_int
      p_htn = p_ht0 + nht - 1
      p_dp0 = p_htn + 1         ! p_dp0 = pointer to the (dynamic) submatrix of dpmo
      p_dpn = p_dp0 + ndp - 1
      p_ev0 = p_dpn + 1         ! p_ev0 = pointer to the (dynamic) eigenvalue of the submatrix of H_int
      p_evn = p_ev0 + nev - 1

c     Double complex scratch ====
      nlc = nat * nqn * nch
      
      p_lc0 = 1
      p_lcn = p_lc0 + nlc - 1

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      cu = protmat(:,3)
      

c     Diagonalize the dipole interaction Hamiltonian within any
c     degenerate manifolds that are ionized.
      chmmtd = chmmt
      pbd = pb
            

      iscalc = .false.
      mfid = 0      
      do i = 1, nch
         is = spinid(i)
         nn = dgn0(i)

         if (.not. iscalc(i)) then
c           ==== Get the degenerate manifold ID ====
            call wfat_get_degen_mfold(nch, orbid, spinid, dgid0, dgn0,
     &           i, mfid, orbid0)
            
            dgtest = (nn > 1)
            isdiag = .true.
            if (dgtest) then
c              ==== Construct H_int within the degenerate MO's ====
               scr(p_ht0:p_ht0+nn**2-1) = 0.0d0
               do k = 1, 3
                  call wfat_submatrix(nb, nb, dpmo(:,:,is,k), nn, 
     &                 nn, orbid0, orbid0, scr(p_dp0))
                  scr(p_ht0:p_ht0+nn**2-1) = scr(p_ht0:p_ht0+nn**2-1) +
     &                                    cu(k)*scr(p_dp0:p_dp0+nn**2-1)
               enddo
               call wfat_isdiag(nn, scr(p_ht0), dthr, isdiag)
               
               if (.not. isdiag) then
cOLD                  call wfat_equal_mmt(nch, nn, mfid, chmmtd, chpos,
cOLD     &                 mthr, ismeq)
cOLD                  if (.not. ismeq) then
cOLD                     do ii = 1, nn
cOLD                        jj = mfid(ii)
cOLD                        write(luout, '("Moment of channel ", i4)') jj
cOLD                        call wfat_dprintmat(1, 3, 'e', chmmtd(:,jj))
cOLD                     enddo
cOLD                     call errquit('WFAT: For the moment, degenerate '//
cOLD     &                    'ionized channels must also have identical '//
cOLD     &                    'dipole moments', 48, UERR)
cOLD                  endif
                  
                  
c                 ==== Diagonalize H_int in this manifold ====
                  call wfat_dsyev_wrap(nn, scr(p_ht0), scr(p_ev0))


c                 ==== Update first moments of ionized orbitals ====
                  call wfat_setorder_polb(nch, nn, mfid, orbid0, is,
     &                 nb, ns, dpmo, scr(p_ht0), pbd)
                  do ii = 1, nn
                     ii_ = mfid(ii)
                     write(luout, '("nondiag: ", 2f10.5, i6, 4es21.12)')
     &                    bt, gm, ii, scr(p_ev0+ii-1), chmmtd(:,ii_)
                  enddo
               endif
            endif

c           ==== Calculate g_\nu and a_\nu ====
            do j = 1, nn
               j_ = mfid(j)
               call wfat_calc_gn(nat, nqn, mid, lmax, bt, gm, 
     &              intl(:,:,:,j_), gn(:,:,j_))
               call wfat_calc_an(nat, nqn, mid, lmax, bt, gm, 
     &              chmmtd(:,j_), intl_j(:,:,:,:,j_),
     &              intl_k1(:,:,:,j_), intl_k2(:,:,:,j_), jn(:,:,j_),
     &              kn(:,:,j_))
            enddo

c           ==== Update g_\nu and a_\nu ====
            if (dgtest) then
               if (.not. isdiag) then
                  call wfat_lincomb_gn(nat, nch, nqn, nn, scr(p_ht0), 
     &                 mfid, scrx(p_lc0), gn)
                  call wfat_lincomb_gn(nat, nch, nqn, nn, scr(p_ht0), 
     &                 mfid, scrx(p_lc0), jn)
                  call wfat_lincomb_gn(nat, nch, nqn, nn, scr(p_ht0), 
     &                 mfid, scrx(p_lc0), kn)
               endif
            endif
         
            do j = 1, nn
               j_ = mfid(j)
               iscalc(j_) = .true.
            enddo
         endif
         
      enddo


c     ==== Calculate structure factors ====
      do i = 1, nch
         call wfat_calc_srfc0(nat, nqn, bt, gm, erg(i), chmmtd(:,i), 
     &        gn(:,:,i))
         call wfat_calc_srfc1(nat, nqn, bt, gm, erg(i), chmmtd(:,i),
     &        jn(:,:,i), kn(:,:,i))
      enddo
      
      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_calc_an(nat, nqn, mid, lmax, bt, gm, mmt,
     &     intl_j, intl_k1, intl_k2, jn, kn)

      implicit none
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
      double complex jn(nat, nqn)              ! [output],
      double complex kn(nat, nqn)              ! [output],

      integer iqn, ic, l, m, ish, i
      double precision protmat(3,3), cu(3), btr, gmr, mu(3), deg2rad
      double precision wgf
      double complex aa, bb, czero, imnum
      parameter ( deg2rad = 0.01745329251994d0 )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      cu = protmat(:,3)      
      btr = bt * deg2rad
      gmr = gm * deg2rad
      

c     IMPORTANT:
c     At this point, if ich corresponds to a degenerate state, mu and alpha must
c     be associated with the interaction-Hamiltonian-diagonalizing eigenstates,
c     which are a linear combination of the original degenerate states.
      mu = matmul(transpose(protmat), -mmt)
      do iqn = 1, nqn
         do ic = 1, nat

            jn(ic, iqn) = czero
            kn(ic, iqn) = czero
            do l = abs(mid(iqn)), lmax
               do m = -l, l
                  ish = l**2 + l + m + 1
                  
                  aa = intl_k1(ish,ic,iqn) +
     &                 mu(3)*intl_k2(ish,ic,iqn)

                  bb = czero
                  do i = 1, 3
                     bb = bb + cu(i)*intl_j(ish,ic,iqn,i)
                  enddo

                  call wfat_wignerf(l, mid(iqn), m, btr, wgf)
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
      subroutine wfat_get_degen_mfold(nch, orbid, spinid, dgid0, dgn0, 
     &     ich, mfid, orbid0)
c     Given ich, :
c     1) mfid will contain the ID's of the ionized channels that are
c        in the same degenerate manifold and has the same spin as ich (it
c        includes ich too) , and
c     2) orbid0 will contain the ID's of the MO's that are in the same
c        degenerate manifold and has the same spin as the MO ich corresponds 
c        to (it includes the MO corresponding to ich too).
c     Hence, for a given ich, mfid and orbid0 constitute two equivalent
c     representations of the same set of degenerate MO's.

      implicit none
#include "errquit.fh"
      integer nch          ! [input],
      integer orbid(nch)   ! [input],
      integer spinid(nch)  ! [input],
      integer dgid0(nch)   ! [input], this array is produced by wfat_detect_degen0.
      integer dgn0(nch)    ! [input], this array is produced by wfat_detect_degen0.
      integer ich          ! [input], an integer whose value satisfies 1 <= ich <= nch.
      integer mfid(1:*)    ! [output], 1 <= mfid(i) <= nch, with i runs from 1 to dgn0(ich).
      integer orbid0(1:*)  ! [output], 1 <= orbid0(i) <= nbas0, with i runs from 1 to dgn0(ich).
      
      integer k, jch, j, j_, jo
      
      
      if (ich > nch) call errquit('WFAT: The argument ich of '//
     &     'wfat_get_degen_mfold must not be larger than nch', 34, UERR)
      
      k = 0
      do jch = 1, nch
         if (dgid0(jch) == dgid0(ich) .and. spinid(jch) == spinid(ich))
     &        then
            k = k + 1
            mfid(k) = jch
         endif
      enddo

      if (k /= dgn0(ich)) call errquit('WFAT: An inconsistency is '//
     &     'detected when identifying degenerate manifold. This is '//
     &     'most likely caused by not all degenerate MO''s within '//
     &     'a given degenerate manifold being included in the '//
     &     'ionized channels.', 23, UERR)

      do j = 1, dgn0(ich)
         j_ = mfid(j)
         jo = orbid(j_)
         orbid0(j) = jo
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
      subroutine wfat_equal_mmt(nch, n, mfid, mmd, ctr, mthr, t)

      implicit none
      integer nch
      integer n
      integer mfid(n)
      double precision mmd(3,nch)
      double precision ctr(3,nch)
      double precision mthr
      logical t

      integer j, j_, i, i_
      double precision d(3), e


      t = .true.

      do j = 1, n
         j_ = mfid(j)
         do i = 1, n
            i_ = mfid(i)

            if (j > i) then
               d = (mmd(:,j_)+ctr(:,j_)) - (mmd(:,i_)+ctr(:,i_))
               e = sqrt( sum(d**2) )
               if (e > mthr) then
                  t = .false.
                  return
               endif
            endif
            
         enddo
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_setorder_polb(nch, n, mfid, orbid0, ispin, nb, ns,
     &     dpmo, cf, pbd)
c     DESCRIPTION:
c     This subroutine reorders the diagonalizing degenerate channels such that
c     they are ordered in the inreasing sum of the diagonal elements of their 
c     respective polarizability tensors. This is done because the (coefficients 
c     of the) diagonalizing degenerate channels change with the orientation angles.
c     Without some sort of standard to fix the ordering of the columns of cf, the
c     angular dependence of, e.g., structure factors and ionization rates will
c     screwed up by ambiguities.
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer nch                             ! [input]
      integer n                               ! [input]
      integer mfid(n)                         ! [input]
      integer orbid0(n)                       ! [input]
      integer ispin                           ! [input]
      integer nb                              ! [input]
      integer ns                              ! [input]
      double precision dpmo(nb, nb, ns, 3)    ! [input]
      double precision cf(n,n)                ! [in/output]
      double precision pbd(3,3,nch)           ! [in/output]

      integer i, pid(n), l_sc, k_sc, j, k, i_sc, i_
      double precision dd(n), pb_(3,3,n)

      
      call wfat_get_newpolb(n, orbid0, ispin, cf, nb, ns, dpmo, pb_)

      do i = 1, n
         dd(i) = 0.0d0
         do j = 1, 3
            dd(i) = dd(i) + (pb_(j,j,i) - pbd(j,j,i))**2
         enddo
      enddo
      
      do i = 1, n
         pid(i) = i
      enddo
      call wfat_dsortx(dd, 1, n, pid)


      if (.not. MA_alloc_get(MT_DBL, n**2, 'cf', l_sc, k_sc))
     &     call errquit('WFAT: Cannot allocate cf', 41, MA_ERR)
      call dcopy(n**2, cf, 1, dbl_mb(k_sc), 1)

      do i = 1, n
         j = pid(i)
         k = mfid(i)
         pbd(:,:,k) = pb_(:,:,j)
         i_sc = k_sc + n*(j-1)
         cf(:,i) = dbl_mb(i_sc:i_sc+n-1)
      enddo

cOLD      do i = 1, n
cOLD         i_ = mfid(i)
cOLD         do k = 1, 3
cOLD            call wfat_submatrix(nb, nb, dpmo(:,:,ispin,k), n, n, orbid0,
cOLD     &           orbid0, dbl_mb(k_sc))
cOLD            call wfat_1e_expval_1eopr(n, cf(:,i), dbl_mb(k_sc),
cOLD     &           mmd(k,i_))
cOLD         enddo
cOLD         mmd(:,i_) = mmd(:,i_) - ctr(:,i_)
cOLD      enddo
      
      
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: cf deallocation failed', 20, MA_ERR)

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_get_newpolb(n, orbid0, ispin, cf, nb, ns, dpmo,
     &     pb)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer n
      integer orbid0(n)
      integer ispin
      double precision cf(n,n)
      integer nb
      integer ns
      double precision dpmo(nb, nb, ns, 3)
      double precision pb(3,3,n)

      integer l_dp, k_dp, is, i, j, i_dp1, i_dp2, k, io

      
      if (.not. MA_alloc_get(MT_DBL, nb*n*3, 'dpo_', l_dp, k_dp))
     &     call errquit('WFAT: Cannot allocate dpo_', 41, MA_ERR)


c     ==== Construct the columns of dpmo corresponding to ====
c     ====   the new degenerate channels (i.e. when the   ====
c     ====     kets are the new degenerate channels)      ====
      is = ispin  
      do i = 1, 3      ! Loop over cartesian coordinates
         do j = 1, n      ! Loop over the new degenerate channels
            i_dp1 = k_dp + nb*n*(i-1) + nb*(j-1)
            i_dp2 = i_dp1 + nb - 1

            dbl_mb(i_dp1:i_dp2) = 0.0d0
            do k = 1, n       ! Loop (sum) over the old degenerate channels
               io = orbid0(k)
               dbl_mb(i_dp1:i_dp2) = dbl_mb(i_dp1:i_dp2) +
     &              cf(k,j)*dpmo(:,io,is,i)
            enddo   
         enddo
      enddo


c     ==== Deliberately set the elements of dpmo between the ====
c     ====  new degenerate channels (i.e. when both the bra  ====
c     ====    and ket are the new degenerate channels) to    ====
c     ====                 Infinity or NaN                   ====           ! 1)
      do i = 1, 3
         do j = 1, n
            do k = 1, n
               io = orbid0(k)
               i_dp1 = k_dp + nb*n*(i-1) + nb*(j-1) + (io-1)
               dbl_mb(i_dp1) = 1/0.0d0
            enddo
         enddo
      enddo
      

c     ==== Compute the new static polarizabilities ====
c     ====     for these degenerate channels       ====
      call wfat_calc_polb_eco(n, orbid0, ispin, dbl_mb(k_dp), pb)


      if (.not. MA_free_heap(l_dp))
     &     call errquit('WFAT: dpo_ deallocation failed', 20, MA_ERR)
      end
c     Notes:
c     1) This is because the elements of dpmo between the new degenerate
c        channels will not be used in the computation of the static
c        polarizabilities. Assigining an Infinity or NaN to these elements
c        helps in identifying when they are mistakenly used in the above
c        computations.
c============================================================
      

c============================================================
c============================================================
      subroutine wfat_calc_lfpolb(protmat, pb, lfpb)

      implicit none
      double precision protmat(3,3)
      double precision pb(3,3)
      double precision lfpb

      call wfat_1e_expval_1eopr(3, protmat(:,3), pb, lfpb)
      lfpb = -2 * lfpb
      
      end
c============================================================
      
      
