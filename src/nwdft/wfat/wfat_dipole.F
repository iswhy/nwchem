c=======================================================
c=======================================================
      subroutine wfat_calc_dpmo(ctr, dpmo)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      double precision ctr(3)   ! [input],
      double precision dpmo(nbas0, nbas0, nspin0, 3)   ! [output],

      integer NMMT
      parameter(NMMT=4)
      integer nmp, ncf, is, i, iq, l_mp, k_mp, i_mp, l_cf, k_cf, i_cf


c     ==== Dipole matrices in GTO basis relative to ctr ====
      nmp = nbas0 * nbas0 * NMMT
      if (.not. MA_alloc_get(MT_DBL, nmp, 'mpmat', l_mp, k_mp))
     &     call errquit('WFAT: mpmat allocation failed', 1,
     &     MA_ERR)

      call int_mpole_all(bas_hd, nbas0, 1, ctr, dbl_mb(k_mp))
cc      write(luout, '("X dipole AO:")')
cc      call wfat_dprintmat(nbas0, nbas0, 'f', dbl_mb(k_mp+nbas0**2))
cc      write(luout, '("Y dipole AO:")')
cc      call wfat_dprintmat(nbas0, nbas0, 'f', dbl_mb(k_mp+2*nbas0**2))
cc      write(luout, '("Z dipole AO:")')
cc      call wfat_dprintmat(nbas0, nbas0, 'f', dbl_mb(k_mp+3*nbas0**2))


c     ==== Get MO coeffs ====
      ncf = nbas0 * nbas0 * nspin0
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vcf', l_cf, k_cf))
     &     call errquit('WFAT: vcf allocation failed', 1, MA_ERR)
      do is = 1, nspin0
         i_cf = k_cf + (is-1)*nbas0*nbas0
         call ga_get(g_mov(is), 1, nbas0, 1, nbas0, dbl_mb(i_cf), 
     &        nbas0)
      enddo
      

c     ==== Dipole matrices in MO basis relative to ctr ====
      do i = 1, 3
         iq = i + 1
         i_mp = k_mp + (iq-1)*nbas0*nbas0
         do is = 1, nspin0
            i_cf = k_cf + (is-1)*nbas0*nbas0
            call wfat_represent_sym(nbas0, nbas0, dbl_mb(i_cf),
     &           dbl_mb(i_mp), dpmo(:,:,is,i))
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_mp))
     &     call errquit('WFAT: mpmat deallocation failed', 38, MA_ERR)
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: vcf deallocation failed', 38, MA_ERR)
      
      end
c=======================================================
      

c===============================================================
c===============================================================
      subroutine wfat_dpmo_dgn_sub(nch, orbid, nsub, maxn, dgs, dgss,
     &     ndgs, nb, ns, dpmo, dps)
c     DESCRIPTION:
c     nsub = The number of degenerate subspaces, a nondegenerate state is considered a
c            degenerate state in a subspace of dimension 1.
c     maxn = The largest dimension of all degenerate subspaces, i.e. maxn=maxval(ndgs).
c     dgs = The channel ID arranged into degenerate subspaces (columns). Each subspace
c           (column) must only have states of the same spin. However, note that dgs has
c           no information about the spins of its columns. This information is not needed
c           because dgs contains the indices of the channels which already encode the spin
c           of that channel. The range of its elements is 1 <= dgs(i,j) <= nch. Its 
c           minimum size ismaxn*nsub.
c     dgss = The spin of the degenerate subspaces.
c     ndgs = The dimension of each degenerate subspace. Its minimum size is nsub.
c     dps = The dipole matrix in molecular frame in each degenerate subspace. The 3
c           indices in the 3rd dimension correspond to x, y, z. Its minimum size is
c           maxn*maxn*6*nsub.

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer nch                              ! [input]
      integer orbid(nch)                       ! [input]
      integer nsub                             ! [input]
      integer maxn                             ! [input]
      integer dgs(maxn,nsub)                   ! [input]
      integer dgss(nsub)                       ! [input]
      integer ndgs(nsub)                       ! [input]
      integer nb                               ! [input]
      integer ns                               ! [input]
      double precision dpmo(nb,nb,ns,3)        ! [input]
      double precision dps(maxn,maxn,3,nsub)   ! [output]

      integer i, is, nn, j, jj, idv(maxn)

      if (maxn < maxval(ndgs)) then
         write(luout, '("maxn = ", i4)') maxn
         write(luout, '("max. ndgs = ", i4)') maxval(ndgs)
         call errquit('WFAT: Incompatibility '//
     &        'between maxn and ndgs detected.', 28, UERR)
      endif
      
      dps = 0.0d0
      do i = 1, nsub
         is = dgss(i)
         nn = ndgs(i)

         do j = 1, nn
            jj = dgs(j,i)         ! 1 <= jj <= nch
            idv(j) = orbid(jj)    ! 1 <= idv(j) <= nb
         enddo
         do j = 1, 3
            call wfat_submatrix(nb, nb, dpmo(:,:,is,j), nn, nn, idv, 
     &           idv, dps(1:nn,1:nn,j,i))
         enddo
      enddo

      end
c===============================================================


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

      
