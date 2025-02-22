c=======================================================
c=======================================================
      subroutine wfat_get_rdm(rdm)

      implicit none
#include "stdio.fh"      
#include "global.fh"      
#include "mafdecls.fh"
#include "errquit.fh"      
c#include "cscf.fh"
#include "wfat_qm_param.fh"
      double precision rdm(nbas0, nbas0, nspin0)

      integer is, io0, io, ir0
      integer max_ne, nrd0
      integer l_r0, k_r0


c cscf      nel(1) = nalpha
c cscf      nel(2) = nbeta

      max_ne = maxval(nel)
      nrd0 = nbas0 * max_ne
      if (.not. MA_alloc_get(MT_DBL, nrd0, 'rdm0', l_r0, k_r0))
     &     call errquit('WFAT: rdm0 allocation failed', 59, MA_ERR)

      do is = 1, nspin0
c        ==== Construct the coefficient array rdm0 ====
         do io = 1, nel(is)
            ir0 = k_r0 + (io-1) * nbas0
            call ga_get(g_mov(is), 1, nbas0, io, io, dbl_mb(ir0),
     &                  nbas0)
         enddo
         call wfat_mult_rdm0(nbas0, nel(is),
     &        dbl_mb(k_r0:k_r0+nbas0*nel(is)-1),
     &        dbl_mb(k_r0:k_r0+nbas0*nel(is)-1),
     &        rdm(:,:,is))
      enddo

      if (.not. ma_free_heap(l_r0))
     &   call errquit('WFAT: rdm0 deallocation failed', 1, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_mult_rdm0(nbas, ne, rdm0_l, rdm0_r, rdm)

      implicit none
      integer nbas, ne
      double precision rdm0_l(nbas, ne), rdm0_r(nbas, ne)
      double precision rdm(nbas, nbas)


      rdm = matmul(rdm0_l, transpose(rdm0_r))

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_dens_1ptbas(nbas, nspin, rdm, maxd, 
     &                            bval, d1bval, d2bval, 
     &                            den, d1den, d2den)

      implicit none
#include "errquit.fh"      
      integer nbas                            ! [input], no. of basis
      integer nspin                           ! [input], no. of spin
      double precision rdm(nbas, nbas, nspin) ! [input], reduced density matrix
      integer maxd                            ! [input], max. derivative order requested.  3)
      double precision bval(nbas)             ! [input], basis
      double precision d1bval(3, nbas)        ! [input], 1st derivatives of basis.    1)
      double precision d2bval(6, nbas)        ! [input], 2nd derivatives of basis.    2)
      double precision den(nspin)             ! [output], density
      double precision d1den(3, nspin)        ! [output], 1st derivatives of density.   1)
      double precision d2den(6, nspin)        ! [output], 2nd derivatives of density.   2)

      integer is, ibl, ibk, ic, ic1, ic2
      integer ff
      double precision den_lk


      if (maxd < 0 .or. maxd > 2) then
         call errquit('WFAT: Argument maxd of wfat_dens_1ptbas '//
     &        'must not be negative and (for the moment) must not '//
     &        'exceed 2.')
      endif


      den = 0.0d0
      d1den = 0.0d0
      d2den = 0.0d0
      do is = 1, nspin

         do ibk = 1, nbas
            do ibl = 1, ibk

               if (ibk == ibl) then
                  ff = 1
               else
                  ff = 2
               endif
               
               if (nspin == 1) then
                  den_lk = 2 * rdm(ibl, ibk, 1)
               else if (nspin == 2) then
                  den_lk = rdm(ibl, ibk, is)
               else
                  call errquit('WFAT: Argument nspin of '//
     &                 'wfat_dens_1ptbas must be either 1or 2.',
     &                 1, UERR)
               endif

c              ==== Density ====
               den(is) = den(is) + ff * den_lk * bval(ibk) * bval(ibl)

c              ==== 1st derivatives of density ====
               if (maxd >= 1) then
                  do ic = 1, 3
                     d1den(ic, is) = d1den(ic, is) +
     &                               ff * den_lk *
     &                               (bval(ibk) * d1bval(ic, ibl) +
     &                               bval(ibl) * d1bval(ic, ibk))
                  enddo
               endif

c              ==== 2nd derivatives of density ====
               if (maxd >= 2) then
                  do ic = 1, 6
                     call wfat_get_i2cart(ic, ic1, ic2)
                     d2den(ic, is) = d2den(ic, is) +
     &                               ff * den_lk * 
     &                               (bval(ibk) * d2bval(ic,ibl) +
     &                               d1bval(ic1,ibk) * d1bval(ic2,ibl) +
     &                               d1bval(ic2,ibk) * d1bval(ic1,ibl) +
     &                               bval(ibl) * d2bval(ic,ibk))
                  enddo
               endif

            enddo
         enddo

      enddo


      end

c     NOTES:
c       1) For a given column, the order of the partial derivatives is the following:
c          x, y, z.
c       2) For a given column, the order of the partial derivatives is the following:
c          xx, yy, zz, xy, yz, xz.
c       3) When maxd=0, only d1bval and d2bval are not used. When maxd=1, d2bval is
c          not used. When maxd=2, bval, d1bval, and d2bval are all used.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_gamma_1pt(nspin, d1dens, d2dens, gamma, d1gamma)
      
      implicit none
#include "wfat_indexing.fh"      
      integer nspin                      !< [Input] The number of spin channels
      double precision d1dens(3, nspin)  !< [Input] The density gradient
      double precision d2dens(6, nspin)  !< [Input] The density curvature
      double precision gamma(3)          !< [Output] The density gradient norm
      double precision d1gamma(3,3)      !< [Output] The first derivative of density gradient norm

      integer id1, id2, id
      

c     ==== Calculate gamma ====
      call calc_rgamma(nspin, 1, d1dens, gamma)


c     ==== Calculate 1st derivatives of gamma ====
      d1gamma = 0.0d0
      d1gamma(1:3, G_AA) = 0.0d0
      do id1 = 1, 3   ! just element identification
         do id2 = 1, 3   ! to be summed
            call wfat_get_cart2i(id1, id2, id)
            d1gamma(id1, G_AA) = d1gamma(id1, G_AA) +
     &                           d2dens(id, 1) * d1dens(id2,1) +
     &                           d1dens(id2,1) * d2dens(id, 1)
         enddo
      enddo

c     ==== If spin-polarzied ====
      if (nspin == 2) then
         d1gamma(1:3, G_AB) = 0.0d0
         do id1 = 1, 3   ! just element identification
            do id2 = 1, 3   ! to be summed
               call wfat_get_cart2i(id1, id2, id)
               d1gamma(id1, G_AB) = d1gamma(id1, G_AB) +
     &                              d2dens(id, 1) * d1dens(id2,2) +
     &                              d1dens(id2,1) * d2dens(id, 2)
            enddo
         enddo

         d1gamma(1:3, G_BB) = 0.0d0
         do id1 = 1, 3   ! just element identification
            do id2 = 1, 3   ! to be summed
               call wfat_get_cart2i(id1, id2, id)
               d1gamma(id1, G_BB) = d1gamma(id1, G_BB) +
     &                              d2dens(id, 2) * d1dens(id2,2) +
     &                              d1dens(id2,2) * d2dens(id, 2)
            enddo
         enddo
      endif

      if (nspin == 1) then
         gamma(2:3) = 0.0d0
         d1gamma(:, 2:3) = 0.0d0
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_tau_1pt(nbas, nspin, rdm, d1bval, tau)

      implicit none
      integer nbas
      integer nspin
      double precision rdm(nbas, nbas, nspin)
      double precision d1bval(3, nbas)
      double precision tau(nspin)

      integer is, j, i, id
      integer t
      double precision dd


      tau = 0.0d0
      do j = 1, nbas
         do i = 1, j
            if (i /= j) then
               t = 2
            else
               t = 1
            endif

            dd = 0.0d0
            do id = 1, 3
               dd = dd + d1bval(id, j) * d1bval(id, i)
            enddo

            do is = 1, nspin
               tau(is) = tau(is) + t * rdm(j, i, is) * dd
            enddo

         enddo
      enddo
      tau = 0.5d0 * tau

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_lapl_1pt(nbas, nspin, d2dens, lapl)

      implicit none
      integer nbas
      integer nspin
      double precision d2dens(6, nspin)
      double precision lapl(nspin)

      integer is, ic, ic1, ic2
      

      lapl = 0.0d0
      do is = 1, nspin
         do ic = 1, 6
            call wfat_get_i2cart(ic, ic1, ic2)

            if (ic1 == ic2) then
               lapl(is) = lapl(is) + d2dens(ic, is)
            endif
         enddo
      enddo

      end
c=======================================================

      
c=======================================================      
c=======================================================
      subroutine wfat_vecdens_numinteg_xy(geom, bas, nz, zp, nx, xlims,
     &           ny, ylims, vcf, g_rh)
c     DESCRIPTION :
c     Numerically calculate the orbital density integrated along x and y axes.
c     The orbital is given by the coefficient vcf.
      
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "bas.fh"
#include "stdio.fh"
#include "global.fh"      
      integer geom                    ! [input], geometry handle using which bas handle was created.
      integer bas                     ! [input], basis handle in which vcf is represented.
      integer nz                      ! [input], number of z points to evaluate the integrated density.
      double precision zp(nz)         ! [input], the z points at which to evaluate the integrated density.
      integer nx                      ! [input], number of integration points along x.
      double precision xlims(2)       ! [input], integration limits along x.
      integer ny                      ! [input], number of integration points along y.
      double precision ylims(2)       ! [input], integration limits along y.
      double precision vcf(1:*)       ! [input], coefficients of MO vector, dimension [nb].
      integer g_rh                    ! [output], the xy-integrated density.
c                                       MUST BE declared such that IT HAS ONLY ONE COLUMN.

      integer nb, nz0
      integer ndid, iz_in, iz_fin, ix, iy, iz, i, col1, col2
      integer l_d, k_d, l_rh, k_rh, l_ov, k_ov


      call ga_sync()

      ndid = ga_nodeid()
      if (.not. bas_numbf(bas, nb))
     &     call errquit('WFAT: Getting nb failed in '//
     &     'wfat_vecdens_numinteg_xy', 1, BASIS_ERR)
      
      
c     ==== Obtain MO-resolved density matrix ====
      if (.not. ma_alloc_get(mt_dbl, nb*nb, 'mo_rdm', l_d, k_d))
     &     call errquit('WFAT: Cannot allocate mo_rdm', 1, MA_ERR)
      call wfat_project_mat(nb, 1, vcf(1:nb), dbl_mb(k_d))
      call wfat_dtranspose(nb, nb, dbl_mb(k_d), dbl_mb(k_d))


c     ==== Obtain range of indices of zp for the present process ====
      call ga_distribution(g_rh, ndid, iz_in, iz_fin, col1, col2)
      if ((col1 /= 1) .or. (col1 /= col2))
     &     call errquit('WFAT: Argument g_rh of wfat_vecdens_numinteg'//
     &     '_xy must be declared such that it has one column.', 1, UERR)


c     ==== Calculate xy-integrated density ====
      if (.not. ma_alloc_get(mt_dbl, nb*nb, 'overlap_z', l_ov, k_ov))
     &     call errquit('WFAT: Cannot allocate overlap_z', 1, MA_ERR)
      nz0 = iz_fin - iz_in + 1
      if (.not. ma_alloc_get(mt_dbl, nz0, 'rho_z', l_rh, k_rh))
     &     call errquit('WFAT: Cannot allocate rho_z', 1, MA_ERR)

      do iz = iz_in, iz_fin
         call wfat_overlap_numinteg_xy(geom, bas, nx, xlims, ny, ylims, 
     &        zp(iz), dbl_mb(k_ov))
         i = iz - iz_in + 1
         dbl_mb(k_rh+i-1) = dot_product(dbl_mb(k_d:k_d+nb*nb-1),
     &                                  dbl_mb(k_ov:k_ov+nb*nb-1))

cforcheck         write(luout, '(" ndid=", i3, " iz=", i6, '//
cforcheck     &        '" k_d=", es14.6, " =k_ov", es14.6, " k_rh=", es14.6)')
cforcheck     &        ndid, iz, sum(abs(dbl_mb(k_d:k_d+nb*nb-1))),
cforcheck     &        sum(abs(dbl_mb(k_ov:k_ov+nb*nb-1))),
cforcheck     &        dbl_mb(k_rh+i-1)
      enddo
      call ga_put(g_rh, iz_in, iz_fin, 1, 1, dbl_mb(k_rh), nz0)


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_d))
     &     call errquit('WFAT: Cannot deallocate mo_rdm', 1, MA_ERR)
      if (.not. ma_free_heap(l_ov))
     &     call errquit('WFAT: Cannot deallocate overlap_z', 1, MA_ERR)
      if (.not. ma_free_heap(l_rh))
     &     call errquit('WFAT: Cannot deallocate rho_z', 1, MA_ERR)

      call ga_sync()

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_rdmas(nb, ns, mov, vpt, asc, ras, rdmas)

      implicit none
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer nb
      integer ns
      double precision mov(nb, nb, ns)    ! [input] all MO's.
      double precision vpt(3)             ! [input] evaluation point.
      double precision asc(3)             ! [input] coordinate of asymptotic center.
      double precision ras                ! [input] radius of asymptotic sphere centered at asc.
      double precision rdmas(nb, nb, ns)  ! [output] asymptotically modified RDM.

      integer is, i, j, io, ieg
      double precision rr, kappa, erg0
      

      do is = 1, ns
         do j = 1, nb
            do i = 1, j

               rdmas(i,j,is) = 0.0d0
               do io = 1, nel(is)
                  rr = sqrt(sum((vpt-asc)**2))
                  if (rr-ras < 0.0d0)
     &                 call errquit('WFAT: Evaluation point is still '//
     &                 'not within asymptotic criteria', 1, UERR)
                  ieg = k_moerg + (is-1)*nbas0 + io - 1
                  erg0 = dbl_mb(ieg)
                  kappa = sqrt(2*abs(erg0))
                  rdmas(i,j,is) = rdmas(i,j,is) +
     &                 mov(j,io,is)*mov(i,io,is)*exp(-kappa*(rr-ras))
               enddo
               rdmas(j,i,is) = rdmas(i,j,is)

            enddo
         enddo
      enddo

      end
c=======================================================
