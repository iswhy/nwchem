c=================================================
c=================================================
      subroutine mewfat_dfterg(geom, bas, ns, gq, mcel, rdm, xcerg)

      implicit none
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "errquit.fh"
      integer geom                ![input], 
      integer bas                 ![input],
      integer ns
      integer gq                  ![input], 
      integer mcel                ![input], 
      double precision rdm(1:*)   ![input], 
      double precision xcerg      ![output], 

      integer ndid
      integer nat, nb, natomp, nptp, max_npt
      integer g_xce
      integer l_icp, l_ipp, l_pt, l_wt
      integer k_icp, k_ipp, k_pt, k_wt
      double precision xce


      call ga_sync()
      ndid = ga_nodeid()
      if (.not. ga_create(mt_dbl, 1, 1, 'xcerg', 0, 0, g_xce))
     &     call errquit('mewfat_dfterg: Cannot create xcerg global '//
     &     'array.', 1, GA_ERR)
      call ga_zero(g_xce)
      

c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
c     ==== Get number of bases ====
      if (.not. bas_numbf(bas, nb))
     &     call errquit('WFAT: Getting nbas failed in ', 53, UERR)
      

c     ==== Obtaining grid ranges for the current process ====
      call wfat_alloc_gridrange(nat, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natomp, nptp)

      
c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gq, max_npt, l_pt, k_pt, l_wt, k_wt)


c     ==== Perform quadrature summation ====
      call mewfat_dfterg_qsum(geom, bas, nb, ns, natomp, int_mb(k_icp),
     &     int_mb(k_ipp), gq, max_npt, dbl_mb(k_pt), dbl_mb(k_wt),
     &     mcel, nat, rdm, xce)
      

c     ==== Accumulate to g_xce ====
      call ga_acc(g_xce, 1, 1, 1, 1, xce, 1, 1.0d0)
      call ga_sync()
      call ga_get(g_xce, 1, 1, 1, 1, xcerg, 1)
      

c     ==== Deallocations ====
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      if (.not. ga_destroy(g_xce))
     &     call errquit('mewfat_dfterg: Cannot destroy xcerg.',
     &     1, GA_ERR)
      
      end
c=================================================

      
c=================================================
c=================================================      
      subroutine mewfat_dfterg_qsum(geom, bas, nb, ns, natp, icp, ipp,
     &     gq, max_npt, qpt, qwt, mcel, nat, rdm, xce)

      implicit none
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
      integer geom                        ! [input], geometry handle.
      integer bas                         ! [input], basis handle.
      integer nb                          ! [input], no. of basis functions in bas.
      integer ns
      integer natp                        ! [input], no. of centers for the current process returned by wfat_get_gridrange.
      integer icp(natp)                   ! [input], centers for the current process returned by wfat_get_gridrange.
      integer ipp(2,natp)                 ! [input], range of points of each centers for the current process returned by wfat_get_gridrange.
      integer gq                          ! [input], grid quality.
      integer max_npt                     ! [input], max. no. of points for the given gq as returned by wfat_get_gridmaxpt.
      double precision qpt(3, max_npt)    ! [output/scratch], array Quadrature points
      double precision qwt(max_npt)       ! [output/scratch], array Quadrature weights
      integer mcel                        ! [input], the m parameter related to the cell function. Usually taken as wf_msm from wfat.fh.
      integer nat                         ! [input], no. of atoms in the molecule.
      double precision rdm(nb,nb,ns)      ! [input], reduced density matrix.
      double precision xce                ! [output],

      character*16 tag
      character*4 lv
      logical t1
      integer ndid, dbas, npt, ic, icent, nb0, nb1, nb2
      integer i, imf, l_b, k_b, i_b0, i_b1, i_b2
      double precision atp(3), zn, v_m(3), xce0, celf, dbl_dum

      parameter (imf = 0)

      logical xc_chktau
      external xc_chktau
      

      ndid = ga_nodeid()
      call ga_sync()


c     ==== Determine DFT XC level ====
      if (xc_chktau()) then
         call errquit('MEWFAT: At the moment, mGGA functional is not '//
     &                'yet supported in WFAT.', 11, UERR)
      endif
      call wfat_get_xclvl(dbas, lv)


c     ==== Allocate basis and its derivatives ====
      nb0 = nb
      nb1 = 1
      nb2 = 1
      if (dbas >= 1) nb1 = 3*nb
      if (dbas >= 2) nb2 = 6*nb
      if (.not. ma_alloc_get(mt_dbl, nb0+nb1+nb2, 'basis', l_b, k_b))
     &     call errquit('MEWFAT: Cannot allocate basis', 1, MA_ERR)
      i_b0 = k_b
      i_b1 = i_b0 + nb0
      i_b2 = i_b1 + nb1
      
      
c     ==== BEGIN QUADRATURE SUMMATION ====
c     ==== Loop over atoms ====
      call ga_sync()
      xce = 0.0d0              ! intg_(nsh,nqn,nat)     nch=1
      do ic = 1, natp
         icent = icp(ic)         ! icp(natp)
         if (.not. geom_cent_get(geom, icent, tag, atp, zn))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
         call wfat_get_gridpt(int(zn), gq, max_npt, npt, qpt, qwt)   ! qpt(3,max_npt)   qwt(max_npt)
         
c        ==== Loop over all quadrature points around the current atom ====
         do i = ipp(1,ic), ipp(2,ic)           ! ipp(2,natp)
            call geom_cent2cent(geom, icent, qpt(:,i), imf, v_m, t1)

c           ==== Calculate basis and its derivatives ====
            call wfat_eval_basis_1pt(geom, bas, dbas, v_m,
     &                               dbl_mb(i_b0), dbl_mb(i_b1),
     &                               dbl_mb(i_b2), dbl_dum)
            
c           ==== Calculate XC energy integrand ====
            call mewfat_xcerg_1ptbas(nb, ns, dbl_mb(i_b0), dbl_mb(i_b1),
     &           dbl_mb(i_b2), rdm, lv, xce0)
            
c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom, icent, qpt(:,i), mcel, celf)
            xce = xce + xce0 * celf * qwt(i)
         enddo
      enddo

      
      call ga_sync()
      if (.not. ma_free_heap(l_b))
     &     call errquit('MEWFAT: Cannot deallocate basis', 1, MA_ERR)
      
      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_xcerg_1ptbas(nb, ns, bval, d1bval, d2bval, rdm,
     &     lv, xce)

      implicit none
#include "stdio.fh"
#include "wfat_indexing.fh"
#include "errquit.fh"
      integer nb                       ! [input], no. of basis.
      integer ns
      double precision bval(nb)        ! [input], basis functions evaluated at the desired point.
      double precision d1bval(3,nb)    ! [input], 1st derivative of basis functions evaluated at the desired point.   1)
      double precision d2bval(6,nb)    ! [input], 2nd derivative of basis functions evaluated at the desired point.   1)
      double precision rdm(nb,nb,ns)   ! [input], reduced density matrix.
      character*4 lv                   ! [input],      1)
      double precision xce             ! [output],

      integer maxd
      double precision dens(ns), d1dens(3,ns), d2dens(6,ns)
      double precision gamma(3), d1gamma(3,3), tau(ns)
      double precision dfndr(ns), dfndr2(NCOL_AMAT2)
      double precision dfndg(3), dfndg2(NCOL_CMAT2)
      double precision dfndt(ns), dfndt2(NCOL_MMAT2)

      
      if (lv /= 'lda' .and. lv /= 'gga' .and. lv /= 'mgga'
     &    .and. lv /= 'none')
     &     call errquit('MEWFAT: Argument lv of mewfat_xcerg_1ptbas '//
     &                  'must be either lda, gga, mgga, or '//
     &                  'none.', 1, UERR)
      if (lv == 'none') then
         xce = 0.0d0
         return
      endif

      
      if (lv == 'lda') maxd = 0
      if (lv == 'gga' .or. lv == 'mgga') maxd = 2
      call wfat_dens_1ptbas(nb, ns, rdm, maxd,
     &                      bval, d1bval, d2bval,
     &                      dens, d1dens, d2dens)

      gamma = 0.0d0
      tau = 0.0d0
      xce = 0.0d0
      dfndr = 0.0d0
      dfndr2 = 0.0d0
      dfndg = 0.0d0
      dfndg2 = 0.0d0
      dfndt = 0.0d0
      dfndt2 = 0.0d0
      if (lv == 'lda') then
         call nwxc_eval_df(ns, 1, dens, gamma, tau, xce, dfndr,
     &                     dfndg, dfndt)         
      else if (lv == 'gga' .or. lv == 'mgga') then
         call wfat_gamma_1pt(ns, d1dens, d2dens, gamma, d1gamma)
ccheck         if (lv == 'mgga') 
ccheck  &           call wfat_tau_1pt(nb, ns, rdm, d1bval, tau)
         call wfat_tau_1pt(nb, ns, rdm, d1bval, tau)
         
         call nwxc_eval_df2(ns, 1, dens, gamma, tau, xce, dfndr, 
     &                      dfndr2, dfndg, dfndg2, dfndt, dfndt2)
      endif

      end
c     Notes:
c       1) d1bval and d2bval are used only when lv is equal to 'gga' or 'mgga'.      
c=================================================
