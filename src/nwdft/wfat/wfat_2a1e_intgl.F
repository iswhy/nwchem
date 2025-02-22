c===============================================================
c===============================================================
      subroutine wfat_2a1e_intgl(rtdb, geom, nqn, nid, mid, lmax_,
     &     gq, zi, celm, sphi, mch, nr, nt, ke, vcf, erg, chmmt, chpos,
     &     pottype, nel, potparam, wfat1, nlgl, nr_m, nt_m, nr_p, nt_p, 
     &     cfx, cfy, cfz, wfat1_exc, fedvr, nfe, ndvr, maxn, xq, wq, 
     &     ismask, mask_param, g_il, g_j, g_k)

      implicit none
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "geom.fh"
      integer rtdb                    ! [input],
      integer geom                    ! [input],
      integer nqn                     ! [input],
      integer nid(nqn)                ! [input],
      integer mid(nqn)                ! [input],
      integer lmax_                   ! [input],
      integer gq                      ! [input],
      double precision zi             ! [input],
      integer celm                    ! [input],
      integer sphi                    ! [input],
      integer mch                     ! [input],
      integer nr                      ! [input],
      integer nt                      ! [input],
      double precision ke             ! [input],
      double precision vcf(nr*nt)     ! [input],
      double precision erg            ! [input],
      double precision chmmt(3)       ! [input],
      double precision chpos(3)       ! [input],
      character*3 pottype             ! [input]
      integer nel                     ! not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)    ! [input], 1)
      logical wfat1                   ! [input],
      integer nlgl                    ! [input], only used when wfat1=.true.
      integer nr_m                    ! [input], only used when wfat1=.true.
      integer nt_m                    ! [input], only used when wfat1=.true.
      integer nr_p                    ! [input], only used when wfat1=.true.
      integer nt_p                    ! [input], only used when wfat1=.true.
      double precision cfx(1:*)       ! [input], only used when wfat1=.true.
      double precision cfy(1:*)       ! [input], only used when wfat1=.true.
cOLD      double precision cfz(nr*nt) ! [input], only used when wfat1=.true.
      double precision cfz(1:*)       ! [input], only used when wfat1=.true.
      logical wfat1_exc(1:*)          ! [input], only used when wfat1=.true., min. size is nqn.
      logical fedvr                   ! [input],
      integer nfe                     ! [input], only used when fedvr=.true.
      integer ndvr(*)                 ! [input], only used when fedvr=.true., size is nfe
      integer maxn                    ! [input], only used when fedvr=.true.
      double precision xq(*)          ! [input], only used when fedvr=.true., size is maxn*nfe
      double precision wq(*)          ! [input], only used when fedvr=.true., size is maxn*nfe
      logical ismask                  ! [input]
      double precision mask_param(2)  ! [input]. Only used when ismask=.true.. Must always be supplied regardless.
      integer g_il                    ! [output],
      integer g_j                     ! [output], only calculated when wfat1=.true.
      integer g_k(2)                  ! [output], only calculated when wfat1=.true.

      character*12 dt_str
      character*16 tag
      logical t1, cgrid(1)
      integer ndid, natom, max_npt, nscr, nsh, nom, ni0, ni_, nil, ninc
      integer nij0, nij_, nij, njc, nik1_0, nik1_, nik1, nk1c, nik2_0
      integer nik2_, nik2, nk2c, natomp, nptp, pr_chk(10), g_nptx, nptx
      integer jpt, pr_count, kk, ozn, znuc, geomc(1)
      integer ic, nom1, nom2, i_om1, i_om2, nnnr, nnnt
      integer icent, npt, ipt_i, ipt_f, ipt, ipt0, i_l, i_r
      integer l_pt, l_wt, l_sc, l_om, l_i0, l_i_, l_il, l_inc, l_icp
      integer l_ipp, l_ij0, l_ij_, l_ij, l_jc, l_ik1_0, l_ik1_, l_ik1
      integer l_k1c, l_ik2_0, l_ik2_, l_ik2, l_k2c, l_om1, l_om2, l_dr
      integer l_dt
      integer k_pt, k_wt, k_sc, k_om, k_i0, k_i_, k_il, k_inc, k_icp
      integer k_ipp, k_ij0, k_ij_, k_ij, k_jc, k_ik1_0, k_ik1_, k_ik1
      integer k_k1c, k_ik2_0, k_ik2_, k_ik2, k_k2c, k_om1, k_om2, k_dr
      integer k_dt
      double precision npt_tot, t_now, zn
      double precision t_last, t_start, celf(1), atp(3), v_c(3)
      double precision v_m(3), vc_psi, wtcelf(1), dt_prog, dt_all
      double precision normch, vc_cpsi(3), ovl1
      double precision polb_
      double complex czero, cone
      integer nch, nch1, imf
      parameter(nch=1, nch1=3, imf=0)
      parameter(czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)),
     &          cone = cmplx(1.0d0, 0.0d0, kind=kind(1.0d0)))


      call ga_sync()
      ndid = ga_nodeid()


c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
      
      
c     ==== Some tests ====
      call wfat_sometests(rtdb, zi, nch, erg, lmax_, chmmt)


c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gq, max_npt, l_pt, k_pt, l_wt, k_wt)


c     ==== Allocate WFAT integral, its integrands, and some other ====
      if (wfat1) then
         if (fedvr) then
            nscr = maxn**2 + maxn +
     &           (nr + nt + nr*nt) +
     &           (nt_p + nr*nt_p)
            if (mch >= 1) nscr = nscr + (nt_m + nr*nt_m)
         else
            nscr = 
     &           (nr + nt + nr*nt) +
     &           (nr_p + nt_p + nr_p*nt_p)
            if (mch >= 1) nscr = nscr + (nr_m + nt_m + nr_m*nt_m)
         endif
      else
         nscr = 
     &        nr +
     &        nt +
     &        nr*nt
         if (fedvr) nscr = nscr + maxn**2 + maxn
      endif
      if (.not. MA_alloc_get(MT_DBL, nscr, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 59, MA_ERR)
      call wfat_alloc_intgl(lmax_, nqn, nch, natom, nsh, nom, ni0, ni_,
     &     nil, ninc, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il,
     &     l_inc, k_inc)


c     ==== Obtaining grid ranges for the current process ====
      if (.not. rtdb_get(rtdb, 'wfat:ozn', mt_int, 1, ozn))
     &     call errquit('WFAT: Cannot read ozn.', 73, RTDB_ERR)
      call wfat_alloc_gridrange(natom+1, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange_c(geom, ozn, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natomp, nptp)
      call wfat_need_cgrid(rtdb, geom, nch, chpos, cgrid)
      call wfat_cgeom_create(geom, ozn, nch, chpos, geomc)
      call wfat_print_gridpt(rtdb, geom, gq, ozn, nch, cgrid, chpos)


c     ==== Calculate WFAT(1) Omega ====
      if (wfat1) then
c        nch1 = the number of channels that contains the correction terms (needed for
c               WFAT(1)), it is equal to 3 because the ionized channel has three correction
c               terms associated to the three cartesian components of the dipole operator.
         call wfat_alloc_intgl_wfat1(lmax_, nqn, nch, nch1, natom,
     &        nptp, nsh, nij0, nij_, nij, njc,
     &        l_ij0, k_ij0, l_ij_, k_ij_, l_ij, k_ij, l_jc, k_jc,
     &        nom1, nik1_0, nik1_, nik1, nk1c, l_om1, k_om1, l_ik1_0, 
     &        k_ik1_0, l_ik1_, k_ik1_, l_ik1, k_ik1, l_k1c, k_k1c,
     &        nom2, nik2_0, nik2_, nik2, nk2c, l_om2, k_om2, l_ik2_0, 
     &        k_ik2_0, l_ik2_, k_ik2_, l_ik2, k_ik2, l_k2c, k_k2c)
         call wfat_omega_wfat1(geom, ozn, nch, nqn, erg, chpos, nid,
     &        mid, lmax_, natomp, nptp, int_mb(k_icp), int_mb(k_ipp),
     &        gq, nlgl, zi, wfat1_exc, dcpl_mb(k_om1), dcpl_mb(k_om2))
      endif
      
      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Computation of WFAT '//
     &        'integrals', 48, .true.)
         write(luout, '("")')
         write(luout, 901)
         write(luout, '("")')
      endif
 901  format(
     &     5x, 'The calculation of WFAT integrals has started. ',
     &     'If this is a ''start'' type job, this part can take ', /,
     &     5x, 'up to several minutes, so please wait. The progress ',
     &     'of this calculation will be shown below.')

      
c     ==== BEGIN QUADRATURE SUMMATION ====
      call ga_sync()
      call wfat_get_intgl_chkpt(nptp, pr_chk)
      call wfat_init_progress(nptp, g_nptx, npt_tot)
      if (ndid == 0) write(luout, '(5x, "The total number of '//
     &     'quadrature grid points is ", i9, ".")') int(npt_tot)

      t_last = util_cpusec()
      t_start = t_last
      jpt = 0
      pr_count = 0
      kk = 1
      dcpl_mb(k_i_ : k_i_+ni_-1) = czero
      dcpl_mb(k_inc : k_inc+ninc-1) = czero
      if (wfat1) then
         dcpl_mb(k_ij_ : k_ij_+nij_-1) = czero
         dcpl_mb(k_ik1_ : k_ik1_+nik1_-1) = czero
         dcpl_mb(k_ik2_ : k_ik2_+nik2_-1) = czero

         dcpl_mb(k_jc : k_jc+njc-1) = czero
         dcpl_mb(k_k1c : k_k1c+nk1c-1) = czero
         dcpl_mb(k_k2c : k_k2c+nk2c-1) = czero
      endif

      
c     ==== Loop over atoms ====
      do ic = 1, natomp
         icent = int_mb(k_icp + ic - 1)
         if (icent == -1) then
            znuc = ozn   ! ozn = charge of the channel center.
         else
            if (.not. geom_cent_get(geom, icent, tag, atp, zn))
     &           call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
            znuc = int(zn)
         endif
         
c        ==== Get all grids associated to the current atom/center ====
         call wfat_get_gridpt(znuc, gq, max_npt, npt, dbl_mb(k_pt),
     &        dbl_mb(k_wt))

         ipt_i = int_mb(k_ipp + (ic-1)*2)
         ipt_f = int_mb(k_ipp + 1 + (ic-1)*2)
c        ==== Loop over all quadrature points around the current atom ====
         do ipt = ipt_i, ipt_f
            jpt = jpt + 1             ! Needed for wfat1.
c           jpt is always incremented regardless of whether the current
c           center is a gridless parabolic channel center. This is because
c           jpt is used to address k_om1 and k_om2, which were calculated
c           at all grids including the grids of a gridless channel center
c           (see wfat_omega_wfat1 above).

            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            if (icent == -1) then
               call wfat_mf_ogrid(nch, chpos, v_c, v_m)
            else
               call geom_cent2cent(geom, icent, v_c, imf, v_m, t1)
            endif
c           v_c = electron coordinate w.r.t. the current nuclei.
c           v_m = electron coordinate in the molecular frame, which is also the
c                 basis center.
            
c           ==== Calculate vc|psi> ====
            if (wfat1) then
               if ((icent == -1 .and. cgrid(1)) .or. icent /= -1) then
                  if (fedvr) then
                     call wfat_2a1efd_vcpsi_ord1(geom, sphi, mch, nfe, 
     &                    ndvr, maxn, xq, wq, nr, nt, nt_m, nt_p, 
     &                    pottype, nel, potparam, zi, vcf, cfx, cfy, 
     &                    cfz, chpos, v_m, dbl_mb(k_sc), vc_psi,
     &                    vc_cpsi)
                  else
                     call wfat_2a1e_vcpsi_ord1(geom, sphi, mch, nr, nt,
     &                    nr_m, nt_m, nr_p, nt_p, ke, pottype, nel, 
     &                    potparam, zi, vcf, cfx, cfy, cfz, chpos, v_m, 
     &                    dbl_mb(k_sc), vc_psi, vc_cpsi)
                  endif
               else
                  vc_psi = 0.0d0
                  vc_cpsi = 0.0d0
               endif
            else
               if ((icent == -1 .and. cgrid(1)) .or. icent /= -1) then
                  if (fedvr) then
                     call wfat_2a1efd_vcpsi(geom, sphi, mch, nfe, ndvr, 
     &                    maxn, xq, wq, nr, nt, pottype, nel, potparam,
     &                    zi, vcf, chpos, v_m, dbl_mb(k_sc), vc_psi)
                  else
                     call wfat_2a1e_vcpsi(geom, sphi, mch, nr, nt, ke, 
     &                    pottype, nel, potparam, zi, vcf, chpos, v_m,
     &                    dbl_mb(k_sc), vc_psi)
                  endif
               else
                  vc_psi = 0.0d0
                  vc_cpsi = 0.0d0
               endif
            endif
            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt_c(icent, geom, nch, v_m, erg, chpos, 
     &           cgrid, nqn, nid, mid, lmax_, int(zi), dcpl_mb(k_om))
            
c           ==== Calculate the integrand of I^{(n,m)}_{l,m'} ====
            call wfat_integrand(icent, nch, nqn, nsh, cgrid, 
     &           dcpl_mb(k_om), vc_psi, dcpl_mb(k_i0))
            if (isnan( abs(dcpl_mb(k_i0)) ))
     &           write(luout, '("isnan", i8, 3es14.6, 6x, 3es21.12E3)')
     &           ic, v_m, dcpl_mb(k_om), vc_psi
            
c           ==== Calculate the integrand of J^{(n,m)x_i}_{l,m'} ====
c           =====           and K^{(n,m)}_{l,m'}                ====
            if (wfat1) then
               i_om1 = k_om1 + (jpt-1)*nsh*nqn*nch
               i_om2 = k_om2 + (jpt-1)*nsh*nqn*nch
               call wfat_integrand_wfat1_c(icent, nch, nqn, nsh, cgrid, 
     &              dcpl_mb(k_om), dcpl_mb(i_om1), dcpl_mb(i_om2), 
     &              vc_cpsi, vc_psi, dcpl_mb(k_ij0), dcpl_mb(k_ik1_0),
     &              dcpl_mb(k_ik2_0))
            endif

c           ==== Apply masks to the integrands ====
            if (ismask) call wfat_mask_intd(geom, natom, nch, nqn, nsh,
     &           v_m, mask_param, dcpl_mb(k_i0))       ! omega(nsh,nqn)   nch=1
            if (ismask .and. wfat1) call wfat_mask_intdj(geom, natom, 
     &           nch, nqn, nsh, v_m, mask_param, dcpl_mb(k_ij0))     ! omega(nsh,nqn)   nch=1
            if (ismask .and. wfat1) then
               call wfat_mask_intd(geom, natom, nch, nqn, nsh, v_m,
     &              mask_param, dcpl_mb(k_ik1_0))           ! omega(nsh,nqn)   nch=1
               call wfat_mask_intd(geom, natom, nch, nqn, nsh, v_m,
     &              mask_param, dcpl_mb(k_ik2_0))           ! omega(nsh,nqn)   nch=1
            endif

c           ==== Quadrature summation ====
            call wfat_celf_1pt_c(geom, nch, geomc, icent, cgrid, chpos,
     &           v_m, celm, celf)
            wtcelf(:) = celf(:) * dbl_mb(k_wt+ipt-1)

c           ==== Quadrature summation ====
            if (icent /= -1) then
               i_l = k_i_ + (icent-1) * nsh*nqn*nch
               call wfat_integrand_acc(nsh, nqn, nch, wtcelf,
     &              dcpl_mb(k_i0), dcpl_mb(i_l))
            else
               call wfat_integrand_acc(nsh, nqn, nch, wtcelf,
     &              dcpl_mb(k_i0), dcpl_mb(k_inc))
            endif
            
            if (wfat1) then
               if (icent /= -1) then
                  i_l = k_ij_ + (icent-1) * nsh*nqn*nch1
                  call wfat_integrandj_acc(nsh, nqn, nch, wtcelf,
     &                 dcpl_mb(k_ij0), dcpl_mb(i_l))
               else
                  call wfat_integrandj_acc(nsh, nqn, nch, wtcelf,
     &                 dcpl_mb(k_ij0), dcpl_mb(k_jc))
               endif

               if (icent /= -1) then
                  i_l = k_ik1_ + (icent-1) * nsh*nqn*nch
                  call wfat_integrand_acc(nsh, nqn, nch, wtcelf,
     &                 dcpl_mb(k_ik1_0), dcpl_mb(i_l))
               else
                  call wfat_integrand_acc(nsh, nqn, nch, wtcelf,
     &                 dcpl_mb(k_ik1_0), dcpl_mb(k_k1c))
               endif

               if (icent /= -1) then
                  i_l = k_ik2_ + (icent-1) * nsh*nqn*nch
                  call wfat_integrand_acc(nsh, nqn, nch, wtcelf,
     &                 dcpl_mb(k_ik2_0), dcpl_mb(i_l))
               else
                  call wfat_integrand_acc(nsh, nqn, nch, wtcelf,
     &                 dcpl_mb(k_ik2_0), dcpl_mb(k_k2c))
               endif
            endif
            
cdebug            write(luout, '("check: ", i8, 10es15.6)')
cdebug     &           jpt, v_m, vc_psi,
cdebug     &           sum(abs(dcpl_mb(k_om1:k_om1+nom1-1))),
cdebug     &           sum(abs(dcpl_mb(k_ik1_0:k_ik1_0+nik1_0-1))),
cdebug     &           sum(abs(dcpl_mb(k_ik2_0:k_ik2_0+nik2_0-1))),
cdebug     &           sum(abs(dcpl_mb(k_ik1_:k_ik1_+nik1_-1))),
cdebug     &           sum(abs(dcpl_mb(k_ik2_:k_ik2_+nik2_-1))),
cdebug     &           wtcelf

            
c           ==== Print progress ====
            pr_count = pr_count + 1
            if (pr_count == pr_chk(kk)) then
               t_now = util_cpusec()
               dt_prog = t_now - t_last
               t_last = t_now
               call wfat_intgl_progress(kk, pr_count, dt_prog, npt_tot, 
     &              g_nptx, nptx)
               kk = kk + 1
            endif
            
         enddo
      enddo
      dt_all = t_now - t_start
      write(dt_str, '(f12.2)') dt_all
      if (ndid == 0) write(luout, '(5x, "The computation of the WFAT '//
     &     'integrals takes ", a, " s.")') trim(adjustl(dt_str))
      

c     ==== Distribute the integrals over grids belonging ====
c     ==== to the channel center for the current process ====
      call wfat_distrib_intgl(nsh, natom, nqn, nch, dcpl_mb(k_inc),
     &     dcpl_mb(k_i_))
      if (wfat1) then
         call wfat_distrib_intgl(nsh, natom, nqn, nch1, dcpl_mb(k_jc),
     &        dcpl_mb(k_ij_))
         call wfat_distrib_intgl(nsh, natom, nqn, nch, dcpl_mb(k_k1c),
     &        dcpl_mb(k_ik1_))
         call wfat_distrib_intgl(nsh, natom, nqn, nch, dcpl_mb(k_k2c),
     &        dcpl_mb(k_ik2_))
      endif

      
c     ==== Reorder elements of intgl ====
      call wfat_reorder_intgl(nsh, natom, nqn, nch,
     &     dcpl_mb(k_i_ : k_i_+ni_-1), dcpl_mb(k_il : k_il+nil-1))
      if (wfat1)
     &     call wfat_reorder_intgl_wfat1(nsh, natom, nqn, nch,
     &     dcpl_mb(k_ij_), dcpl_mb(k_ij),
     &     dcpl_mb(k_ik1_), dcpl_mb(k_ik1),
     &     dcpl_mb(k_ik2_), dcpl_mb(k_ik2))


c     ==== Accumulate to global intgl (g_il) ====
      call ga_zero(g_il)
      call ga_acc(g_il, 1, nil, 1, 1, dcpl_mb(k_il), nil, cone)
      call ga_sync()

      if (wfat1) then
         call ga_zero(g_j)
         call ga_acc(g_j, 1, nij, 1, 1, dcpl_mb(k_ij), nij, cone)
         call ga_sync()
         
         call ga_zero(g_k(1))
         call ga_acc(g_k(1), 1, nik1, 1, 1, dcpl_mb(k_ik1), nik1, cone)
         call ga_sync()
         
         call ga_zero(g_k(2))
         call ga_acc(g_k(2), 1, nik2, 1, 1, dcpl_mb(k_ik2), nik2, cone)
         call ga_sync()
      endif
      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_sc))
     &     call errquit('WFAT: Cannot deallocate scratch', 18, MA_ERR)
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      call wfat_cgeom_destroy(nch, geomc)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_dealloc_intgl(l_om, l_i0, l_i_, l_il, l_inc)
      if (wfat1)
     &     call wfat_dealloc_intgl_wfat1(l_ij0, l_ij_, l_ij, l_jc,
     &     l_om1, l_ik1_0, l_ik1_, l_ik1, l_k1c,
     &     l_om2, l_ik2_0, l_ik2_, l_ik2, l_k2c)
      call wfat_terminate_progress(g_nptx)
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_orbcent(rtdb, geom, orbmmt0, orbcent, orbmmt)

      implicit none
#include "geom.fh"
#include "errquit.fh"
      integer rtdb
      integer geom                 ! [input],
      double precision orbmmt0(3)  ! [input], Unshifted electronic dipole moment.
      double precision orbcent(3)  ! [output],
      double precision orbmmt(3)   ! [output]. Shifted electronic dipole moment.
      
      character*16 t
      logical umanpos(1)
      integer i, nc
      double precision r(3), q, totzn, manpos(3,1)


      if (.not. geom_ncent(geom, nc))
     &     call errquit('WFAT: Cannot get ncenters.', 40, GEOM_ERR)

c     Orbital center (orbcent) is given by D_{ion}/(tot. charge)
c     where D_{ion} and tot. charge are the total dipole and total
c     charge of cation. For 2a1e system (one electron initially), the
c     cation is a bare nuclei. Thus the above two quantities
c     correspond to the total dipole and charge of the nuclei.
      orbcent = 0.0d0
      do i = 1, nc
         if (.not. geom_cent_get(geom, i, t, r, q))
     &        call errquit('WFAT: geom center info failed', 88,
     &        GEOM_ERR)
         orbcent = orbcent + q * r
      enddo

      if (.not. geom_nuc_charge(geom, totzn))
     &     call errquit('WFAT: Getting totzn failed', 5, GEOM_ERR)
      call wfat_man_chpos(rtdb, 1, umanpos, manpos)
      if (umanpos(1)) then
         orbcent = manpos(:,1)
      else
         orbcent = orbcent / totzn
      endif
      orbmmt = orbmmt0 - orbcent
      
      end
c===============================================================
