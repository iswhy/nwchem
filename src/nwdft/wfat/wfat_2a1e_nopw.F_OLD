c===============================================================
c===============================================================
      subroutine wfat_2a1e_task_nopw(rtdb, geom, nqn, nid, mid, gq, zi,
     &     celm, sphi, mch, nr, nt, orbid, vcf, ke, erg, chmmt0, 
     &     pottype, nel, potparam, wfat1, nlgl, nr_m, nt_m, nr_p, nt_p,
     &     cfx, cfy, cfz)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "geom.fh"
      integer rtdb
      integer geom
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      integer gq
      double precision zi          ! [input],
      integer celm                 ! [input],
      integer sphi
      integer mch
      integer nr
      integer nt
      integer orbid
      double precision vcf(nr*nt)
      double precision ke
      double precision erg
      double precision chmmt0(3)
      character*3 pottype
      integer nel                ! not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)  ! 1)
      logical wfat1
      integer nlgl                 ! [input], only used when wfat1=.true.
      integer nr_m                 ! [input], only used when wfat1=.true.
      integer nt_m                 ! [input], only used when wfat1=.true.
      integer nr_p                 ! [input], only used when wfat1=.true.
      integer nt_p                 ! [input], only used when wfat1=.true.
      double precision cfx(1:*)    ! [input], only used when wfat1=.true.
      double precision cfy(1:*)    ! [input], only used when wfat1=.true.
      double precision cfz(nr*nt)  ! [input], only used when wfat1=.true.

      character*255 fname, fnamej, fnamek
      integer nat, g_il, g_j, g_k, un, unj, unk, nbt, ngm, nintg, nintj
      integer nintk, ig, ib, i, i_gn, i_jn, i_kn
      integer l_in, l_gn, l_inj, l_jn, l_ink, l_kn
      integer k_in, k_gn, k_inj, k_jn, k_ink, k_kn
      double precision chmmt(3), chpos(3), bt_start, bt_last, dbt
      double precision gm_start, gm_last, dgm, gm, bt
      parameter(un=85, unj=86, unk=87)
      

c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 94, GEOM_ERR)

      write(luout, '("nopw dipole = ", 3es16.8)') chmmt0
      
c     ==== Obtaining info on beta ====
      call wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)
c      dbtr = dbt * deg2rad


c     ==== Obtaining info on gamma ====
      call wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)
c      dgmr = dgm * deg2rad

      
      nintg = nat * nqn
      if (.not. ga_create(mt_dcpl, nintg, 1, 'intgl', 0, 0, g_il))
     &     call errquit('WFAT: Cannot create intgl global array.',
     &     77, GA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &     call errquit('WFAT: Cannot allocate intgl.', 24, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nqn, 'gn', l_gn, k_gn))
     &     call errquit('WFAT: Cannot allocate gn.', 24, MA_ERR)

      if (wfat1) then
         nintj = nat * nqn
         if (.not. ga_create(mt_dcpl, nintj, 1, 'intgl_j', 0, 0, g_j))
     &        call errquit('WFAT: Cannot create intgl_j global array.',
     &        77, GA_ERR)
         if (.not. MA_alloc_get(MT_DCPL, nintj, 'intgl_j',l_inj, k_inj))
     &        call errquit('WFAT: Cannot allocate intgl_j.', 24, MA_ERR)
         if (.not. MA_alloc_get(MT_DCPL, nqn, 'jn', l_jn, k_jn))
     &        call errquit('WFAT: Cannot allocate jn.', 24, MA_ERR)
         nintk = nat * nqn
         if (.not. ga_create(mt_dcpl, nintk, 1, 'intgl_k', 0, 0, g_k))
     &        call errquit('WFAT: Cannot create intgl_k global array.',
     &        77, GA_ERR)
         if (.not. MA_alloc_get(MT_DCPL, nintk, 'intgl_k',l_ink, k_ink))
     &        call errquit('WFAT: Cannot allocate intgl_k.', 24, MA_ERR)
         if (.not. MA_alloc_get(MT_DCPL, nqn, 'kn', l_kn, k_kn))
     &        call errquit('WFAT: Cannot allocate kn.', 24, MA_ERR)
      endif
      

      call util_file_name('wfats.nopw', .false., .false., fname)
      open(un, file=fname, action='write')
      if (wfat1) then
         call util_file_name('wfatj.nopw', .false., .false., fnamej)
         open(unj, file=fnamej, action='write')
         call util_file_name('wfatk.nopw', .false., .false., fnamek)
         open(unk, file=fnamek, action='write')
      endif
         
      
c     ==== Calculate and print structure factor and ion yield ====
      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)
         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)
            write(un, '(2es12.4)', advance='no') bt, gm
            if (wfat1) then
               write(unj, '(2es12.4)', advance='no') bt, gm
               write(unk, '(2es12.4)', advance='no') bt, gm
            endif
            
            
            call wfat_2a1e_intgl_nopw(rtdb, geom, nqn, nid, mid, bt, gm,
     &           gq, zi, celm, sphi, mch, nr, nt, ke, orbid, vcf, erg, 
     &           chmmt0, pottype, nel, potparam, wfat1, nlgl, nr_m, 
     &           nt_m, nr_p, nt_p, cfx, cfy, cfz, chpos, chmmt, g_il,
     &           g_j, g_k)

            write(luout, '("g_il:")')
            call ga_print(g_il)
            call ga_get(g_il, 1, nintg, 1, 1, dcpl_mb(k_in), nintg)
            call wfat_2a1e_srfc_nopw(bt, gm, erg, chmmt, nat, nqn,
     &           dcpl_mb(k_in), dcpl_mb(k_gn))
            if (wfat1) then
               call ga_get(g_j, 1, nintj, 1, 1, dcpl_mb(k_inj), nintj)
               call wfat_2a1e_srfc_nopw(bt, gm, erg, chmmt, nat, nqn,
     &              dcpl_mb(k_inj), dcpl_mb(k_jn))
               write(luout, '("g_k:")')
               call ga_print(g_k)
               call ga_get(g_k, 1, nintk, 1, 1, dcpl_mb(k_ink), nintk)
               call wfat_2a1e_srfc_nopw(bt, gm, erg, chmmt, nat, nqn,
     &              dcpl_mb(k_ink), dcpl_mb(k_kn))
               write(luout, '("kk__", f10.5, 4es20.8)')
     &              bt, dcpl_mb(k_ink), dcpl_mb(k_ink+1)
            endif
            
            

            do i = 1, nqn
               i_gn = k_gn + i - 1
               write(un, '(5x, 3es20.8)', advance='no')
     &              real(dcpl_mb(i_gn)), imag(dcpl_mb(i_gn)),
     &              (abs(dcpl_mb(i_gn)))**2
            enddo
            if (wfat1) then
               do i = 1, nqn
                  i_jn = k_jn + i - 1
                  write(unj, '(5x, 3es20.8)', advance='no')
     &                 real(dcpl_mb(i_jn)), imag(dcpl_mb(i_jn)),
     &                 (abs(dcpl_mb(i_jn)))**2
                  i_kn = k_kn + i - 1
                  write(unk, '(5x, 3es20.8)', advance='no')
     &                 real(dcpl_mb(i_kn)), imag(dcpl_mb(i_kn)),
     &                 (abs(dcpl_mb(i_kn)))**2
               enddo
            endif
            
            write(un, '("")')
            if (wfat1) then
               write(unj, '("")')
               write(unk, '("")')
            endif
            
         enddo
         write(un, '("")')
         if (wfat1) then
            write(unj, '("")')
            write(unk, '("")')
         endif
      enddo
      close(un)
      close(unj)
      close(unk)


      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 72, GA_ERR)
      if (.not. MA_free_heap(l_in))
     &     call errquit('WFAT: intgl deallocation failed', 24, MA_ERR)
      if (.not. MA_free_heap(l_gn))
     &     call errquit('WFAT: gn deallocation failed', 24, MA_ERR)

      if (wfat1) then
         if (.not. ga_destroy(g_j))
     &        call errquit('WFAT: Cannot destroy intgl_j.', 72, GA_ERR)
         if (.not. MA_free_heap(l_inj))
     &        call errquit('WFAT: intgl_j deallocation failed',
     &        24, MA_ERR)
         if (.not. MA_free_heap(l_jn))
     &        call errquit('WFAT: jn deallocation failed', 24, MA_ERR)

         if (.not. ga_destroy(g_k))
     &        call errquit('WFAT: Cannot destroy intgl_k.', 72, GA_ERR)
         if (.not. MA_free_heap(l_ink))
     &        call errquit('WFAT: intgl_k deallocation failed',
     &        24, MA_ERR)
         if (.not. MA_free_heap(l_kn))
     &        call errquit('WFAT: kn deallocation failed', 24, MA_ERR)
      endif
      
      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_intgl_nopw(rtdb, geom, nqn, nid, mid, bt, gm,
     &     gq, zi, celm, sphi, mch, nr, nt, ke, orbid, vcf, erg, chmmt0,
     &     pottype, nel, potparam, wfat1, nlgl, nr_m, nt_m, nr_p, nt_p, 
     &     cfx, cfy, cfz, chpos, chmmt, g_il, g_j, g_k)

      implicit none
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "geom.fh"
      integer rtdb                 ! [input],
      integer geom                 ! [input],
      integer nqn                  ! [input],
      integer nid(nqn)             ! [input],
      integer mid(nqn)             ! [input],
      double precision bt          ! [input], beta orientation angle in degrees.
      double precision gm          ! [input], gamma orientation angle in degrees.
      integer gq                   ! [input],
      double precision zi          ! [input],
      integer celm                 ! [input],
      integer sphi                 ! [input],
      integer mch                  ! [input],
      integer nr                   ! [input],
      integer nt                   ! [input],
      double precision ke          ! [input],
      integer orbid                ! [input],
      double precision vcf(nr*nt)  ! [input],
      double precision erg         ! [input],
      double precision chmmt0(3)   ! [input],
      character*3 pottype
      integer nel                  ! not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2) ! 1)
      logical wfat1                ! [input],
      integer nlgl                 ! [input], only used when wfat1=.true.
      integer nr_m                 ! [input], only used when wfat1=.true.
      integer nt_m                 ! [input], only used when wfat1=.true.
      integer nr_p                 ! [input], only used when wfat1=.true.
      integer nt_p                 ! [input], only used when wfat1=.true.
      double precision cfx(1:*)    ! [input], only used when wfat1=.true.
      double precision cfy(1:*)    ! [input], only used when wfat1=.true.
      double precision cfz(nr*nt)  ! [input], only used when wfat1=.true.
      double precision chpos(3)    ! [output],
      double precision chmmt(3)    ! [output],
      integer g_il                 ! [output],
      integer g_j                  ! [output], only calculated when wfat1=.true.
      integer g_k               ! [output], only calculated when wfat1=.true.

      character*12 dt_str
      character*16 tag
      logical ismask_g, ismask_j, ismask_k, t1
      integer ndid, natom, max_npt, nscr, ni0, ni_, nil
      integer nij0, nij_, nij, nik0, nik_, nik
      integer natomp, nptp, pr_chk(10), g_nptx, nptx, jpt, pr_count, kk
      integer ic, i_tu0, i_tu1
      integer icent, npt, ipt_i, ipt_f, ipt, ipt0, i_l, i_r
      integer l_pt, l_wt, l_sc, l_i0, l_i_, l_il, l_icp, l_ipp, l_tu0
      integer l_tu1, l_ij0, l_ij_, l_ij, l_ik0, l_ik_, l_ik, l_dr, l_dt
      integer k_pt, k_wt, k_sc, k_i0, k_i_, k_il, k_icp, k_ipp, k_tu0
      integer k_tu1, k_ij0, k_ij_, k_ij, k_ik0, k_ik_, k_ik, k_dr, k_dt
      double precision mask_g(2), mask_j(2), mask_k(2), npt_tot, t_now
      double precision t_last, t_start, celf, atp(3), znuc, v_c(3)
      double precision v_m(3), vc_psi, vc_cpsi(3), vc_cpsi_, wtcelf
      double precision dt_prog, dt_all, normch, iondip(3), protmat(3,3)
      double precision cu(3), sumvc, sumvcc
      double complex czero
      integer nch, imf, nsh
      parameter(nch=1, nsh=1, imf=0)
      parameter(czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)))
      

      call ga_sync()
      ndid = ga_nodeid()
      
      call ga_zero(g_il)
      if (wfat1) then
         call ga_zero(g_j)
         call ga_zero(g_k)
      endif


c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)

      
      call wfat_passrot(bt, gm, protmat)
      cu = protmat(:,3)
      

      ismask_g = .false.
      if (rtdb_get(rtdb, 'wfat:mask_g', mt_dbl, 2, mask_g))
     &     ismask_g = .true.
      ismask_j = .false.
      if (wfat1 .and. rtdb_get(rtdb, 'wfat:mask_j', mt_dbl, 2, mask_j))
     &     ismask_j = .true.
      ismask_k = .false.
      if (wfat1 .and. rtdb_get(rtdb, 'wfat:mask_k', mt_dbl, 2, mask_k))
     &     ismask_k = .true.

      if (ismask_g) write(luout, '("Masking g")')
      if (ismask_j) write(luout, '("Masking j")')
      if (ismask_k) write(luout, '("Masking k")')
      

c     ==== Obtain properties of the ionized channel ====
      call wfat_2a1e_orbcent(geom, chmmt0, chpos, chmmt)
      iondip = (/ 0.0d0, 0.0d0, 0.0d0 /)
      call wfat_print_orb_property(nch, orbid, 1, erg, chpos, chmmt,
     &     iondip)


c     ==== Some checks in wfat1 is used ====
      if (wfat1) then
clater         ovl1 = sum(vcf*cfz)
clater         if (.not. ma_alloc_get(MT_DBL, nr*nr, 'dr_dvr', l_dr, k_dr))
clater     &        call errquit('WFAT: Cannot allocate dr_dvr', 22, MA_ERR)
clater         if (.not. ma_alloc_get(MT_DBL, nt*nt, 'dt_dvr', l_dt, k_dt))
clater     &        call errquit('WFAT: Cannot allocate dt_dvr', 25, MA_ERR)
clater         call wfat_2a1e_dpo('z', mch, 0, nr, nt, nr, nt, sphi, sphi, ke,
clater     &        vcf, cfz, .false., dbl_mb(k_dr), .false., dbl_mb(k_dt), 
clater     &        .false., polb_)
clater         polb_ = -2 * polb_
clater         write(luout, '("<wfn0|wfn_z>, pol_zz = ", 2es21.10)')
clater     &        ovl1, polb_
clater         if (.not. ma_free_heap(l_dr))
clater     &        call errquit('WFAT: Cannot deallocate dr_dvr', 1, MA_ERR)
clater         if (.not. ma_free_heap(l_dt))
clater     &        call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)
      endif

      
c     ==== Obtaining grid ranges for the current process ====
      call wfat_alloc_quadgrid(gq, max_npt, l_pt, k_pt, l_wt, k_wt)
      call wfat_alloc_gridrange(natom, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natomp, nptp)

      
c     ==== Allocate WFAT integral, its integrands, and some other ====
      if (wfat1) then
         nscr = 
     &        nr +
     &        nt +
     &        nr*nt +
     &        nr_m +
     &        nt_m +
     &        nr_m*nt_m +
     &        nr_p +
     &        nt_p +
     &        nr_p*nt_p
      else
         nscr = 
     &        nr +
     &        nt +
     &        nr*nt
      endif
      write(luout, '("nscr = ", 3i5)') nscr, nr, nt
      if (.not. MA_alloc_get(MT_DBL, nscr, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 59, MA_ERR)
      call wfat_alloc_intgl_nopw(nqn, nch, natom, ni0, ni_, nil, l_i0,
     &     k_i0, l_i_, k_i_, l_il, k_il)
      if (.not. MA_alloc_get(MT_DBL, nqn*nptp, 'tau0', l_tu0, k_tu0))
     &     call errquit('WFAT: tau0 allocation failed', 59, MA_ERR)
      call wfat_collect_tau0_OLD(geom, nch, nqn, erg, chpos, nid, mid, 
     &     bt, gm, natomp, nptp, int_mb(k_icp), int_mb(k_ipp), gq, zi,
     &     dcpl_mb(k_tu0))

      if (wfat1) then
         call wfat_alloc_intgl_wfat1_nopw(nqn, nch, natom,
     &        nij0, nij_, nij, l_ij0, k_ij0, l_ij_, k_ij_, l_ij, k_ij,
     &        nik0, nik_, nik, l_ik0, k_ik0, l_ik_, k_ik_, l_ik, k_ik)
         if (.not. MA_alloc_get(MT_DBL, nqn*nptp, 'tau1', l_tu1, k_tu1))
     &        call errquit('WFAT: tau1 allocation failed', 59, MA_ERR)
         call wfat_collect_tau1_OLD(geom, nch, nqn, erg, chpos, nid, 
     &        mid, bt, gm, natomp, nptp, int_mb(k_icp), int_mb(k_ipp), 
     &        gq, zi, chmmt, dcpl_mb(k_tu1))
      endif

      

      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Computation the WFAT '//
     &        'integrals', 30, .true.)
         write(luout, '("")')
         write(luout, 901)
         write(luout, '("")')
      endif
 901  format(
     &     5x, 'The calculation of the WFAT integrals has started. ',
     &     'If this is a ''start'' type job, this part can take ', /,
     &     5x, 'up to several minutes, so please wait. The progress ',
     &     'of this calculation will be shown below.')


      
c     ==== BEGIN QUADRATURE SUMMATION ====
c     ==== Loop over atoms ====
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
      if (wfat1) then
         dcpl_mb(k_ij_ : k_ij_+nij_-1) = czero
         dcpl_mb(k_ik_ : k_ik_+nik_-1) = czero
      endif

      sumvc = 0.0d0
      sumvcc = 0.0d0
      do ic = 1, natomp
         icent = int_mb(k_icp + ic - 1)
         if (.not. geom_cent_get(geom, icent, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)

         call wfat_get_gridpt(int(znuc), gq, max_npt, npt,
     &                        dbl_mb(k_pt : k_pt+3*max_npt-1),
     &                        dbl_mb(k_wt : k_wt+max_npt-1))

         ipt_i = int_mb(k_ipp + (ic-1)*2)
         ipt_f = int_mb(k_ipp + 1 + (ic-1)*2)
c        ==== Loop over all quadrature points around the current atom ====
         do ipt = ipt_i, ipt_f
            jpt = jpt + 1
            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            call geom_cent2cent(geom, icent, v_c, imf, v_m, t1)

            
c           ==== Calculate vc|psi> ====
            if (wfat1) then
               call wfat_2a1e_vcpsi_ord1(geom, sphi, mch, nr, nt, nr_m,
     &              nt_m, nr_p, nt_p, ke, pottype, nel, potparam, zi, 
     &              vcf, cfx, cfy, cfz, chpos, v_m, dbl_mb(k_sc),
     &              vc_psi, vc_cpsi)
               call wfat_2a1e_vcpsi(geom, sphi, mch, nr, nt, ke, 
     &              pottype, nel, potparam, zi, vcf, chpos, v_m,
     &              dbl_mb(k_sc), vc_psi)
               vc_cpsi_ = sum(cu * vc_cpsi)
c               vc_cpsi_ = vc_cpsi(3)
               sumvc = sumvc + vc_psi
               sumvcc = sumvcc + vc_cpsi_
            else
               call wfat_2a1e_vcpsi(geom, sphi, mch, nr, nt, ke, 
     &              pottype, nel, potparam, zi, vcf, chpos, v_m,
     &              dbl_mb(k_sc), vc_psi)
               sumvc = sumvc + vc_psi
            endif

            
c           ==== Calculate the integrand of I^{(n,m)}_{l,m'} ====
            i_tu0 = k_tu0 + (jpt-1)*nsh*nqn*nch
            call wfat_integrand(nch, nqn, nsh, dcpl_mb(i_tu0), vc_psi, 
     &           dcpl_mb(k_i0))
            if (isnan( abs(dcpl_mb(k_i0)) )) then
               write(luout, '("isnan", i8, 3es14.6, 6x, 3es21.12E3)')
     &              ic, v_m, dcpl_mb(k_tu0), vc_psi
            endif
            
            
c           ==== Calculate the integrand of J^{(n,m)x_i}_{l,m'} ====
c           =====           and K^{(n,m)}_{l,m'}                ====
            if (wfat1) then
               i_tu1 = k_tu1 + (jpt-1)*nqn*nch
               call wfat_integrand(nch, nqn, nsh, dcpl_mb(i_tu1), 
     &              vc_psi, dcpl_mb(k_ik0))
               call wfat_integrand(nch, nqn, nsh, dcpl_mb(i_tu0), 
     &              vc_cpsi_, dcpl_mb(k_ij0))
            endif

            
c           ==== Apply masks to the integrands ====
            if (ismask_g)
     &           call wfat_mask_intd(geom, natom, nch, nqn, nsh,
     &           v_m, mask_g, dcpl_mb(k_i0))
            if (ismask_j .and. wfat1)
     &           call wfat_mask_intd(geom, natom, nch, nqn, nsh, 
     &           v_m, mask_j, dcpl_mb(k_ij0))
            if (ismask_k .and. wfat1) 
     &           call wfat_mask_intd(geom, natom, nch, nqn, nsh,
     &           v_m, mask_k, dcpl_mb(k_ik0))


c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom, icent, v_c, celm, celf)
            wtcelf = celf * dbl_mb(k_wt+ipt-1)

            i_l = k_i_ + (icent-1) * nsh*nqn*nch
            i_r = i_l + nsh*nqn*nch - 1
            dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                         dcpl_mb(k_i0:k_i0+ni0-1) *
     &                         wtcelf                        !! 1)
            
            if (wfat1) then
               i_l = k_ij_ + (icent-1) * nsh*nqn*nch
               i_r = i_l + nsh*nqn*nch - 1
               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                            dcpl_mb(k_ij0:k_ij0+nij0-1) *
     &                            wtcelf

               i_l = k_ik_ + (icent-1) * nsh*nqn*nch
               i_r = i_l + nsh*nqn*nch - 1
               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                            dcpl_mb(k_ik0:k_ik0+nik0-1) *
     &                            wtcelf
            endif
            
c            jpt = jpt + 1
c            write(luout, '("check: ", i8, 9es15.6)')
c     &           jpt, v_m, vc_psi, dcpl_mb(k_i0), dcpl_mb(k_i_), wtcelf

            
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
      write(luout, '("sumvc = ", f10.4, 2es16.8)') bt, sumvc, sumvcc

      write(luout, '("intgl_ = ", f10.4, 4es16.8)')
     &     bt, dcpl_mb(k_i_), dcpl_mb(k_i_+1)
      if (wfat1) write(luout, '("intgl_j_ = ", f10.4, 4es16.8)')
     &     bt, dcpl_mb(k_ik_), dcpl_mb(k_ik_+1)
      
      
      
c     ==== Reorder elements of intgl ====
      call wfat_reorder_intgl(nsh, natom, nqn, nch,
     &     dcpl_mb(k_i_ : k_i_+ni_-1), dcpl_mb(k_il : k_il+nil-1))
      if (wfat1) then
         call wfat_reorder_intgl(nsh, natom, nqn, nch,
     &        dcpl_mb(k_ij_ : k_ij_+nij_-1), dcpl_mb(k_ij : k_ij+nij-1))
         call wfat_reorder_intgl(nsh, natom, nqn, nch,
     &        dcpl_mb(k_ik_ : k_ik_+nik_-1), dcpl_mb(k_ik : k_ik+nik-1))
      endif
      
      
c     ==== Accumulate to global intgl (g_il) ====
      call ga_acc(g_il, 1, nil, 1, 1, dcpl_mb(k_il), nil, 1.0d0)
      call ga_sync()
      call ga_norm1(g_il, normch)      

      if (wfat1) then
         call ga_acc(g_j, 1, nij, 1, 1, dcpl_mb(k_ij), nij, 1.0d0)
c         call ga_zero(g_j)
         call ga_sync()
         call ga_norm1(g_j, normch)
         if (ndid==0) write(luout, '("intgl_j norm = ", es14.6)') normch

         call ga_acc(g_k, 1, nik, 1, 1, dcpl_mb(k_ik), nik, 1.0d0)
c         call ga_zero(g_k)
         call ga_sync()
         call ga_norm1(g_k, normch)
         if (ndid==0) write(luout, '("intgl_k norm = ", es14.6)') normch
      endif      


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: scratch deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_tu0))
     &     call errquit('WFAT: tau0 deallocation failed', 61, MA_ERR)
      if (wfat1) then
         if (.not. MA_free_heap(l_tu1))
     &        call errquit('WFAT: tau1 deallocation failed', 61, MA_ERR)
      endif
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_dealloc_intgl_nopw(l_i0, l_i_, l_il)
      if (wfat1)
     &     call wfat_dealloc_intgl_wfat1_nopw(l_ij0, l_ij_, l_ij,
     &     l_ik0, l_ik_, l_ik)
      call wfat_terminate_progress(g_nptx)

      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_alloc_intgl_nopw(nqn, nch, nat, ni0, ni_, nil,
     &     l_i0, k_i0, l_i_, k_i_, l_il, k_il)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nqn                       ! [input],
      integer nch                       ! [input],
      integer nat                       ! [input],
      integer ni0, ni_, nil             ! [output],
      integer l_i0, k_i0                ! [output],
      integer l_i_, k_i_                ! [output],
      integer l_il, k_il                ! [output],


      ni0 = nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, ni0, 'intgl0', l_i0, k_i0))
     &     call errquit('WFAT: intgl0 allocation failed', 59, MA_ERR)
      ni_ = nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, ni_, 'intgl_', l_i_, k_i_))
     &     call errquit('WFAT: intgl_ allocation failed', 59, MA_ERR)
      nil = nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nil, 'intgl', l_il, k_il))
     &     call errquit('WFAT: intgl allocation failed', 59, MA_ERR)

      end
c===============================================================
      

c=======================================================
c=======================================================
      subroutine wfat_dealloc_intgl_nopw(l_i0, l_i_, l_il)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_i0, l_i_, l_il     ! [input]
      
      if (.not. MA_free_heap(l_i0))
     &     call errquit('WFAT: intgl0 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_i_))
     &     call errquit('WFAT: intgl_ deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_il))
     &     call errquit('WFAT: intgl deallocation failed', 61, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_alloc_intgl_wfat1_nopw(nqn, nch, nat,
     &     nij0, nij_, nij, l_ij0, k_ij0, l_ij_, k_ij_, l_ij, k_ij,
     &     nik0, nik_, nik, l_ik0, k_ik0, l_ik_, k_ik_, l_ik, k_ik)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nqn                         ! [input],
      integer nch                         ! [input],
      integer nat                         ! [input], the total number of atoms in the system.
      integer nij0, nij_, nij        ! [output],
      integer l_ij0, k_ij0                ! [output],
      integer l_ij_, k_ij_                ! [output],
      integer l_ij, k_ij                  ! [output],
      integer nik0, nik_, nik   ! [output],
      integer l_ik0, k_ik0            ! [output],
      integer l_ik_, k_ik_              ! [output],
      integer l_ik, k_ik                ! [output],

      
c     ==== J ====
      nij0 = nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nij0, 'intd_j0', l_ij0, k_ij0))
     &     call errquit('WFAT: intd_j0 allocation failed', 59, MA_ERR)

      nij_ = nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nij_, 'intl_j_', l_ij_, k_ij_))
     &     call errquit('WFAT: intl_j_ allocation failed', 59, MA_ERR)

      nij = nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nij, 'intl_j', l_ij, k_ij))
     &     call errquit('WFAT: intl_j allocation failed', 59, MA_ERR)


c     ==== K ====
      nik0 = nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nik0, 'intd_k0', l_ik0, k_ik0))
     &     call errquit('WFAT: intd_k0 allocation failed', 59, MA_ERR)

      nik_ = nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nik_, 'intl_k_', l_ik_, k_ik_))
     &     call errquit('WFAT: intl_k_ allocation failed', 59, MA_ERR)

      nik = nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nik, 'intl_k', l_ik, k_ik))
     &     call errquit('WFAT: intl_k allocation failed', 59, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dealloc_intgl_wfat1_nopw(l_ij0, l_ij_, l_ij,
     &     l_ik0, l_ik_, l_ik)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_ij0, l_ij_, l_ij     ! [input]
      integer l_ik0, l_ik_, l_ik

      if (.not. MA_free_heap(l_ij0))
     &     call errquit('WFAT: intd_j0 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ij_))
     &     call errquit('WFAT: intl_j_ deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ij))
     &     call errquit('WFAT: intl_j deallocation failed', 61, MA_ERR)

      if (.not. MA_free_heap(l_ik0))
     &     call errquit('WFAT: intd_k0 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ik_))
     &     call errquit('WFAT: intl_k_ deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ik))
     &     call errquit('WFAT: intl_k deallocation failed', 61, MA_ERR)

      end
c=======================================================
            

c=======================================================
c=======================================================
      subroutine wfat_2a1e_srfc_nopw(bt, gm, erg, chmmt, nat, nqn,
     &     intgl, gn)

      implicit none
      double precision bt              ! [input], beta orientation angle in degrees.
      double precision gm              ! [input], gamma orientation angle in degrees.
      double precision erg
      double precision chmmt(3)
      integer nat
      integer nqn
      double complex intgl(nat, nqn)
      double complex gn(nqn)

      integer i
      double precision kp, mu(3), protmat(3,3)

      kp = sqrt(2 * abs(erg))
      call wfat_passrot(bt, gm, protmat)
      mu = matmul(transpose(protmat), -chmmt)

      do i = 1, nqn
         gn(i) = sum(intgl(:,i)) * exp(-kp*mu(3))
      enddo

      end
c=======================================================
      
