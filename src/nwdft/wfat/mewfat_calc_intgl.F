c=================================================
c=================================================
      subroutine mewfat_calc_intgl(rtdb, geom, bas, nb, gq, sx, sy, n0, 
     &     n1, ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, chpos, ionerg, 
     &     mcel, nqn, nid, mid, lmax_, smh, ismask, mask_param, wphase,
     &     rmode, g_il)

      implicit none
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input],
      integer geom
      integer bas
      integer nb
      integer gq               ! [input],
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)             ! [input], no. of electrons in cation.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients returned by mewfat_calc_detcoef. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient returned by mewfat_calc_detcoef.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient returned by mewfat_calc_detcoef. The 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients returned by mewfat_calc_detcoef. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input], charge of the parent ion.
      double precision chpos(3) ! [input],
      double precision ionerg   ! [input],
      integer mcel                                    ! [input], the m parameter related to the cell function. Usually taken as wf_msm from wfat.fh.
      integer nqn               ! [input],
      integer nid(nqn)            ! [input],
      integer mid(nqn)            ! [input],
      integer lmax_               ! [input],
      double precision smh      ! [input],
      logical ismask                                  ! [input], 
      double precision mask_param(2)                  ! [input]. Only used when ismask=.true.. Must always be supplied regardless.
      double complex wphase
      logical rmode             ! [input],
      integer g_il                ! [output],

      logical cgrid
      integer nat, ndid, natp, nptp, max_npt, nsh, nom, ni0, ni_, nil
      integer ninc, ozn, geomc
      integer l_icp, l_ipp, l_pt, l_wt, l_om, l_i0, l_i_, l_il, l_inc
      integer k_icp, k_ipp, k_pt, k_wt, k_om, k_i0, k_i_, k_il, k_inc
      double precision normch
      double complex cone

      parameter ( cone = cmplx(1.0d0, 0.0d0, kind=kind(1.0d0)) )


      call ga_sync()
      ndid = ga_nodeid()
      
c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
      

c     ==== Obtaining grid ranges for the current process ====
      if (.not. rtdb_get(rtdb, 'wfat:ozn', mt_int, 1, ozn))
     &     call errquit('WFAT: Cannot read ozn.', 73, RTDB_ERR)
      call wfat_alloc_gridrange(nat+1, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange_c(geom, ozn, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natp, nptp)
      call wfat_need_cgrid(rtdb, geom, 1, chpos, cgrid)
      call wfat_cgeom_create(geom, ozn, 1, chpos, geomc)
      call wfat_print_gridpt(rtdb, geom, gq, ozn, 1, cgrid, chpos)


c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gq, max_npt, l_pt, k_pt, l_wt, k_wt)
      
      
c     ==== Perform quadrature summation ====
      call wfat_alloc_intgl(lmax_, nqn, 1, nat, nsh, nom, ni0, ni_,
     &     nil, ninc, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il,
     &     l_inc, k_inc)
      call mewfat_calc_intgl_qsum(rtdb, geom, bas, nb, natp, 
     &     int_mb(k_icp), int_mb(k_ipp), gq, sx, sy, n0, n1, ocmo0,
     &     ocmo1, pcf, rcf, qcf, scf, zi, chpos, ionerg, max_npt, 
     &     dbl_mb(k_pt), dbl_mb(k_wt), cgrid, geomc, ozn, mcel, nsh, 
     &     nqn, nid, mid, lmax_, nat, smh, ismask, mask_param, wphase,
     &     dcpl_mb(k_om), dcpl_mb(k_i_), dcpl_mb(k_inc), dcpl_mb(k_i0),
     &     dcpl_mb(k_il))


c     ==== Accumulate to global intgl (g_il) ====
      call ga_zero(g_il)
      call ga_acc(g_il, 1, nil, 1, 1, dcpl_mb(k_il), nil, cone)
      call ga_sync()
      
      call ga_norm1(g_il, normch)
cfordebug      if (ndid==0) write(luout, '("intgl norm = ", es14.6)') normch


c     ==== Deallocations ====
      call wfat_cgeom_destroy(1, geomc)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      call wfat_dealloc_intgl(l_om, l_i0, l_i_, l_il, l_inc)

      end
c=================================================
      

c=================================================
c=================================================
      subroutine mewfat_calc_intgl_qsum(rtdb, geom, bas, nb, natp, icp,
     &     ipp, gq, sx, sy, n0, n1, ocmo0, ocmo1, pcf, rcf, qcf, scf, 
     &     zi, chpos, ionerg, max_npt, qpt, qwt, cgrid, geomc, ozn, 
     &     mcel, nsh, nqn, nid, mid, lmax_, nat, smh, ismask, 
     &     mask_param, wphase, omega, intg_, cintg, intg0, intg)

      implicit none
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
      integer rtdb                                    ! [input], RTDB handle.
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of basis functions in bas.
      integer natp                                    ! [input], no. of centers for the current process returned by wfat_get_gridrange.
      integer icp(natp)                               ! [input], centers for the current process returned by wfat_get_gridrange.
      integer ipp(2,natp)                             ! [input], range of points of each centers for the current process returned by wfat_get_gridrange.
      integer gq                                      ! [input], grid quality.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients returned by mewfat_calc_detcoef. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient returned by mewfat_calc_detcoef.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient returned by mewfat_calc_detcoef. The 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients returned by mewfat_calc_detcoef. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input], charge of the parent ion.
      double precision chpos(3)                       ! [input], center position of the channel returned by wfat_get_gridrange.
      double precision ionerg                         ! [input], ionization energy of the channel returned by wfat_get_gridrange.
      integer max_npt                                 ! [input], max. no. of points for the given gq as returned by wfat_get_gridmaxpt.
      double precision qpt(3, max_npt)                ! [output/scratch], array Quadrature points
      double precision qwt(max_npt) ! [output/scratch], array Quadrature weights
      logical cgrid
      integer geomc
      integer ozn
      integer mcel                                    ! [input], the m parameter related to the cell function. Usually taken as wf_msm from wfat.fh.
      integer nsh                                     ! [input], no. of spherical harmonics.
      integer nqn                                     ! [input], no. of parabolic quantum number pairs.
      integer nid(nqn)                                ! [input], n parabolic quantum numbers.
      integer mid(nqn)                                ! [input], m parabolic quantum numbers.
      integer lmax_                                   ! [input], max. angular momentum quantum number.
      integer nat                                     ! [input], no. of atoms in the molecule.
      double precision smh                            ! [input], smoothing parameter.
      logical ismask                                  ! [input], 
      double precision mask_param(2)                  ! [input], Only used when ismask=.true.. Must always be supplied regardless.
      double complex wphase                           ! [input], any wave function phase
      double complex omega(nsh, nqn)                  ! [output/scratch],
      double complex intg_(nsh, nqn, nat)             ! [output/scratch],
      double complex cintg(nsh, nqn)                  ! [output/scratch],
      double complex intg0(nsh, nqn)                  ! [output/scratch],
      double complex intg(nsh, nqn, nat)              ! [output], 

      character*5 si_type
      character*12 dt_str
      character*16 tag
      logical t1, rs, isdft
      integer ndid, nsc, l_sc, k_sc, l_rd0, k_rd0, npt, ic, icent
      integer i, imf, nsc_rs, nptp, pr_chk(10), pr_count, kk
      integer g_nptx, g_npt_tot, znuc
      double precision nptx, npt_tot, t_now, t_last, t_start, dt_all
      double precision atp(3), zn, v_m(3,1), vcpsi, celf(1), dt_prog
      double complex czero

      parameter (imf = 0)
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )

      logical nwxc_has_cam, nwxc_has_xc
      external nwxc_has_cam, nwxc_has_xc


      ndid = ga_nodeid()
      call ga_sync()

      
c     ==== Calculate neutral RDM in case of DFT ====
      rs = nwxc_has_cam()
      isdft = nwxc_has_xc()
      if (isdft) then
         if (.not. ma_alloc_get(mt_dbl, nb*nb*2, 'rdm0', l_rd0, k_rd0))
     &        call errquit('MEWFAT: Cannot allocate rdm0', 1, MA_ERR)
         call mewfat_calc_rdm0(nb, n0, ocmo0, dbl_mb(k_rd0))

         if (.not. rtdb_cget(rtdb, 'mewfat:si', 1, si_type))
     &        call errquit('MEWFAT: Failed reading SI from RTDB.', 1,
     &        RTDB_ERR)
      endif
      if (.not. isdft .and. rs)
     &     call errquit('MEWFAT: Range-separation can only be invoked'//
     &     ' when DFT is used.', 16, UERR)

      
c     ==== Scratch array ====
      if (rs) then
         nsc_rs = nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)
      else
         nsc_rs = 0
      endif
      nsc = 21*nb + n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) +
     &      sum(n1*n0) + nsc_rs                ! 1)
      if (.not. ma_alloc_get(mt_dbl, nsc, 'scratch', l_sc, k_sc))
     &     call errquit('MEWFAT: Cannot allocate scratch', 1, MA_ERR)
            

      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Computation of WFAT '//
     &        'integrals', 50, .true.)
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
c     ==== Loop over atoms ====
      nptp = 0
      do ic = 1, natp
         nptp = nptp + ipp(2,ic) - ipp(1,ic) + 1
      enddo
      call wfat_get_intgl_chkpt(nptp, pr_chk)
      call wfat_init_progress(nptp, g_nptx, npt_tot)
      if (ndid == 0) write(luout, '(5x, "The total number of '//
     &     'quadrature grid points is ", i9, ".")') int(npt_tot)

      t_last = util_cpusec()
      t_start = t_last
      pr_count = 0
      kk = 1
      intg_ = czero             ! intg_(nsh,nqn,nat)     nch=1
      cintg = czero
      do ic = 1, natp
         icent = icp(ic)        ! icp(natp)
         if (icent == -1) then
            znuc = ozn   ! ozn = charge of the channel center.
         else
            if (.not. geom_cent_get(geom, icent, tag, atp, zn))
     &           call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
            znuc = int(zn)
         endif
         
c        ==== Get all grids associated to the current atom/center ====
         call wfat_get_gridpt(znuc, gq, max_npt, npt, qpt, qwt)   ! qpt(3,max_npt)   qwt(max_npt)
c        The grid will always be calculated regardless of whether the current
c        center is a gridless parabolic channel center. This is because its
c        omission will be taken care of inside the subsequent subroutines.

         
c        ==== Loop over all quadrature points around the current atom ====
         do i = ipp(1,ic), ipp(2,ic)           ! ipp(2,natp)
            if (icent == -1) then
               call wfat_mf_ogrid(1, chpos, qpt(:,i), v_m)
            else
               call geom_cent2cent(geom, icent, qpt(:,i), imf, v_m(:,1),
     &              t1)
cforlater               v_m(1,2:) = v_m(1,1)
cforlater               v_m(2,2:) = v_m(2,1)
cforlater               v_m(3,2:) = v_m(3,1)
            endif

c           ==== Calculate vc|psi> ====
            if (isdft) then
               call mewfat_dft_vcpsi_1pt_c(icent, geom, bas, nb, sx, sy,
     &              n0, n1, v_m, cgrid, ocmo0, ocmo1, dbl_mb(k_rd0), 
     &              pcf, rcf, qcf, scf, zi, chpos, smh, rs, si_type, 
     &              dbl_mb(k_sc), vcpsi)
            else
               call mewfat_hf_vcpsi_1pt_c(icent, geom, bas, nb, sx, sy,
     &              n0, n1, v_m, cgrid, ocmo0, ocmo1, pcf, rcf, qcf, 
     &              scf, zi, chpos, smh, dbl_mb(k_sc), vcpsi)
            endif
            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt_c(icent, geom, 1, v_m, ionerg, chpos,
     &           cgrid, nqn, nid, mid, lmax_, int(zi), omega)  ! adjust v_m dimension      ! omega(nsh,nqn)   nch=1
            
c           ==== Calculate the integrand of I^(n,m)_(l,m') ====
            intg0 = conjg(omega) * vcpsi           ! intg0(nsh,nqn)   nch=1
            if (ismask) call wfat_mask_intd(geom, nat, 1, nqn, nsh, 
     &           v_m, mask_param, intg0) ! adjust v_m dimension      ! omega(nsh,nqn)   nch=1

c           ==== Cell function ====
            call wfat_celf_1pt_c(geom, 1, geomc, icent, cgrid, chpos,
     &           v_m, mcel, celf)
            
c           ==== Quadrature summation ====
            if (icent /= -1) then
               call wfat_integrand_acc(nsh, nqn, 1, celf*qwt(i),
     &              intg0, intg_(:,:,icent))
            else
               call wfat_integrand_acc(nsh, nqn, 1, celf*qwt(i),
     &              intg0, cintg)
            endif


cfordebug            if (ndid==0) then
cfordebug               write(luout, '("s1", i3, i8, 2es20.10, 3es20.10)')
cfordebug     &              ic, i, omega(1,1), vcpsi, celf, qwt(i)
cfordebug            endif

            
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
      call wfat_distrib_intgl(nsh, nat, nqn, 1, cintg, intg_)

      
c     ==== Reorder elements of intgl ====
      call wfat_reorder_intgl(nsh, nat, nqn, 1, intg_, intg)    ! intg(nsh,nqn,nat)    nch=1
      intg = wphase * intg
      call ga_sync()


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_sc))
     &     call errquit('MEWFAT: Cannot deallocate scratch', 1, MA_ERR)
      if (isdft) then
         if (.not. ma_free_heap(l_rd0))
     &        call errquit('MEWFAT: Cannot deallocate rdm0', 8, MA_ERR)
      endif
      call wfat_terminate_progress(g_nptx)
      
      end
c     NOTES:
c     1) This is the min. value for the length of scratch. This value is based on the 
c        info next to scr argument declaration in mewfat_potential.F.
c=================================================


c=================================================
c=================================================
      subroutine mewfat_alloc_ocmo(nb, n0, n1, nocmo0, nocmo1,
     &     l_mo0, k_mo0, l_mo1, k_mo1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nb              ! [input],
      integer n0(2)           ! [input],
      integer n1(2)           ! [input],
      integer nocmo0          ! [output],
      integer nocmo1          ! [output],
      integer l_mo0, k_mo0    ! [output],
      integer l_mo1, k_mo1    ! [output],

      
      nocmo0 = nb * sum(n0)
      if (.not. ma_alloc_get(mt_dbl, nocmo0, 'mo0', l_mo0, k_mo0))
     &     call errquit('MEWFAT: Cannot allocate mo0', 1, MA_ERR)
      nocmo1 = nb * sum(n1)
      if (.not. ma_alloc_get(mt_dbl, nocmo1, 'mo1', l_mo1, k_mo1))
     &     call errquit('MEWFAT: Cannot allocate mo1', 1, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_dealloc_ocmo(l_mo0, l_mo1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_mo0, l_mo1       ! [input]
      
      if (.not. MA_free_heap(l_mo0))
     &     call errquit('WFAT: mo0 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_mo1))
     &     call errquit('WFAT: mo1 deallocation failed', 61, MA_ERR)

      end
c=================================================
      

c=================================================
c=================================================
      subroutine mewfat_alloc_detcoef(sx, sy, n0, cx, cy, npcf, nrcf, 
     &     nqcf, nscf, l_pcf, k_pcf, l_rcf, k_rcf, l_qcf, k_qcf,
     &     l_scf, k_scf)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer sx              ! [input],
      integer sy              ! [input],
      integer n0(2)           ! [input], no. of electrons in the neutral.
      integer cx              ! [output],
      integer cy              ! [output],
      integer npcf            ! [output], 
      integer nrcf            ! [output], 
      integer nqcf            ! [output], 
      integer nscf            ! [output], 
      integer l_pcf, k_pcf    ! [output], 
      integer l_rcf, k_rcf    ! [output], 
      integer l_qcf, k_qcf    ! [output], 
      integer l_scf, k_scf    ! [output], 


      npcf = n0(sx)
      if (.not. ma_alloc_get(mt_dbl, npcf, 'pcf', l_pcf, k_pcf))
     &     call errquit('MEWFAT: Cannot allocate pcf', 1, MA_ERR)
      nrcf = 1
      if (.not. ma_alloc_get(mt_dbl, nrcf, 'rcf', l_rcf, k_rcf))
     &     call errquit('MEWFAT: Cannot allocate pcf', 1, MA_ERR)
      cx = max(1, n0(sx)-1)
      nqcf = n0(sx) * n0(sx) * cx
      if (.not. ma_alloc_get(mt_dbl, nqcf, 'qcf', l_qcf, k_qcf))
     &     call errquit('MEWFAT: Cannot allocate qcf', 1, MA_ERR)
      cy = max(1, n0(sy))
      nscf = cy * cy
      if (.not. ma_alloc_get(mt_dbl, nscf, 'scf', l_scf, k_scf))
     &     call errquit('MEWFAT: Cannot allocate pcf', 1, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_dealloc_detcoef(l_pcf, l_rcf, l_qcf, l_scf)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_pcf, l_rcf, l_qcf, l_scf


      if (.not. MA_free_heap(l_pcf))
     &     call errquit('WFAT: pcf deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_rcf))
     &     call errquit('WFAT: rcf deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_qcf))
     &     call errquit('WFAT: qcf deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_scf))
     &     call errquit('WFAT: scf deallocation failed', 61, MA_ERR)

      end
c=================================================
