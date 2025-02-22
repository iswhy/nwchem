c=======================================================
c=======================================================
c==== Imam ====      
      subroutine wfat_calc_intgl(nch, orbid, spinid, nqn, nid, mid,
     &                           lmax_, gqual, smh, erg, chmmt,
     &                           intgl)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "errquit.fh"
#include "geom.fh"
c     #include "cdft.fh"
#include "wfat.fh"
#include "wfat_qm_param.fh"      
      integer nch                     ! [input]
      integer orbid(nch)              ! [input]
      integer spinid(nch)             ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax_                   ! [input]
      integer gqual                   ! [input]
      double precision smh            ! [input]
      double precision erg(nch)       ! [output]
      double precision chmmt(3, nch)  ! [output]
      double complex intgl(1:*)       ! [output], dimension: (1+lmax_)**2 x natom x nqn x nch

      character*16 tag
      logical t1
      integer natom, nbas, ncf, nrd, max_npt, nvc, nsh, nom, ni0, ni_
      integer npt, ncpos, nscr
      integer i, j, icent, ipt, jpt, ich, iqn, ish, i_vc
      integer l_ps, l_cf, l_rd, l_pt, l_wt, l_vc, l_om
      integer l_i0, l_i_, l_sc
      integer k_ps, k_cf, k_rd, k_pt, k_wt, k_vc, k_om
      integer k_i0, k_i_, k_sc
      integer is, io, ipt0, i0l, i0r, oml, omr, i_l, i_r, imf
      double precision atp(3), znuc, v_c(3), v_m(3), celf
      double complex czero
      
      parameter (imf = 0)
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )



c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom_hd, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
c     ==== Get number of bases ====
      if (.not. bas_numbf(bas_hd, nbas))
     &     call errquit('WFAT: Getting nbas failed in ', 53, UERR)

      
      !!==== Obtain properties of each channel ====!!
      ncpos = 3 * nch
      if (.not. MA_alloc_get(MT_DBL, ncpos, 'chpos', l_ps, k_ps))
     &     call errquit('WFAT: chpos allocation failed', 59, MA_ERR)
      ncf = nbas * nch
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vch', l_cf, k_cf))
     &     call errquit('WFAT: vch allocation failed', 59, MA_ERR)

      write(luout, '("herep1")')
      call wfat_orb_property(nch, orbid, spinid, dbl_mb(k_cf), erg,
     &                       dbl_mb(k_ps), chmmt)
c     dbl_mb(k_ps : k_ps+ncpos-1) = 0.0d0
c      call wfat_get_channel_coef(nch, orbid, spinid, dbl_mb(k_cf))
      write(luout, '("sum vch = ", e14.5)')
     &      sum(dbl_mb(k_cf : k_cf+ncf-1))

c     ==== Get the reduced density matrix ====
      nrd = nbas * nbas * nspin0
      if (.not. MA_alloc_get(MT_DBL, nrd, 'rdm', l_rd, k_rd))
     &     call errquit('WFAT: rdm allocation failed', 59, MA_ERR)
      write(luout, '("herep2")')
      call wfat_get_rdm(dbl_mb(k_rd:k_rd+nrd-1))
      

c     ==== Allocate grid-related arrays ====
      call wfat_get_gridmaxpt(gqual, max_npt)
c     == Radial points ==
      if (.not. MA_alloc_get(MT_DBL, 3*max_npt, 'gpt', l_pt, k_pt))
     &     call errquit('WFAT: gpt allocation failed', 59, MA_ERR)
c     == Angular points ==
      if (.not. MA_alloc_get(MT_DBL, max_npt, 'gwt', l_wt, k_wt))
     &     call errquit('WFAT: gwt allocation failed', 59, MA_ERR)


c     ==== Allocate WFAT integral, its integrands, and some other ====
      nscr = 20*nbas + 4*nspin0 + nbas * (2*nbas+1+nspin0) + 4*nch
      if (.not. MA_alloc_get(MT_DBL, nscr, 'vscr', l_sc, k_sc))
     &     call errquit('WFAT: vscr allocation failed', 59, MA_ERR)
      nvc = nch
      if (.not. MA_alloc_get(MT_DBL, nvc, 'vcpsi', l_vc, k_vc))
     &     call errquit('WFAT: vcpsi allocation failed', 59, MA_ERR)
      nsh = (1 + lmax_)**2
      nom = nsh * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nom, 'omega', l_om, k_om))
     &     call errquit('WFAT: omega allocation failed', 59, MA_ERR)
      ni0 = nom
      if (.not. MA_alloc_get(MT_DCPL, ni0, 'intgl0', l_i0, k_i0))
     &     call errquit('WFAT: intgl0 allocation failed', 59, MA_ERR)
      ni_ = nsh * nqn * nch * natom
      if (.not. MA_alloc_get(MT_DCPL, ni_, 'intgl_', l_i_, k_i_))
     &     call errquit('WFAT: intgl_ allocation failed', 59, MA_ERR)      
      

c     ==== BEGIN QUADRATURE SUMMATION ====
c     ==== Loop over atoms ====
      jpt = 0
      dcpl_mb(k_i_ : k_i_+ni_-1) = czero
      do icent = 1, natom
         if (.not. geom_cent_get(geom_hd, icent, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)

         call wfat_get_gridpt(int(znuc), gqual, max_npt, npt,
     &                        dbl_mb(k_pt : k_pt+3*max_npt-1),
     &                        dbl_mb(k_wt : k_wt+max_npt-1))

c        ==== Loop over all quadrature points around the current atom ====
         do ipt = 1, npt
            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            call geom_cent2cent(geom_hd, icent, v_c, imf, v_m, t1)
cc            write(luout, '("        v_m = ", 3e14.5)') v_m

c           ==== Calculate vc|psi> ====
            call errquit('The input list of wfat_dft_vcpsi_1pt has '//
     &           'been modified', 52, UERR)
            call wfat_dft_vcpsi_1pt(geom_hd, bas_hd, nbas, nspin0, nch,
     &           dbl_mb(k_cf:k_cf+ncf-1), 
     &           dbl_mb(k_cf1:k_cf1+ncf-1), spinid,
     &           dbl_mb(k_ps:k_ps+ncpos-1), v_m, 
     &           dbl_mb(k_rd:k_rd+nrd-1), zion, smh,
     &           dbl_mb(k_sc:k_sc+nscr-1),
     &           dbl_mb(k_vc:k_vc+nvc-1))
            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt(geom_hd, v_m, nch, erg,
     &           dbl_mb(k_ps:k_ps+ncpos-1), nqn,
     &           nid, mid, lmax_, int(zion),
     &           dcpl_mb(k_om:k_om+nom-1))
            
c           ==== Calculate the integrand of I^(n,m)_(l,m') ====
            do i = 1, nch
               i0l = k_i0 + (i-1)*nsh*nqn
               i0r = i0l + nsh*nqn - 1
               oml = k_om + (i-1)*nsh*nqn
               omr = oml + nsh*nqn - 1
               i_vc = k_vc + i - 1
               if (abs(dbl_mb(i_vc)) > 1.0d-12) then
                  dcpl_mb(i0l:i0r) = conjg(dcpl_mb(oml:omr)) *
     &                               dbl_mb(i_vc)
               else
                  dcpl_mb(i0l:i0r) = 0.0d0
               endif
            enddo

c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom_hd, icent, v_c, wf_msm, celf)
            i_l = k_i_ + (icent-1) * nsh*nqn*nch
            i_r = i_l + nsh*nqn*nch - 1
            if (abs(celf) > 1.0d-12) then
               
               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                            dcpl_mb(k_i0:k_i0+ni0-1) *
     &                            celf *
     &                            dbl_mb(k_wt+ipt-1)     !! 1)
            endif

            jpt = jpt + 1
            write(luout, '("        sum vcpsi = ", e14.5)')
     &           sum(dbl_mb(k_vc:k_vc+nvc-1))
            write(luout, '("        sum intgl0 = ", 2e14.5)')
     &           sum(dcpl_mb(k_i0:k_i0+ni0-1))
            write(luout, '("        sum intgl_ = ", 2e14.5)')
     &           sum(dcpl_mb(i_l:i_r))
            write(luout, '(i7, " sum omega = ", 2e14.5)')
     &           jpt, sum(dcpl_mb(k_om:k_om+nom-1))
            if (abs(imag(sum(dcpl_mb(k_i0:k_i0+ni0-1)))) > 1.0d-8) then
c               call errquit('WFAT: Big imaginary part.', 2, UERR)
            endif
            
         enddo
      enddo
      
      
c     ==== Reorder elements of wf_intgl ====
      write(luout, '("here5")')
      call wfat_reorder_intgl(nsh, natom, nqn, nch,
     &     dcpl_mb(k_i_ : k_i_+ni_-1), intgl)

      j = 0
      do ich = 1, nch
         do iqn = 1, nqn
            do icent = 1, natom
               do ish = 1, nsh
                  j = j + 1
                  write(luout, '(i4, ")", 4i4, 2e14.5)')
     &                 j, ich, iqn, icent, ish, intgl(j)
               enddo
            enddo
         enddo
      enddo
      
                  
      
      write(luout, '("sum CHAN_INT = ", 2e14.5)') sum(intgl(1:ni_))
      write(luout, '("max of real part = ", e14.5)')
     &      maxval(abs(real(intgl(1:ni_))))
      write(luout, '("max of imag part = ", e14.5)')
     &      maxval(abs(imag(intgl(1:ni_))))


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_ps))
     &     call errquit('WFAT: chpos deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: vch deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_rd))
     &     call errquit('WFAT: rdm deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_pt))
     &     call errquit('WFAT: Grid points deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_wt))
     &     call errquit('WFAT: Grid weights deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc|psi> deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_om))
     &     call errquit('WFAT: omega deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_i0))
     &     call errquit('WFAT: intgl0 deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_i_))
     &     call errquit('WFAT: intgl_ deallocation failed',
     &     61, MA_ERR)
      
      end subroutine wfat_calc_intgl
!! 1) The quadrature summation is performed simultaneously for all channels (orbitals),
!!    each of which is given by 3rd index of wf_intgl. It's possible to do this because
!!    the cell function, the quadrature points and weights are the same for all channels.
!!    Note that it may not be possible to perform this simultaneous summation if the
!!    channels may or may not be assigned its own quadrature points. In the current
!!    implementation, all channels don't have quadrature points, so that the quadrature
!!    points (and weights) solely come from the atoms.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_reorder_intgl(nsh, ncent, nqn, nch, intgl_, intgl)

      implicit none
      integer nsh
      integer ncent
      integer nqn
      integer nch
      double complex intgl_(nsh, nqn, nch, ncent)
      double complex intgl(nsh, ncent, nqn, nch)

      integer ich, iqn, icent, i_lm
      
      
      do ich = 1, nch
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  intgl(i_lm, icent, iqn, ich) =
     &                     intgl_(i_lm, iqn, ich, icent)
               enddo
            enddo
         enddo
      enddo

      end
c=======================================================      
      
      

      
c!!=======================================================
c!!=======================================================
csubroutine wfat_check_big_intgl(icent, x_c, y_c, z_c, ii, chpos, &
c                                x_ch, y_ch, z_ch, iorb, &
c                                veffpsi, vspsi, vcpsi, max_omega, cellf, max_ingtl)
c
cc  use, intrinsic :: iso_c_binding
c  implicit none
c  integer icent, ii, iorb
c  double precision x_c, y_c, z_c, chpos(3), x_ch, y_ch, z_ch, &
c                                veffpsi, vspsi, vcpsi, max_omega, cellf, max_ingtl
c
c  
c  print *, " "
c  write(*,'("Checking for super big values:")')
c  write(*,'("  icent = ", i3)') icent
c  write(*,'("  Point and index = ", 3e15.6, i8)') x_c, y_c, z_c, ii
c  write(*,'("  Orbital center = ", 3e15.6)') chpos
c  write(*,'("  Distance from orb. center = ", e15.6)') sqrt(x_ch**2 + y_ch**2 + z_ch**2)
c  write(*,'("  orbital no. = ", i3)')  iorb
c  write(*,'("  veffpsi = ", e15.6)') veffpsi
c  write(*,'("  vspsi = ", e15.6)')   vspsi
c  write(*,'("  vcpsi = ", e15.6)')   vcpsi
c  write(*,'("  omega = ", e15.6)')  max_omega
c  write(*,'("  cell fun. = ", e15.6)') cellf
c  write(*,'("  max. abs. wf_intgl0 = ", e15.6)') max_ingtl
c  stop
c  
cend subroutine wfat_check_big_intgl
c!!=======================================================
