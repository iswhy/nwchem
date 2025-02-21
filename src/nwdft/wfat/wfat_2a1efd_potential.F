c===============================================================
c===============================================================
      subroutine wfat_2a1efd_vcpsi(geom, sphi, m, nfe, ndvr, maxn, xq, 
     &     wq, nr, nt, pottype, nel, potparam, zi, vcf, chpos, vpt, scr,
     &     vc_psi)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer geom
      integer sphi
      integer m
      integer nfe
      integer ndvr(nfe)
      integer maxn
      double precision xq(maxn,nfe)
      double precision wq(maxn,nfe)
      integer nr         ! nr=sum(ndvr)-nfe*2+nfe-1
      integer nt
      character*3 pottype
      integer nel                ! not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)  ! 1)
      double precision zi
      double precision vcf(nr*nt)
      double precision chpos(3)
      double precision vpt(3)
      double precision scr(maxn**2+maxn+nr+nt+nr*nt)
      double precision vc_psi

      character*1 sp
      integer ptm, pr, pt, pb, i
      double precision r, t, theta, deg2rad, pot, dr(3), s, sc
      double precision phi, wfn
      parameter(deg2rad=0.01745329251994d0)

      
c     ptm = pointer to temp array for wfat_eval0_fedvr
c     pr = pointer to radial DVR basis
c     pt = pointer to polar DVR basis
c     pb = pointer to the global basis
      ptm = 1
      pr = ptm + maxn**2 + maxn
      pt = pr + nr
      pb = pt + nt

      
c     The radial and polar DVR bases are defined in the molecular frame.
c     This is also the frame where the nuclei's positions are defined.
c     Since vpt contains the coordinate of the evaluation point in the
c     molecular frame, r and theta can directly be extracted from vpt
c     without prior shifting.

c     ==== Radial part ====
      r = sqrt(sum(vpt*vpt))
      call wfat_eval0_fedvr(nfe, ndvr, maxn, xq, wq, nr, r, scr(ptm),
     &     scr(pr))
      scr(pr:pr+nr-1) = scr(pr:pr+nr-1) / r
cc      write(luout, '("sum pr = ", 2es14.6)') r, sum(scr(pr:pr+nr-1))
cc      call wfat_dprintmat(1, nr, 'e', scr(pr))
cc      call errquit('stopping inside wfat_fd_pot.', 3, UERR)
      
c     ==== Theta part ====
      t = vpt(3) / r
      theta = acos(t) / deg2rad
      if (theta > 180.0d0 + 1.0e-12)  theta = 180.0d0
      if (theta <   0.0d0 - 1.0e-12)  theta =   0.0d0
      call wfat_eval_dvrt(m, nt, theta, scr(pt))

c     ==== The (r,\theta) part of the wave function ====
      call wfat_kron_prod(nr, 1, nt, 1, scr(pr), scr(pt), scr(pb))

c     ==== Phi part ====
      phi = atan2(vpt(2), vpt(1))
      if (sphi == 1) sp = 'c'
      if (sphi == 2) sp = 's'
      call wfat_2a1e_eval_phi(sp, m, phi, sc)
      

c     ==== The wavefunction evaluated at vpt ====
      wfn = 0.0d0
      do i = 1, nr*nt
         if (abs(vcf(i)) >= 1.0d-9) wfn = wfn + vcf(i) * scr(pb+i-1)
      enddo
      wfn = wfn * sc

      
      call wfat_2a1e_potential(geom, vpt, pottype, nel, potparam, pot)
      dr = vpt - chpos
      s = zi / sqrt(sum(dr*dr))
      vc_psi = -(pot + s) * wfn

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1efd_vcpsi_ord1(geom, sphi, m, nfe, ndvr, maxn,
     &     xq, wq, nr, nt, nt_m, nt_p, pottype, nel, potparam, zi, vcf, 
     &     cfx, cfy, cfz, chpos, vpt, scr, vc_psi, vc_cpsi)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer geom                    ! [input]
      integer sphi                    ! [input]
      integer m                       ! [input]
      integer nfe
      integer ndvr(nfe)
      integer maxn
      double precision xq(maxn,nfe)
      double precision wq(maxn,nfe)
      integer nr                      ! [input]
      integer nt                      ! [input]
      integer nt_m                    ! [input]
      integer nt_p                    ! [input]
      character*3 pottype
      integer nel                     ! not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)    ! [input], 1)
      double precision zi             ! [input]
      double precision vcf(nr*nt)     ! [input]
      double precision cfx(1:*)       ! [input]
      double precision cfy(1:*)       ! [input]
      double precision cfz(nr*nt)     ! [input]
      double precision chpos(3)       ! [input]
      double precision vpt(3)         ! [input]
      double precision scr(1:*)       ! [output],   1)
      double precision vc_psi         ! [output]
      double precision vc_cpsi(3)     ! [output]

      character*1 sp, sp_
      integer ptm, pr, pt, pb, pt_m, pb_m, pt_p, pb_p, n_m, n_p, i, ii
      double precision r, t, theta, deg2rad, pot, dr(3), s, sc
      double precision phi, scx_m, scx_p, scy_m, scy_p, wfnx_m, wfnx_p
      double precision wfn, wfny_m, wfny_p, wfnz
      parameter(deg2rad=0.01745329251994d0)

c     1) If m = 0, the minimum size of scr is (nr+nt+nr*nt) + (nt_p+nr*nt_p).
c        If m >= 1, the minimum size of scr is the above value for m = 0 plus 
c        nt_m+nr*nt_m.

      
      if (m >= 1) then
         n_m = nr * nt_m
      else if (m == 0) then
         n_m = 0
      endif
      n_p = nr * nt_p


c     pr = pointer to radial DVR basis of the m spectrum
c     pt = pointer to polar DVR basis of the m spectrum
c     pb = pointer to the global basis of the m spectrum
      ptm = 1
      pr = ptm + maxn**2 + maxn
      pt = pr + nr
      pb = pt + nt
      
c     pt_p = pointer to polar DVR basis of the m+1 spectrum
c     pb_p = pointer to the global basis of the m+1 spectrum
      pt_p = pb + nr*nt
      pb_p = pt_p + nt_p

c     pt_m = pointer to polar DVR basis of the m-1 spectrum
c     pb_m = pointer to the global basis of the m-1 spectrum
      if (m >= 1) then
         pt_m = pb_p + nr*nt_p
         pb_m = pt_m + nt_m
      else
         pt_m = 0
         pb_m = 0
      endif

      
c     The radial and polar DVR bases are defined in the molecular frame.
c     This is also the frame where the nuclei's positions are defined.
c     Since vpt contains the coordinate of the evaluation point in the
c     molecular frame, r and theta can directly be extracted from vpt
c     without prior shifting.

c     ==== Radial part ====
      r = sqrt(sum(vpt*vpt))
      call wfat_eval0_fedvr(nfe, ndvr, maxn, xq, wq, nr, r, scr(ptm),
     &     scr(pr))
c      write(luout, '("nfe_ = ", 3i6, 2es14.6, i4)')
c     &     nfe, sum(ndvr), maxn, sum(xq), sum(wq), nr
      scr(pr:pr+nr-1) = scr(pr:pr+nr-1) / r
c      write(luout, '("sum pr = ", 2es14.6)') r, sum(scr(pr:pr+nr-1))
c      call wfat_dprintmat(1, nr, 'e', scr(pr))
      
c     ==== Theta part ====
      t = vpt(3) / r
      theta = acos(t) / deg2rad
      if (theta > 180.0d0 + 1.0e-12)  theta = 180.0d0
      if (theta <   0.0d0 - 1.0e-12)  theta =   0.0d0
      call wfat_eval_dvrt(m, nt, theta, scr(pt))
      if (m >= 1)
     &     call wfat_eval_dvrt(m-1, nt_m, theta, scr(pt_m))
      call wfat_eval_dvrt(m+1, nt_p, theta, scr(pt_p))
      
c     ==== The (r,\theta) part of the wave function ====
      call wfat_kron_prod(nr, 1, nt, 1, scr(pr), scr(pt), scr(pb))
      if (m >= 1)
     &     call wfat_kron_prod(nr, 1, nt_m, 1, scr(pr), scr(pt_m),
     &     scr(pb_m))
      call wfat_kron_prod(nr, 1, nt_p, 1, scr(pr), scr(pt_p),
     &     scr(pb_p))
      
c     ==== Phi part ====
      phi = atan2(vpt(2), vpt(1))
      if (sphi == 1) then
         sp = 'c'
         sp_ = 's'
      else if (sphi == 2) then
         sp = 's'
         sp_ = 'c'
      endif
      call wfat_2a1e_eval_phi(sp, m, phi, sc)
      if (m >= 1) call wfat_2a1e_eval_phi(sp, m-1, phi, scx_m)
      call wfat_2a1e_eval_phi(sp, m+1, phi, scx_p)
      if (m >= 1) call wfat_2a1e_eval_phi(sp_, m-1, phi, scy_m)
      call wfat_2a1e_eval_phi(sp_, m+1, phi, scy_p)

c     ==== The potential term ====
      call wfat_2a1e_potential(geom, vpt, pottype, nel, potparam, pot)
      dr = vpt - chpos
      s = zi / sqrt(sum(dr*dr))

      
c     ==== The 0th order wavefunction evaluated at vpt ====
      wfn = 0.0d0
      do i = 1, nr*nt
         if (abs(vcf(i)) >= 1.0d-9) wfn = wfn + vcf(i) * scr(pb+i-1)
      enddo
      wfn = wfn * sc

      vc_psi = -(pot + s) * wfn
      if (isnan(abs(vc_psi))) write(luout, '("isnan vc", 7es18.7)')
     &     pot, s, wfn, sum(vcf), sum(scr(pr:pr+nr-1)),
     &     sum(scr(pt:pt+nt-1)), r
c      call errquit('stopping inside wfat_fd_pot.', 3, UERR)
      
      
c     ==== The m-1 part of the 1st order correction evaluated at vpt ====
      if (m >= 1) then
         wfnx_m = 0.0d0
         wfny_m = 0.0d0
         do i = 1, n_m
            if (abs(cfx(i)) >= 1.0d-9)
     &           wfnx_m = wfnx_m + cfx(i) * scr(pb_m+i-1)
            if (abs(cfy(i)) >= 1.0d-9)
     &           wfny_m = wfny_m + cfy(i) * scr(pb_m+i-1)
         enddo
         wfnx_m = wfnx_m * scx_m
         wfny_m = wfny_m * scy_m
      else if (m == 0) then
         wfnx_m = 0.0d0
         wfny_m = 0.0d0
      endif

c     ==== The m+1 part of the 1st order correction evaluated at vpt ====
      wfnx_p = 0.0d0
      wfny_p = 0.0d0
      do i = n_m+1, n_m+n_p
         ii = i - n_m
         if (abs(cfx(i)) >= 1.0d-9)
     &        wfnx_p = wfnx_p + cfx(i) * scr(pb_p+ii-1)
         if (abs(cfy(i)) >= 1.0d-9)
     &        wfny_p = wfny_p + cfy(i) * scr(pb_p+ii-1)
      enddo
      wfnx_p = wfnx_p * scx_p
      wfny_p = wfny_p * scy_p

c     ==== The z part of the 1st order correction evaluated at vpt ====
      wfnz = 0.0d0
      do i = 1, nr*nt
         if (abs(cfz(i)) >= 1.0d-9) wfnz = wfnz + cfz(i) * scr(pb+i-1)
      enddo
      wfnz = wfnz * sc

c     ==== V_c \psi^{(1)} ====
      vc_cpsi(1) = -(pot + s) * (wfnx_m + wfnx_p)
      vc_cpsi(2) = -(pot + s) * (wfny_m + wfny_p)
      vc_cpsi(3) = -(pot + s) * wfnz
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1efd_collect_vcpsi(geom, sphi, m, nfe, ndvr, 
     &        maxn, xq, wq, nr, nt, nptp, pottype, nel, potparam, zi, 
     &        vcf, chpos, vpt, vc_psi)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                  ! [input],
      integer sphi                  ! [input],
      integer m                     ! [input],
      integer nfe                   ! [input],
      integer ndvr(*)               ! [input], size is nfe
      integer maxn                  ! [input],
      double precision xq(*)        ! [input], size is maxn*nfe
      double precision wq(*)        ! [input], size is maxn*nfe
      integer nr                    ! [input],
      integer nt                    ! [input],
      integer nptp                  ! [input],
      character*3 pottype           ! [input],
      integer nel                   ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)  ! [input], 1)
      double precision zi           ! [input],
      double precision vcf(nr*nt)   ! [input],
      double precision chpos(3)     ! [input],
      double precision vpt(3,nptp)  ! [input],
      double precision vc_psi(nptp) ! [output],

      integer nscr, l_sc, k_sc, i

      nscr = nr + nt + nr*nt + maxn**2 + maxn
      if (.not. MA_alloc_get(MT_DBL, nscr, 'vscr', l_sc, k_sc))
     &     call errquit('WFAT: vscr allocation failed', 59, MA_ERR)
      
      do i = 1, nptp
         call wfat_2a1efd_vcpsi(geom, sphi, m, nfe, ndvr, maxn, xq, 
     &        wq, nr, nt, pottype, nel, potparam, zi, vcf, chpos, 
     &        vpt(:,i), dbl_mb(k_sc), vc_psi(i))
      enddo

      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)
      
      end
c===============================================================
      
      
c===============================================================
c===============================================================
      subroutine wfat_2a1efd_collect_vcpsi_ord1(geom, sphi, m, nfe, 
     &     ndvr, maxn, xq, wq, nr, nt, nt_m, nt_p, nptp, pottype, nel, 
     &     potparam, zi, vcf, cfx, cfy, cfz, chpos, vpt, vc_psi, 
     &     vc_cpsi)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                      ! [input],
      integer sphi                      ! [input],
      integer m                         ! [input],
      integer nfe                       ! [input],
      integer ndvr(*)                   ! [input], size is nfe
      integer maxn                      ! [input],
      double precision xq(*)            ! [input], size is maxn*nfe
      double precision wq(*)            ! [input], size is maxn*nfe
      integer nr                        ! [input],
      integer nt                        ! [input],
      integer nt_m                      ! [input],
      integer nt_p                      ! [input],
      integer nptp                      ! [input],
      character*3 pottype               ! [input],
      integer nel                       ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)      ! [input], 1)
      double precision zi               ! [input],
      double precision vcf(nr*nt)       ! [input],
      double precision cfx(1:*)         ! [input],
      double precision cfy(1:*)         ! [input],
      double precision cfz(nr*nt)       ! [input],
      double precision chpos(3)         ! [input],
      double precision vpt(3,nptp)      ! [input],
      double precision vc_psi(nptp)     ! [input],
      double precision vc_cpsi(nptp, 3) ! [output],

      integer nscr, l_sc, k_sc, i

      nscr = maxn**2 + maxn +
     &     (nr + nt + nr*nt) +
     &     (nt_p + nr*nt_p)
      if (m >= 1) nscr = nscr + (nt_m + nr*nt_m)
      if (.not. MA_alloc_get(MT_DBL, nscr, 'vscr', l_sc, k_sc))
     &     call errquit('WFAT: vscr allocation failed', 59, MA_ERR)
      
      do i = 1, nptp
         call wfat_2a1efd_vcpsi_ord1(geom, sphi, m, nfe, ndvr, maxn, xq,
     &        wq, nr, nt, nt_m, nt_p, pottype, nel, potparam, zi, vcf, 
     &        cfx, cfy, cfz, chpos, vpt(:,i), dbl_mb(k_sc), vc_psi(i), 
     &        vc_cpsi(i,:))
      enddo

      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)

      end
c===============================================================
      
