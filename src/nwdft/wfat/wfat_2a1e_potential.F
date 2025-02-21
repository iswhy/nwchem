c===============================================================
c===============================================================
      subroutine wfat_2a1e_potential(geom, vpt, pottype, nel, param,
     &     pot)

      implicit none
#include "errquit.fh"
#include "geom.fh"
      integer geom
      double precision vpt(3)
      character*3 pottype
      integer nel                ! not used when pottype = 'cou', but still has to be supplied a value.
      double precision param(2)  ! 1)
      double precision pot

c     1) When pottype = 'cou', param(1) = soft Coulomb parameter, while param(2) is unused.
c        When pottype = 'eff', param(1:2) = the two screening parameters.

      character*16 t
      integer i, nc
      double precision dr(3), atp(3), q, dist

      
      if (pottype /= 'cou' .and. pottype /= 'eff') then
         call errquit('WFAT: Undefined value assigned to the input '//
     &        'argument pottype of wfat_2a1e_potential (pottype = '//
     &        pottype//').')
      endif
      
      if (.not. geom_ncent(geom, nc))
     &     call errquit('WFAT: Cannot get ncenters.', 28, GEOM_ERR)

      pot = 0.0d0
      do i = 1, nc
         if (.not. geom_cent_get(geom, i, t, atp, q))
     &        call errquit('WFAT: geom center info failed.', 45,
     &        GEOM_ERR)
         dr = vpt - atp
         if (pottype == 'cou') then
            dist = sqrt(sum(dr*dr) + param(1))
            pot = pot + (-q)/dist
         else if (pottype == 'eff') then
            dist = sqrt(sum(dr*dr))
            pot = pot - (nel - (nel-1)*
     &           (1-1.0d0/(param(2)/param(1)*(exp(param(1)*dist)-1)+1)))
     &           / dist
         endif
      enddo

ccc        == Ne ==
ccc         u = 1.704d0
ccc         v = 2.810d0
cc
ccc        == Ar ==
cc         u = 0.933d0
cc         v = 3.600d0
cc
ccc        == Kr ==
cc         u = 1.340d0
cc         v = 4.311d0
cc         
      

      end
c===============================================================
      
      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_vcpsi(geom, sphi, m, nr, nt, ke, pottype, 
     &     nel, potparam, zi, vcf, chpos, vpt, scr, vc_psi)

      implicit none
#include "stdio.fh"
      integer geom
      integer sphi
      integer m
      integer nr
      integer nt
      double precision ke
      character*3 pottype
      integer nel                ! not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)  ! 1)
      double precision zi
      double precision vcf(nr*nt)
      double precision chpos(3)
      double precision vpt(3)
      double precision scr(nr+nt+nr*nt)
      double precision vc_psi

      character*1 sp
      integer pr, pt, pb
      double precision r, t, theta, deg2rad, pot, dr(3), s, a, b, kp, sc
      double precision phi
      parameter(deg2rad=0.01745329251994d0)


c     pr = pointer to radial DVR basis
c     pt = pointer to polar DVR basis
c     pb = pointer to the global basis
      pr = 1
      pt = 1 + nr
      pb = 1 + nr + nt

      
c     The radial and polar DVR bases are defined in the molecular frame.
c     This is also the frame where the nuclei's positions are defined.
c     Since vpt contains the coordinate of the evaluation point in the
c     molecular frame, r and theta can directly be extracted from vpt
c     without prior shifting.
      a = 1.0d0
      b = ke * 50.0d0
      kp = ke

c     ==== Radial part ====
      r = sqrt(sum(vpt*vpt))
      call wfat_eval_dvrr_mask(m, nr, ke, r, a, b, zi, kp, scr(pr))

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
      
      
      call wfat_2a1e_potential(geom, vpt, pottype, nel, potparam, pot)
      dr = vpt - chpos
      s = zi / sqrt(sum(dr*dr))
      vc_psi = -(pot + s) * sum(vcf(:) * scr(pb:pb+nr*nt-1)) * sc

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_vcpsi_ord1(geom, sphi, m, nr, nt, nr_m, nt_m,
     &     nr_p, nt_p, ke, pottype, nel, potparam, zi, vcf, cfx, cfy, 
     &     cfz, chpos, vpt, scr, vc_psi, vc_cpsi)

      implicit none
      integer geom                    ! [input]
      integer sphi                    ! [input]
      integer m                       ! [input]
      integer nr                      ! [input]
      integer nt                      ! [input]
      integer nr_m                    ! [input]
      integer nt_m                    ! [input]
      integer nr_p                    ! [input]
      integer nt_p                    ! [input]
      double precision ke             ! [input]
      character*3 pottype             ! [input]
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
      integer pr, pt, pb, pr_m, pt_m, pb_m, pr_p, pt_p, pb_p, n_m, n_p
      double precision r, t, theta, deg2rad, pot, dr(3), s, a, b, kp, sc
      double precision phi, scx_m, scx_p, scy_m, scy_p, wfnx_m, wfnx_p
      double precision wfny_m, wfny_p, wfnz
      parameter(deg2rad=0.01745329251994d0)

c     1) If m = 0, the minimum size of scr is nr + nt + nr*nt + nr_p + nt_p + nr_p*nt_p.
c        If m >= 1, the minimum size of scr is the above value for m = 0 plus 
c        nr_m + nt_m + nr_m*nt_m.

      
      if (m >= 1) then
         n_m = nr_m * nt_m
      else if (m == 0) then
         n_m = 0
      endif
      n_p = nr_p * nt_p


c     pr = pointer to radial DVR basis of the m spectrum
c     pt = pointer to polar DVR basis of the m spectrum
c     pb = pointer to the global basis of the m spectrum
      pr = 1
      pt = pr + nr
      pb = pt + nt
      
c     pr_p = pointer to radial DVR basis of the m+1 spectrum
c     pt_p = pointer to polar DVR basis of the m+1 spectrum
c     pb_p = pointer to the global basis of the m+1 spectrum
      pr_p = pb + nr*nt
      pt_p = pr_p + nr_p
      pb_p = pt_p + nt_p

c     pr_m = pointer to radial DVR basis of the m-1 spectrum
c     pt_m = pointer to polar DVR basis of the m-1 spectrum
c     pb_m = pointer to the global basis of the m-1 spectrum
      if (m >= 1) then
         pr_m = pb_p + nr_p*nt_p
         pt_m = pr_m + nr_m
         pb_m = pt_m + nt_m
      else
         pr_m = 0
         pt_m = 0
         pb_m = 0
      endif
      
      
c     The radial and polar DVR bases are defined in the molecular frame.
c     This is also the frame where the nuclei's positions are defined.
c     Since vpt contains the coordinate of the evaluation point in the
c     molecular frame, r and theta can directly be extracted from vpt
c     without prior shifting.
      a = 1.0d0
      b = ke * 50.0d0
      kp = ke

c     ==== Radial part ====
      r = sqrt(sum(vpt*vpt))
      call wfat_eval_dvrr_mask(m, nr, ke, r, a, b, zi, kp, scr(pr))
      if (m >= 1)
     &     call wfat_eval_dvrr_mask(m-1, nr_m, ke, r, a, b, zi, kp,
     &     scr(pr_m))
      call wfat_eval_dvrr_mask(m+1, nr_p, ke, r, a, b, zi, kp,
     &     scr(pr_p))
      

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
     &     call wfat_kron_prod(nr_m, 1, nt_m, 1, scr(pr_m), scr(pt_m),
     &     scr(pb_m))
      call wfat_kron_prod(nr_p, 1, nt_p, 1, scr(pr_p), scr(pt_p),
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
      vc_psi = -(pot + s) * sum(vcf(:) * scr(pb:pb+nr*nt-1)) * sc
      
c     ==== The m-1 part of the 1st order correction evaluated at vpt ====
      if (m >= 1) then
         wfnx_m = sum(cfx(1:n_m) * scr(pb_m:pb_m+n_m-1)) * scx_m
         wfny_m = sum(cfy(1:n_m) * scr(pb_m:pb_m+n_m-1)) * scy_m
      else if (m == 0) then
         wfnx_m = 0.0d0
         wfny_m = 0.0d0
      endif

c     ==== The m+1 part of the 1st order correction evaluated at vpt ====
      wfnx_p = sum(cfx(n_m+1:n_m+n_p) * scr(pb_p:pb_p+n_p-1)) * scx_p
      wfny_p = sum(cfy(n_m+1:n_m+n_p) * scr(pb_p:pb_p+n_p-1)) * scy_p

c     ==== The z part of the 1st order correction evaluated at vpt ====
      wfnz = sum(cfz(:) * scr(pb:pb+nr*nt-1)) * sc

c     ==== V_c \psi^{(1)} ====
      vc_cpsi(1) = -(pot + s) * (wfnx_m + wfnx_p)
      vc_cpsi(2) = -(pot + s) * (wfny_m + wfny_p)
      vc_cpsi(3) = -(pot + s) * wfnz
      
      end
c===============================================================
      

c===============================================================
c===============================================================
      subroutine wfat_2a1e_collect_vcpsi(geom, sphi, m, nr, nt, nptp, 
     &     ke, pottype, nel, potparam, zi, vcf, chpos, vpt, vc_psi)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                  ! [input],
      integer sphi                  ! [input],
      integer m                     ! [input],
      integer nr                    ! [input],
      integer nt                    ! [input],
      integer nptp                  ! [input],
      double precision ke           ! [input],
      character*3 pottype           ! [input],
      integer nel                   ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)  ! [input], 1)
      double precision zi           ! [input],
      double precision vcf(nr*nt)   ! [input],
      double precision chpos(3)     ! [input],
      double precision vpt(3,nptp)  ! [input],
      double precision vc_psi(nptp) ! [output],

      integer nscr, l_sc, k_sc, i

      nscr = 
     &        nr +
     &        nt +
     &        nr*nt
      if (.not. MA_alloc_get(MT_DBL, nscr, 'vscr', l_sc, k_sc))
     &     call errquit('WFAT: vscr allocation failed', 59, MA_ERR)
      
      do i = 1, nptp
         call wfat_2a1e_vcpsi(geom, sphi, m, nr, nt, ke, pottype, nel, 
     &        potparam, zi, vcf, chpos, vpt(:,i), dbl_mb(k_sc), 
     &        vc_psi(i))
      enddo

      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)

      end
c===============================================================
      

c===============================================================
c===============================================================
      subroutine wfat_2a1e_collect_vcpsi_ord1(geom, sphi, m, nr, nt,
     &     nr_m, nt_m, nr_p, nt_p, nptp, ke, pottype, nel, potparam, zi,
     &     vcf, cfx, cfy, cfz, chpos, vpt, vc_psi, vc_cpsi)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                     ! [input],
      integer sphi                     ! [input],
      integer m                        ! [input],
      integer nr                       ! [input],
      integer nt                       ! [input],
      integer nr_m                     ! [input],
      integer nt_m                     ! [input],
      integer nr_p                     ! [input],
      integer nt_p                     ! [input],
      integer nptp                     ! [input],
      double precision ke              ! [input],
      character*3 pottype              ! [input],
      integer nel                      ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)     ! [input], 1)
      double precision zi              ! [input],
      double precision vcf(nr*nt)      ! [input],
      double precision cfx(1:*)        ! [input],
      double precision cfy(1:*)        ! [input],
      double precision cfz(nr*nt)      ! [input],
      double precision chpos(3)        ! [input],
      double precision vpt(3,nptp)     ! [input],
      double precision vc_psi(nptp)    ! [output],
      double precision vc_cpsi(nptp,3) ! [output],

      integer nscr, l_sc, k_sc, i

      nscr = 
     &     (nr + nt + nr*nt) +
     &     (nr_p + nt_p + nr_p*nt_p)
      if (m >= 1) nscr = nscr + (nr_m + nt_m + nr_m*nt_m)
      if (.not. MA_alloc_get(MT_DBL, nscr, 'vscr', l_sc, k_sc))
     &     call errquit('WFAT: vscr allocation failed', 59, MA_ERR)
      
      do i = 1, nptp
         call wfat_2a1e_vcpsi_ord1(geom, sphi, m, nr, nt, nr_m, nt_m, 
     &        nr_p, nt_p, ke, pottype, nel, potparam, zi, vcf, cfx, cfy,
     &        cfz, chpos, vpt(:,i), dbl_mb(k_sc), vc_psi(i),
     &        vc_cpsi(i,:))
      enddo

      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)
      
      end
c===============================================================
      
