c=======================================================
c=======================================================
      subroutine wfat_calc_tmd(kx, ky, nch, erg, nqn, nid, mid,
     &                         ef, srfc, tmd)

      implicit none
#include "geom.fh"
c     #include "cdft.fh"
#include "wfat_qm_param.fh"      
#include "errquit.fh"      
      double precision kx
      double precision ky
      integer nch
      double precision erg(nch)
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      double precision ef
      double complex srfc(1:*)    ! [input] srfc = strucutre factor, dimension [natom x nqn x nch]
      double complex tmd(1:*)     ! [output] tmd = transverse momentum distribution, dimension [natom x nch]

      integer natom, ntmd
      double precision kk, phk, pi_const
      double complex imnum
      parameter ( pi_const = 3.14159265359d0 )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )


c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom_hd, natom))
     &     call errquit('WFAT: Getting natom failed', 1, GEOM_ERR)      

      kk = kx*kx + ky*ky
      phk = atan2(ky, kx)

      call wfat_asymp_xi(kk/ef, phk, nch, erg, nqn, nid, mid,
     &                   natom, zion, ef, srfc(1), tmd(1))

      ntmd = natom * nch
      tmd(1:ntmd) = (sqrt(2.0d0))**3 * pi_const * imnum / sqrt(ef) *
     &              tmd(1:ntmd)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_asymp_xi(xi, phi, nch, erg, nqn, nid, mid,
     &                         natom, zion, ef, srfc, asx)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"      
      double precision xi
      double precision phi
      integer nch
      double precision erg(nch)
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      integer natom
      double precision zion
      double precision ef
      double complex srfc(natom, nqn, nch)
      double complex asx(natom, nch)

      integer nadb, ich, iqn, n, m, l_ad, k_ad
      double precision ffc
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


      nadb = nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nadb, 'adb', l_ad, k_ad))
     &     call errquit('WFAT: adb allocation failed', 1, MA_ERR)
      call wfat_adbas(xi, phi, nch, erg, nqn, nid, mid, dcpl_mb(k_ad))

      do ich = 1, nch
         asx(:,ich) = czero
         do iqn = 1, nqn
            n = nid(iqn)
            m = mid(iqn)
            call wfat_field_factor(ef, erg(ich), zion, n, m, ffc)
    
            !!== Atom-resolved asymptotic orbital ==!!
            asx(:,ich) = asx(:,ich) + sqrt(ffc) * srfc(:,iqn,ich) *
     &                   dcpl_mb(k_ad + (ich-1)*nqn + iqn - 1)
         end do
      end do

      if (.not. MA_free_heap(l_ad))
     &     call errquit('WFAT: adb deallocation failed', 2, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_adbas(xi, phi, nch, erg, nqn, nid, mid, adb)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"      
      integer nch, nqn, nid(nqn), mid(nqn)
      double precision xi, phi, erg(nch)
      double complex adb(nqn, nch)
    
      integer ich, iqn, l_la, k_la
      integer n, m, nlag
      double precision kp, s(1), pi_const
      double complex imnum

      parameter ( pi_const = 3.14159265359d0 )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )
      


      nlag = maxval(nid) + 1
      if (.not. MA_alloc_get(MT_DBL, nlag, 'lag', l_la, k_la))
     &     call errquit('WFAT: lag allocation failed', 1, MA_ERR)
      
      do ich = 1, nch
         kp = sqrt(2*abs(erg(ich)))
         do iqn = 1, nqn

            s = kp * xi
            n = nid(iqn)
            m = mid(iqn)
            call wfat_lagfun_norm(1, n, dble(abs(m)), s, dbl_mb(k_la))

            adb(iqn, ich) = sqrt(kp) * s(1)**(dble(abs(m))/2) *
     &                      exp(-s(1)/2) * dbl_mb(k_la + n) * 
     &                      exp(imnum*m*phi) / sqrt(2*pi_const)

         enddo
      enddo

      if (.not. MA_free_heap(l_la))
     &     call errquit('WFAT: lag deallocation failed',
     &     2, MA_ERR)
      

      end
c=======================================================


c================================================
c================================================
      subroutine wfat_get_kmax(ef, nch, erg, nqn, nid, mid, kmax)
      
      implicit none
#include "wfat.fh"      
      integer nch, nqn, nid(nqn), mid(nqn)
      double precision ef, erg(nch)
      double precision kmax
      
      integer ich
      double precision kmax_, kp
      
      
      kmax = 0.0d0
      do ich = 1, nch
         kp = sqrt(2*abs(erg(ich)))
      
         if (maxval(nid) == 0 .and. maxval(abs(mid)) == 0) then
            kmax_ = sqrt(2 * ef / kp * -log(KMAX00))
         else
            kmax_ = sqrt(ef * maxval(abs(mid) + 2*nid) / kp)
         end if
      
         if (kmax_ > kmax) kmax = kmax_
      end do
      kmax = KMAX_FACTOR * kmax
        
      end
c================================================
      
