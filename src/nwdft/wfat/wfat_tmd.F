c=======================================================
c=======================================================
      subroutine wfat_calc_tmd0(rtdb, beta, gamma, nat, nch, nqn, nid,
     &     mid, erg, zi, ef, gn)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer rtdb
      double precision beta
      double precision gamma
      integer nat
      integer nch
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      double precision erg(nch)
      double precision zi
      double precision ef
      double complex gn(nat, nqn, nch)

      integer l_fn, k_fn, i_fn, n, m, i, iqn, nkx, nky, ix, iy
      double precision kx0, kxn, kx, dkx, ky0, kyn, ky, dky, pk, ph
      double precision pi_const, phi0, pk_part(nch)
      double complex ak, czero, imnum

      parameter ( pi_const = 3.14159265359d0 )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)),
     &            imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)))      
      

      call wfat_get_momenta(rtdb, nkx, kx0, kxn, nky, ky0, kyn)
      

      if (.not. MA_alloc_get(MT_DCPL, nqn*nch, 'fn', l_fn, k_fn))
     &     call errquit('WFAT: fn allocation failed', 93, MA_ERR)

      do i = 1, nch
         do iqn = 1, nqn
            n = nid(iqn)
            m = mid(iqn)
            i_fn = k_fn + (i-1)*nqn + (iqn-1)
            call wfat_get_ascf0(beta, gamma, n, m, zi, erg(i),
     &           sum(gn(:,iqn,i)), ef, dcpl_mb(i_fn))
         enddo
      enddo

      
      dkx = (kxn-kx0) / (nkx-1)
      kx = kx0 - dkx
      dky = (kyn-ky0) / (nky-1)
      ky = ky0 - dky
      do ix = 1, nkx
         kx = kx + dkx
         do iy = 1, nky
            ky = ky + dky

            pk = 0.0d0
            do i = 1, nch
               ak = czero
               do iqn = 1, nqn
                  n = nid(iqn)
                  m = mid(iqn)

                  ph = atan2(ky,kx)
                  call wfat_calc_phi0xi((kx**2+ky**2)/ef, n, m, erg(i),
     &                 phi0)
                  i_fn = k_fn + (i-1)*nqn + (iqn-1)
                  ak = ak + dcpl_mb(i_fn) * phi0 * exp(imnum*m*ph)
     &                 / sqrt(2*pi_const)
               enddo
               ak = sqrt(8/ef) * pi_const * imnum * ak
               pk_part(i) = (abs(ak))**2
               pk = pk + pk_part(i)
            enddo

c            call wfat_print_tmd(nch, pk_part, pk)
         enddo
      enddo


      if (.not. MA_free_heap(l_fn))
     &     call errquit('WFAT: fn deallocation failed', 89, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_tmd1(rtdb, beta, gamma, nat, nch, nqn, nid,
     &     mid, erg, chmmt, pb, zi, ef, wfat0, on, gn, jn, kn)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer rtdb
      double precision beta
      double precision gamma
      integer nat
      integer nch
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      double precision erg(nch)
      double precision chmmt(3, nch)
      double precision pb(3, 3, nch)
      double precision zi
      double precision ef
      logical wfat0(nqn, nch)
      double precision on(2, nqn, nch)
      double complex gn(nat, nqn, nch)
      double complex jn(nat, nqn, nch)
      double complex kn(nat, nqn, nch)

      integer l_fn0, k_fn0, i_fn0, l_fn1, k_fn1, i_fn1, n, m, i, iqn
      integer nkx, nky, ix, iy
      double precision kx0, kxn, kx, dkx, ky0, kyn, ky, dky, pk, ph
      double precision pi_const, phi0, phi1, pk_part(nch)
      double complex ak, czero, imnum, fn

      parameter ( pi_const = 3.14159265359d0 )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)),
     &            imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)))      
      

      call wfat_get_momenta(rtdb, nkx, kx0, kxn, nky, ky0, kyn)
      

      if (.not. MA_alloc_get(MT_DCPL, nqn*nch, 'fn0', l_fn0, k_fn0))
     &     call errquit('WFAT: fn0 allocation failed', 93, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nqn*nch, 'fn1', l_fn1, k_fn1))
     &     call errquit('WFAT: fn1 allocation failed', 93, MA_ERR)

      do i = 1, nch
         do iqn = 1, nqn
            n = nid(iqn)
            m = mid(iqn)

            i_fn1 = k_fn1 + (i-1)*nqn + (iqn-1)
            call wfat_get_ascf1(beta, gamma, n, m, zi, erg(i),
     &           chmmt(:,i), pb(:,:,i), on(:,iqn,i), sum(gn(:,iqn,i)), 
     &           sum(jn(:,iqn,i)), sum(kn(:,iqn,i)), ef, dcpl_mb(i_fn1))
            i_fn0 = k_fn0 + (i-1)*nqn + (iqn-1)
            call wfat_get_ascf0(beta, gamma, n, m, zi, erg(i),
     &           sum(gn(:,iqn,i)), ef, dcpl_mb(i_fn0))
         enddo
      enddo

      
      dkx = (kxn-kx0) / (nkx-1)
      kx = kx0 - dkx
      dky = (kyn-ky0) / (nky-1)
      ky = ky0 - dky
      do ix = 1, nkx
         kx = kx + dkx
         do iy = 1, nky
            ky = ky + dky

            pk = 0.0d0
            do i = 1, nch
               ak = czero
               do iqn = 1, nqn
                  n = nid(iqn)
                  m = mid(iqn)

                  ph = atan2(ky,kx)
                  call wfat_calc_phi0xi((kx**2+ky**2)/ef, n, m, erg(i),
     &                 phi0)                  
                  
                  if (wfat0(iqn,i)) then
                     i_fn0 = k_fn0 + (i-1)*nqn + (iqn-1)
                     fn = dcpl_mb(i_fn0)
                     phi1 = 0.0d0
                  else
                     i_fn1 = k_fn1 + (i-1)*nqn + (iqn-1)
                     fn = dcpl_mb(i_fn1)
                     call wfat_calc_phi1xi((kx**2+ky**2)/ef, n, m, 
     &                    erg(i), beta, gamma, chmmt(:,i), phi1)
                  endif
                  
                  ak = ak + fn * (phi0+ef*phi1) * exp(imnum*m*ph)
     &                 / sqrt(2*pi_const)
               enddo
               ak = sqrt(8/ef) * pi_const * imnum * ak
               pk_part(i) = (abs(ak))**2
               pk = pk + pk_part(i)
            enddo

c            call wfat_print_tmd(nch, pk_part, pk)
         enddo
      enddo


      if (.not. MA_free_heap(l_fn0))
     &     call errquit('WFAT: fn0 deallocation failed', 89, MA_ERR)
      if (.not. MA_free_heap(l_fn1))
     &     call errquit('WFAT: fn1 deallocation failed', 89, MA_ERR)
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_get_momenta(rtdb, nkx, kx0, kxn, nky, ky0, kyn)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb           ! [input]
      integer nkx            ! [input]
      double precision kx0   ! [output]
      double precision kxn   ! [output]
      integer nky            ! [output]
      double precision ky0   ! [output]
      double precision kyn   ! [output]


      if (.not. rtdb_get(rtdb, 'wfat:kx_start', mt_dbl, 1, kx0))
     &     call errquit('WFAT: Failed reading kx_start from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:kx_last', mt_dbl, 1, kxn))
     &     call errquit('WFAT: Failed reading kx_last from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:nkx', mt_int, 1, nkx))
     &     call errquit('WFAT: Failed reading nkx from rtdb',
     &     1, INPUT_ERR)

      if (.not. rtdb_get(rtdb, 'wfat:ky_start', mt_dbl, 1, ky0))
     &     call errquit('WFAT: Failed reading ky_start from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:ky_last', mt_dbl, 1, kyn))
     &     call errquit('WFAT: Failed reading ky_last from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:nky', mt_int, 1, nky))
     &     call errquit('WFAT: Failed reading nky from rtdb',
     &     1, INPUT_ERR)
      
      end
c=======================================================
      
