c============================================================
c============================================================
      subroutine wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid, mid,
     &     ef, erg, zi, gn, jn, kn, on, chmmt, pb, ffc, wfat0, prate, 
     &     rate)

      implicit none
#include "stdio.fh"
#include "global.fh"
      double precision bt                    ! [input],
      double precision gm                    ! [input],
      integer nat                            ! [input],
      integer nqn                            ! [input],
      integer nch                            ! [input],
      integer nid(nqn)                       ! [input],
      integer mid(nqn)                       ! [input],
      double precision ef                    ! [input],
      double precision erg(nch)              ! [input],
      double precision zi                    ! [input],
      double complex gn(nat, nqn, nch)       ! [input],
      double complex jn(nat, nqn, nch)       ! [input],
      double complex kn(nat, nqn, nch)       ! [input],
      double precision on(2, nqn, nch)       ! [input],
      double precision chmmt(3,nch)          ! [input],
      double precision pb(3,3,nch)           ! [input],
      double precision ffc(nqn, nch)         ! [input],
      logical wfat0(nqn,nch)                 ! [input],
      double precision prate(nqn,nch)        ! [output],
      double precision rate                  ! [output],

      integer i, iqn, n, m
      double precision dummy


      do i = 1, nch
         do iqn = 1, nqn
            n = nid(iqn)
            m = mid(iqn)

            if (wfat0(iqn,i)) then
               call wfat_calc_rate(nat, 1, 1, gn(:,iqn,i), ffc(iqn,i),
     &              prate(iqn,i), dummy)
            else
               call wfat_wfat1_rate(bt, gm, nat, 1, 1, n, m, ef, erg(i),
     &              zi, gn(:,iqn,i), jn(:,iqn,i), kn(:,iqn,i), 
     &              on(:,iqn,i), chmmt(:,i), pb(:,:,i), ffc(iqn,i), 
     &              prate(iqn,i), dummy)
            endif

         enddo
      enddo

      rate = sum(prate)
       
      end
c============================================================

      
c============================================================
c============================================================
      subroutine wfat_wfat1_rate(bt, gm, nat, nqn, nch, nid, mid,
     &     ef, erg, zi, gn, jn, kn, on, chmmt, pb, ffc, p_rt, rt)
c     DESCRIPTION:
c     gn :
c       $\exp(-\kappa\mu_z) g_\nu$ where $g_\nu$ are structure factors for all
c       channels and quantum numbers resolved into atomic contributions based
c       on Becke partition. The exponential prefactor with the z component of
c       dipole vector in the argument is
c     jn and kn :
c       First-order WFAT equivalents of gn.
      implicit none
#include "stdio.fh"
#include "global.fh"
      double precision bt                    ! [input],
      double precision gm                    ! [input],
      integer nat                            ! [input],
      integer nqn                            ! [input],
      integer nch                            ! [input],
      integer nid(nqn)                       ! [input],
      integer mid(nqn)                       ! [input],
      double precision ef                    ! [input],
      double precision erg(nch)              ! [input],
      double precision zi                    ! [input],
      double complex gn(nat, nqn, nch)       ! [input],
      double complex jn(nat, nqn, nch)       ! [input],
      double complex kn(nat, nqn, nch)       ! [input],
      double precision on(2, nqn, nch)       ! [input],
      double precision chmmt(3,nch)          ! [input],
      double precision pb(3,3,nch)           ! [input],
      double precision ffc(nqn, nch)         ! [input],
      double precision p_rt(nqn,nch)         ! [output],
      double precision rt                    ! [output],
      
      integer i, iqn, n, m
      double precision on_, protmat(3,3), mu(3), pbl(3,3), kp, gmm
      double precision bt0, bt1_1, bt1_2, bt1, aa, bb_
      double complex gn_, an_

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      

      do i = 1, nch
         mu = matmul(transpose(protmat), -chmmt(:,i))
c         call wfat_represent_sym(3, 3, protmat, pb, pbl) ! pbl = polarizability tensor in lab frame.
         call wfat_represent_sym(3, 3, protmat, pb(:,:,i), pbl)   ! pbl = polarizability tensor in lab frame.
cOLD         call wfat_calc_lfpolb(protmat, pb(:,:,i), lpbz)
         kp = sqrt( 2 * abs(erg(i)) )
         
         do iqn = 1, nqn
            n = nid(iqn)
            m = mid(iqn)

            gmm = (1.0d0 - m**2) / 4.0d0
            bt0 = zi - kp*(n + (abs(m)+1)/2.0d0)
            bt1_1 = -1.0d0/4/kp**2
     &              * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
            bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))
            bt1 = bt1_1 + mu(3)*bt1_2
            aa = -2*bt1/kp - (gmm-2*mu(3)*bt0)/kp**3 - 3*bt0**2/kp**5
            bb_ = -kp*pbl(3,3) - mu(3)**2/kp + mu(3)/kp**2 +
     &            4*mu(3)*bt0/kp**3 - (10+18*gmm+3*gmm**2)/24/kp**3 -
     &            (9-6*gmm)*bt0/4/kp**4 - (49+2*gmm)*bt0**2/8/kp**5 +
     &            3*bt0**3/2/kp**6 - bt0**4/8/kp**7
c           bb_ is equal to \tilde B_\nu defined in eq. (18b) of this paper: http://dx.doi.org/10.1103/PhysRevA.93.033412.
            
cfordebug            write(luout, '("B_nu - 2 Re(a_nu) = ", f15.7)') bb_
            on_ = on(1,iqn,i) + mu(3)*on(2,iqn,i)
            gn_ = sum(gn(:,iqn,i))    ! sum over atomic contributions
            an_ = sum(jn(:,iqn,i)) + sum(kn(:,iqn,i))    ! sum over atomic contributions
c           an_ is equal to h_\nu defined in eq. (12) of this paper: https://doi.org/10.1103/4xqd-yhdt.
            p_rt(iqn,i) = ffc(iqn,i) *
     &           (abs(gn_)**2 * 
     &           (1+aa*ef*log(ef/4/kp**2)+bb_*ef+2*ef*on_) +
     &           2*ef*real(conjg(gn_)*an_))
         enddo
      enddo
      
      rt = sum(p_rt)
      
      end
c============================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_wfat1_yield(bt, gm, zi, tmax, nt, las, nch, 
     &                 erg, nqn, nid, mid, nat, gn, jn, kn, gnn, 
     &                 jnn, knn, on, chmmt, pb, wfat0, yield)
      
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "wfat_field_param.fh"
      double precision bt           ! [input],
      double precision gm           ! [input],
      double precision zi
      double precision tmax
      integer nt
      double precision las(nt)
      integer nch
      double precision erg(nch)
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      integer nat
      double complex gn(nat, nqn, nch)       ! [input],
      double complex jn(nat, nqn, nch)       ! [input],
      double complex kn(nat, nqn, nch)       ! [input],
      double complex gnn(nat, nqn, nch)      ! [input],
      double complex jnn(nat, nqn, nch)      ! [input],
      double complex knn(nat, nqn, nch)      ! [input],
      double precision on(2, nqn, nch)       ! [input],
      double precision chmmt(3,nch)          ! [input],
      double precision pb(3,3,nch)           ! [input],
      logical wfat0(nqn,nch)                 ! [input],
      double precision yield                 ! [output],
      
      integer i, j, k, npi, l_pi, k_pi, nfc, l_fc, k_fc, i_fc, ich, iqn
      double precision t0, dt, t, rate
      parameter (t0=0.0d0)

      
      npi = nqn * nch
      if (.not. ma_alloc_get(mt_dbl, npi, 'partial_ion', l_pi, k_pi))
     &     call errquit('WFAT: Cannot allocate partial_ion', 1, MA_ERR)
      nfc = nqn * nch
      if (.not. MA_alloc_get(MT_DBL, nfc, 'ffc', l_fc, k_fc))
     &     call errquit('WFAT: ffc allocation failed', 2, MA_ERR)
      
      
      dt = (tmax-t0)/(nt-1)

      yield = 0.0d0
      do i = 1, nt
         t = t0 + (i-1)*dt
         
c        ==== Calculate field factors ====
         do ich = 1, nch
            do iqn = 1, nqn
               i_fc = k_fc + (ich-1)*nqn + (iqn-1)
               if (las(i) >= 1.0d-13) then
                  call wfat_field_factor(las(i), erg(ich), zi, 
     &                 nid(iqn), mid(iqn), dbl_mb(i_fc))
               else if (las(i) <= -1.0d-13) then
                  call wfat_field_factor(abs(las(i)), erg(ich), zi, 
     &                 nid(iqn), mid(iqn), dbl_mb(i_fc))
               else
                  dbl_mb(i_fc) = 0.0d0
               endif
            enddo
         enddo
         
c        ==== Calculate rates ====
         dbl_mb(k_pi:k_pi+npi-1) = 0.0d0
         if (las(i) >= 1.0d-13) then
            call wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid, mid,
     &           las(i), erg, zi, gn, jn, kn, on, chmmt, pb, 
     &           dbl_mb(k_fc), wfat0, dbl_mb(k_pi), rate)
         else if (las(i) <= -1.0d-13) then
            call wfat_calc_wfat1_rate(180.0d0-bt, 180.0d0+gm, nat, nqn,
     &           nch, nid, mid, abs(las(i)), erg, zi, gnn, jnn, knn, on,
     &           chmmt, pb, dbl_mb(k_fc), wfat0, dbl_mb(k_pi),
     &           rate)
         else
            dbl_mb(k_pi : k_pi+npi-1) = 0.0d0
            rate = 0.0d0
         endif
         
c        ==== Calculate yield ====
         yield = yield + rate * dt
      enddo

      yield = 1.0d0 - exp(-yield)

      if (.not. ma_free_heap(l_pi))
     &     call errquit('WFAT: Cannot deallocate partial_ion',
     &     1, MA_ERR)
      if (.not. MA_free_heap(l_fc))
     &     call errquit('WFAT: ffc deallocation failed',
     &     3, MA_ERR)
      
      end
c=======================================================

      
      
      
