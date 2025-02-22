c============================================================
c============================================================
      subroutine wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid, mid,
     &     ef, erg, zi, gn, jn, kn, on, chmmt, pb, ffc, wfat0, interm, 
     &     prate, rate)

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
      logical interm                         ! [input],
      double precision prate(nqn,nch)        ! [output],
      double precision rate                  ! [output],

      integer i, iqn, n, m
      double precision dummy


      if (ga_nodeid() == 0 .and. interm) then
         write(luout, '("")')
         write(luout, '("ATTENTION: Printing intermediate '//
     &        'quantities. A_nu, B_nu, and some other intermediate '//
     &        'quantities will be printed below.")')
      endif
       
      
      do i = 1, nch
         do iqn = 1, nqn
            n = nid(iqn)
            m = mid(iqn)

            if (wfat0(iqn,i)) then
               call wfat_calc_rate(nat, 1, 1, gn(:,iqn,i), 
     &              ffc(iqn,i), prate(iqn,i), dummy)
            else
               if (ga_nodeid() == 0 .and. interm) then
                  write(luout, '("Channel = ", i3)') i
                  write(luout, '("para. quantum no. (n,m) = ", 2i3)')
     &                 n, m
               endif
               call wfat_wfat1_rate(bt, gm, nat, 1, 1, n, m, ef, erg(i),
     &              zi, gn(:,iqn,i), jn(:,iqn,i), kn(:,iqn,i), 
     &              on(:,iqn,i), chmmt(:,i), pb(:,:,i), ffc(iqn,i), 
     &              interm, prate(iqn,i), dummy)
            endif

         enddo
      enddo

      rate = sum(prate)

      if (ga_nodeid() == 0 .and. interm) write(luout, '("")')
       
      end
c============================================================

      
c============================================================
c============================================================
      subroutine wfat_wfat1_rate(bt, gm, nat, nqn, nch, nid, mid,
     &     ef, erg, zi, gn, jn, kn, on, chmmt, pb, ffc, interm, p_rt,
     &     rt)

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
      logical interm                         ! [input],
      double precision p_rt(nqn,nch)         ! [output],
      double precision rt                    ! [output],
      
      integer i, iqn, n, m
      double precision on_, protmat(3,3), mu(3), lpbz, kp, gmm
      double precision bt0, bt1_1, bt1_2, bt1, aa, bb_, bb2
      double complex gn_, an_, an2

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      

      do i = 1, nch
         mu = matmul(transpose(protmat), -chmmt(:,i))
         call wfat_calc_lfpolb(protmat, pb(:,:,i), lpbz)
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
            bb_ = -kp*lpbz - mu(3)**2/kp + mu(3)/kp**2 +
     &            4*mu(3)*bt0/kp**3 - (10+18*gmm+3*gmm**2)/24/kp**3 -
     &            (9-6*gmm)*bt0/4/kp**4 - (49+2*gmm)*bt0**2/8/kp**5 +
     &            3*bt0**3/2/kp**6 - bt0**4/8/kp**7

cfordebug            write(luout, '("B_nu - 2 Re(a_nu) = ", f15.7)') bb_
            on_ = on(1,iqn,i) + mu(3)*on(2,iqn,i)
            gn_ = sum(gn(:,iqn,i))
            an_ = sum(jn(:,iqn,i)) + sum(kn(:,iqn,i))
            p_rt(iqn,i) = ffc(iqn,i) *
     &           (abs(gn_)**2 * 
     &           (1+aa*ef*log(ef/4/kp**2)+bb_*ef+2*ef*on_) +
     &           2*ef*real(conjg(gn_)*an_))


c           ==== These lines below are only for debugging ====
            an2 = on_ + an_/gn_
            bb2 = bb_ + 2 * real(an2)
            if (ga_nodeid() == 0 .and. interm) then
               write(luout, '("(beta, gamma) = ", 2f15.6)') bt, gm
               write(luout, '("   a_nu = ", 2es15.7)') an2
               write(luout, '("   A_nu = ", es15.7)') aa
               write(luout, '("   B_nu = ", es15.7)') bb2
               write(luout, '("   alpha_zz = ", es15.7)') lpbz
            endif
            
         enddo
      enddo
      
      rt = sum(p_rt)
      
      end
c============================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_wfat1_yield(bt, gm, zi, nlas, ftlas, alas, 
     &                 wlas, plas, tlas, dlas, nt, tmax, extf, nch, erg,
     &                 ergf, nqn, nid, mid, nat, gn, jn, kn, gnn, jnn, 
     &                 knn, on, chmmt, pb, wfat0, yield)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "wfat_field_param.fh"
      double precision bt           ! [input],
      double precision gm           ! [input],
      double precision zi
      integer nlas
      character(len=FLD_TP_LEN) ftlas(nlas)       ! Function type of envelope, 'gaussian' or 'hann'         
      double precision alas(nlas)   ! Peak electric field, in a.u. of e. field.                              
      double precision wlas(nlas)   ! angular frequency, in a.u. (Hartree).                                      
      double precision plas(nlas)   ! Carrier-envelope phase, rad                          
      double precision tlas(nlas)   ! Delay (temporal position of the envelope peak), a.u. of time
      double precision dlas(nlas)   ! Field FWHM, a.u. of time                                       
      integer nt                    
      double precision tmax         ! Time window width, a.u. of time
      logical extf
      integer nch
      double precision erg(nch)
      double precision ergf(nch)
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
      double precision t0, dt, t, las, rate
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
         
c        ==== Evaluate laser field ====
         if (extf) then
            call wfat_eval_extlaser(i, las)
         else
            call wfat_eval_laser(nlas, ftlas, alas, wlas, plas, tlas, 
     &           dlas, t, las)
         endif
         
c        ==== Calculate field factors ====
         do ich = 1, nch
            do iqn = 1, nqn
               i_fc = k_fc + (ich-1)*nqn + (iqn-1)
               if (las >= 1.0d-13) then
                  call wfat_field_factor(las, ergf(ich), zi, nid(iqn),
     &                 mid(iqn), dbl_mb(i_fc))
               else if (las <= -1.0d-13) then
                  call wfat_field_factor(abs(las), ergf(ich), zi, 
     &                 nid(iqn), mid(iqn), dbl_mb(i_fc))
               else
                  dbl_mb(i_fc) = 0.0d0
               endif
            enddo
         enddo
         
c        ==== Calculate rates ====
         dbl_mb(k_pi:k_pi+npi-1) = 0.0d0
         if (las >= 1.0d-13) then
            call wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid, mid,
     &           las, erg, zi, gn, jn, kn, on, chmmt, pb, dbl_mb(k_fc),
     &           wfat0, .false., dbl_mb(k_pi), rate)
         else if (las <= -1.0d-13) then
            call wfat_calc_wfat1_rate(180.0d0-bt, 180.0d0+gm, nat, nqn,
     &           nch, nid, mid, abs(las), erg, zi, gnn, jnn, knn, on, 
     &           chmmt, pb, dbl_mb(k_fc), wfat0, .false., dbl_mb(k_pi),
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

      
      
      
