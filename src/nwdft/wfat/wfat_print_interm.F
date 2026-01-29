c============================================================
c============================================================
      subroutine wfat_initprint_interm(nch, orbid0, spinid,
     &     nqn, nid, mid, wfat0)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer nch, orbid0(nch), spinid(nch), nqn, nid(nqn), mid(nqn)
      logical wfat0(nqn,nch)

      character*6 itm_suff

      if (WF_ITM_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_interm cannot be '//
     &     'called when the previous call to this routine has not '//
     &     'been terminated by calling wfat_endprint_interm.',
     &     1, UERR)


      itm_suff = 'interm'
      call util_file_name(itm_suff, .false., .false., WF_ITM_FNAME)
      open(WF_ITM_UNIT, file=WF_ITM_FNAME, action='write')
      WF_ITM_PRINT_READY = .true.
      
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout,
     &     '(5x, "Printing WFAT(1) intermediate quantities to ", a, a)')
     &     ''''//trim(adjustl(WF_ITM_FNAME))//'''',
     &     ' Please wait ...'

      call wfat_print_itmline(nch, nqn, wfat0)
      call wfat_print_itmcoln(nch, nqn, wfat0)
      call wfat_print_itmorbtitle(nch, orbid0, spinid, nqn, wfat0)
      call wfat_print_itmqnttitle(nch, nqn, nid, mid, wfat0)
      call wfat_print_itmitmtitle(nch, nqn, wfat0)
      call wfat_print_itmline(nch, nqn, wfat0)
      
      end
c============================================================

      
c============================================================
c============================================================
      subroutine wfat_print_interm(bt, gm, nat, nqn, nch, nid, mid,
     &     erg, zi, gn, jn, kn, on, chmmt, pb, wfat0, gap)
c     DESCRIPTION:
c     gn :
c       $\exp(-\kappa\mu_z) g_\nu$ where $g_\nu$ are structure factors for all
c       channels and quantum numbers resolved into atomic contributions based
c       on Becke partition. The exponential prefactor with the z component of
c       dipole vector in the argument is
c     jn and kn :
c       First-order WFAT equivalents of gn.
      implicit none
#include "wfat.fh"
#include "stdio.fh"
#include "global.fh"
      double precision bt                    ! [input],
      double precision gm                    ! [input],
      integer nat                            ! [input],
      integer nqn                            ! [input],
      integer nch                            ! [input],
      integer nid(nqn)                       ! [input],
      integer mid(nqn)                       ! [input],
      double precision erg(nch)              ! [input],
      double precision zi                    ! [input],
      double complex gn(nat, nqn, nch)       ! [input],
      double complex jn(nat, nqn, nch)       ! [input],
      double complex kn(nat, nqn, nch)       ! [input],
      double precision on(2, nqn, nch)       ! [input],
      double precision chmmt(3,nch)          ! [input],
      double precision pb(3,3,nch)           ! [input],
      logical wfat0(nqn,nch)                 ! [input],
      logical gap                            ! [input],
      
      integer i, iqn, n, m, j
      double precision on_, protmat(3,3), mu(3), pbl(3,3), kp, gmm
      double precision bt0, bt1_1, bt1_2, bt1, aa, bb_, bb2
      double precision itm(WF_NITM)
      double complex gn_, an_, an2

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      

c     ==== beta and gamma ====
      write(WF_ITM_UNIT, '(1x, 2f15.6)', advance='no') bt, gm
      do i = 1, nch
         mu = matmul(transpose(protmat), -chmmt(:,i))
c         call wfat_represent_sym(3, 3, protmat, pb, pbl) ! pbl = polarizability tensor in lab frame.
         call wfat_represent_sym(3, 3, protmat, pb(:,:,i), pbl)   ! pbl = polarizability tensor in lab frame.
cOLD         call wfat_calc_lfpolb(protmat, pb(:,:,i), lpbz)
         kp = sqrt( 2 * abs(erg(i)) )
         
         do iqn = 1, nqn
            if (wfat0(iqn,i)) cycle
            write(WF_ITM_UNIT, WF_QNTGAP, advance='no')    ! horizontal gap between channel and paraqn combination.
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
            on_ = on(1,iqn,i) + mu(3)*on(2,iqn,i)
            gn_ = sum(gn(:,iqn,i))    ! sum over atomic contributions
            an_ = sum(jn(:,iqn,i)) + sum(kn(:,iqn,i))    ! sum over atomic contributions
c           an_ is equal to h_\nu defined in eq. (12) of this paper: https://doi.org/10.1103/4xqd-yhdt.
            an2 = on_ + an_/gn_
c           an2 is equal to a_\nu defined in eq. (60) of this paper: http://dx.doi.org/10.1103/PhysRevA.93.033412.
            bb2 = bb_ + 2 * real(an2)
c           bb2 is equal to B_\nu defined immediately below eq. (24) of this paper: http://dx.doi.org/10.1103/PhysRevA.93.033412.

            itm = (/ real(an2), imag(an2), aa, bb_, bb2, pbl(3,3), on_,
     &           real(conjg(gn_)*sum(jn(:,iqn,i))),
     &           real(conjg(gn_)*sum(kn(:,iqn,i))),
     &           real(gn_), imag(gn_) /)
            do j = 1, WF_NITM
               write(WF_ITM_UNIT, '(1x, es15.7)', advance='no') itm(j)
            enddo
            
         enddo
      enddo

      write(WF_ITM_UNIT, '(" ")') 
      if (gap) write(WF_ITM_UNIT, '(" ")')
      call flush(WF_ITM_UNIT)
      
      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_print_itmcoln(nch, nqn, wfat0)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer nch, nqn
      logical wfat0(nqn,nch)

      character*10 cols
      integer col, i, j, k

      if (.not. WF_ITM_PRINT_READY) 
     &     call errquit('WFAT: wfat_initprint_interm must first be '//
     &     'called before calling wfat_print_itmcoln.', 1, UERR)

      write(WF_ITM_UNIT, '("#", 2a15)', advance='no')
     &     'col. #1', 'col. #2'
      col = 3
      do i = 1, nch
         do j = 1, nqn
            if (wfat0(j,i)) cycle
            write(WF_ITM_UNIT, WF_QNTGAP, advance='no')    ! horizontal gap between channel and paraqn combination.
            do k = 1, WF_NITM               
               write(cols, '(i10)') col
               write(WF_ITM_UNIT, '(1x, a15)', advance='no')
     &              'col. #'//trim(adjustl(cols))
               col = col + 1
            enddo
         enddo
      enddo
      write(WF_ITM_UNIT, '("")')

      end
c============================================================
      

c============================================================
c============================================================
      subroutine wfat_print_itmorbtitle(nch, oid0, sid, nqn, wfat0)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer nch, oid0(nch), sid(nch), nqn
      logical wfat0(nqn,nch)

      integer io, is, i, j, k

      if (.not. WF_ITM_PRINT_READY) 
     &     call errquit('WFAT: wfat_initprint_interm must first be '//
     &     'called before calling wfat_print_itmorbtitle.', 1, UERR)

      
      write(WF_ITM_UNIT, '("#", 2a15)', advance='no')
     &     'beta (deg)', 'gamma (deg)'
      do i = 1, nch
         io = oid0(i)
         is = sid(i)
         do j = 1, nqn
            if (wfat0(j,i)) cycle
            write(WF_ITM_UNIT, WF_QNTGAP, advance='no')    ! horizontal gap between channel and paraqn combination.
            do k = 1, WF_NITM
               write(WF_ITM_UNIT, '(1x)', advance='no')
               call wfat_print_orbtitle_fx(WF_ITM_UNIT, is, io, '(a15)')
            enddo
         enddo
      enddo
      write(WF_ITM_UNIT, '("")')


      end
c============================================================
      

c============================================================
c============================================================
      subroutine wfat_print_itmqnttitle(nch, nqn, nid, mid, wfat0)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer nch, nqn, nid(nqn), mid(nqn)
      logical wfat0(nqn,nch)

      integer n, m, i, j, k

      if (.not. WF_ITM_PRINT_READY) 
     &     call errquit('WFAT: wfat_initprint_interm must first be '//
     &     'called before calling wfat_print_itmqnttitle.', 1, UERR)
      
      write(WF_ITM_UNIT, '("#", 2a15)', advance='no') ' ', ' '
      do i = 1, nch
         do j = 1, nqn
            if (wfat0(j,i)) cycle
            n = nid(j)
            m = mid(j)
            write(WF_ITM_UNIT, WF_QNTGAP, advance='no')    ! horizontal gap between channel and paraqn combination.
            do k = 1, WF_NITM
               write(WF_ITM_UNIT, '(1x)', advance='no')
               call wfat_print_qnttitle_fx(WF_ITM_UNIT, n, m, '(a15)')
            enddo
         enddo
      enddo
      write(WF_ITM_UNIT, '("")')

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_print_itmitmtitle(nch, nqn, wfat0)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer nch, nqn
      logical wfat0(nqn,nch)

      character*15 itm_title(WF_NITM)
      integer i, j, k

      if (.not. WF_ITM_PRINT_READY) 
     &     call errquit('WFAT: wfat_initprint_interm must first be '//
     &     'called before calling wfat_print_itmitmtitle.', 1, UERR)

      itm_title(1) = 'Re[a_nu]'
      itm_title(2) = 'Im[a_nu]'
      itm_title(3) = 'A_nu'
      itm_title(4) = 'tilde B_nu'
      itm_title(5) = 'B_nu'
      itm_title(6) = 'alpha_zz'
      itm_title(7) = 'o_nu'
      itm_title(8) = 'Re[g*.j]'
      itm_title(9) = 'Re[g*.k]'
      itm_title(10) = 'Re[g]'
      itm_title(11) = 'Im[g]'
      write(WF_ITM_UNIT, '("#", 2a15)', advance='no') ' ', ' '
      do i = 1, nch
         do j = 1, nqn
            if (wfat0(j,i)) cycle
            write(WF_ITM_UNIT, WF_QNTGAP, advance='no')    ! horizontal gap between channel and paraqn combination.
            do k = 1, WF_NITM
               write(WF_ITM_UNIT, '(1x)', advance='no')
               write(WF_ITM_UNIT, '(a15)', advance='no')
     &              adjustr(itm_title(k))
            enddo
         enddo
      enddo
      write(WF_ITM_UNIT, '("")')

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_print_itmline(nch, nqn, wfat0)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer nch, nqn
      logical wfat0(nqn,nch)

      integer i, j

      if (.not. WF_ITM_PRINT_READY) 
     &     call errquit('WFAT: wfat_initprint_interm must first be '//
     &     'called before calling wfat_print_itmline.', 1, UERR)

      write(WF_ITM_UNIT, '("#")', advance='no')
      call wfat_print_repeat(WF_ITM_UNIT, 2*15, 0, '-')
      do i = 1, nch
         do j = 1, nqn
            if (wfat0(j,i)) cycle
            call wfat_print_repeat(WF_ITM_UNIT, 3, 0, '-')   ! 3 is the length of WF_QNTGAP.
            call wfat_print_repeat(WF_ITM_UNIT, WF_NITM*(1+15), 0, '-')
         enddo
      enddo
      write(WF_ITM_UNIT, '("")')

      end
c============================================================
      

c=======================================================
c=======================================================
      subroutine wfat_endprint_interm()

      implicit none
#include "errquit.fh"
#include "wfat.fh"
#include "stdio.fh"      

      if (.not. WF_ITM_PRINT_READY) then
         call errquit('WFAT: wfat_initprint_interm must first '//
     &        'be called before calling wfat_endprint_interm.',
     &        1, UERR)
      else
         WF_ITM_PRINT_READY=.false.
      endif
      
      write(luout, '(/, 3x, "*** INFO ***")')      
      write(luout, '(5x, "WFAT(1) intermediate quantities have '//
     &     'been printed in ", a)')
     &     ''''//trim(adjustl(WF_ITM_FNAME))//'''.'
      
      close(WF_ITM_UNIT)
      
      end
c=======================================================

      
