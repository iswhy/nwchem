c=======================================================
c=======================================================
      subroutine wfat_initprint_srfc1(nch, orbid0, spinid,
     &           nqn, nid, mid)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer nch, orbid0(nch), spinid(nch), nqn, nid(nqn), mid(nqn)

      character*6 wd_str
      character*5 wf1_suff(2)
      integer i, WF_SF_LINEWD


      if (WF_SF1_PRINT_READY(1) .or. WF_SF1_PRINT_READY(2))
     &   call errquit('WFAT: wfat_initprint_srfc1 cannot be '//
     &   'called when the previous call to this routine has not '//
     &   'been terminated by calling wfat_endprint_srfc1.',
     &   1, UERR)

      wf1_suff = (/ 'wfatj', 'wfatk' /)
      do i = 1, 2
         call util_file_name(wf1_suff(i), .false., .false.,
     &        WF_SF1_FNAME(i))
         open(WF_SF1_UNIT(i), file=WF_SF1_FNAME(i), action='write')
      enddo
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout,
     &     '(5x, "Printing J_n to ", a, " and K_n to ", a, a)')
     &     ''''//trim(adjustl(WF_SF1_FNAME(1)))//'''',
     &     ''''//trim(adjustl(WF_SF1_FNAME(2)))//'''',
     &     ' Please wait ...'

         
      do i = 1, 2
c        ==== Column number and title ====
         WF_SF_LINEWD = 2*15 + nch*6 + nch*nqn*(3+3*15)
         write(WF_SF1_UNIT(i), '(a)', advance='no') '#'
         call wfat_print_repeat(WF_SF1_UNIT(i), WF_SF_LINEWD, 0, '-')
         write(WF_SF1_UNIT(i), '(1x)')
         call wfat_print_sfcoln(WF_SF1_UNIT(i), nch, nqn)
         call wfat_print_sftitle(WF_SF1_UNIT(i), nch, orbid0, spinid, 
     &        nqn, nid, mid, 'orb')
         call wfat_print_sftitle(WF_SF1_UNIT(i), nch, orbid0, spinid, 
     &        nqn, nid, mid, 'qnt')
         call wfat_print_sftitle(WF_SF1_UNIT(i), nch, orbid0, spinid,
     &        nqn, nid, mid, 'par')
         write(WF_SF1_UNIT(i), '(a)', advance='no') '#'
         call wfat_print_repeat(WF_SF1_UNIT(i), WF_SF_LINEWD, 0, '-')
         write(WF_SF1_UNIT(i), '(1x)')
         
         WF_SF1_PRINT_READY(i) = .true.
      enddo
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_srfc1(nat, bt, gm, nch, nqn, jn, kn,
     &     gap)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "geom.fh"      
      integer nat
      double precision bt, gm
      integer nch, nqn
      double complex jn(nat, nqn, nch)
      double complex kn(nat, nqn, nch)
      logical gap

      integer j, ich, iqn
      double complex jn_tot, kn_tot
      

      if ( .not. (WF_SF1_PRINT_READY(1) .and. WF_SF1_PRINT_READY(2)) )
     &     call errquit('WFAT: wfat_initprint_srfc1 must first '//
     &     'be called before calling wfat_print_srfc1.', 1, UERR)


      do j = 1, 2
         write(WF_SF1_UNIT(j), '(a)', advance='no') ' '
         write(WF_SF1_UNIT(j), '(2es15.5)', advance='no') bt, gm
      enddo
      

      do ich = 1, nch
         write(WF_SF1_UNIT(1), WF_ORBGAP, advance='no')
         write(WF_SF1_UNIT(2), WF_ORBGAP, advance='no') 
         do iqn = 1, nqn
            write(WF_SF1_UNIT(1), WF_QNTGAP, advance='no')
            write(WF_SF1_UNIT(2), WF_QNTGAP, advance='no') 

            jn_tot = sum(jn(:, iqn, ich))
            write(WF_SF1_UNIT(1), WF_SF_VALFMT, advance='no')
     &           jn_tot, abs(jn_tot)**2
            kn_tot = sum(kn(:, iqn, ich))
            write(WF_SF1_UNIT(2), WF_SF_VALFMT, advance='no')
     &           kn_tot, abs(kn_tot)**2
         enddo
      enddo

      do j = 1, 2
         write(WF_SF1_UNIT(j), '(" ")') 
         if (gap) write(WF_SF1_UNIT(j), '(" ")')
      enddo

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_endprint_srfc1()

      implicit none
#include "errquit.fh"
#include "wfat.fh"
#include "stdio.fh"

      integer i
      

      do i = 1, 2
         if (.not. WF_SF1_PRINT_READY(i)) then
            call errquit('WFAT: wfat_initprint_srfc1 must first '//
     &           'be called before calling wfat_endprint_srfc1.',
     &           1, UERR)
         else
            WF_SF1_PRINT_READY(i)=.false.
         endif
         
         write(luout, '(/, 3x, "*** INFO ***")')      
         write(luout, '(5x, "WFAT(1) Structure factor has been '//
     &        'printed in ", a)')
     &        ''''//trim(adjustl(WF_SF1_FNAME(i)))//'''.'
         
         close(WF_SF1_UNIT(i))
      enddo
      
      end
c=======================================================      
