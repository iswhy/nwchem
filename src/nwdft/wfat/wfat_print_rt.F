c=======================================================
c=======================================================
      subroutine wfat_initprint_rate(nch, orbid0, spinid,
     &           nqn, nid, mid)

      implicit none
#include "wfat.fh"      
      integer nch, orbid0(nch), spinid(nch), nqn, nid(nqn), mid(nqn)
      

      call wfat_initprint_rate0(WF_RT_UNIT, nch, orbid0, spinid,
     &                          nqn, nid, mid)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_initprint_rate0(iunit, nch, orbid0, spinid,
     &           nqn, nid, mid)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer iunit
      integer nch, orbid0(nch), spinid(nch), nqn, nid(nqn), mid(nqn)

      character*6 wd_str
      character*255 fname
      character*16 h
      integer dspn, dorb0, dnid, dmid, WF_RT_LINEWD, ang_wd


      if (iunit==WF_RT_UNIT .and. WF_RT_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_rate0 cannot be '//
     &     'called when the previous call to this routine has not '//
     &     'been terminated by calling wfat_endprint_rate.',
     &     1, UERR)
      

      call util_file_name('wfatr', .false., .false., fname)
      WF_RT_FNAME = fname
      h = ' '
      ang_wd = 2
      open(iunit, file=fname, action='write')
      
         
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Printing", a, " ion rate to ", a)') 
     &     trim(h),
     &     ''''//trim(adjustl(fname))//'''. Please wait ...'

      
c     ==== Column number and title ====
      WF_RT_LINEWD = ang_wd*15 + nch*6 + nch*nqn*(3+14) + 6+14
      write(iunit, '(a)', advance='no') '#'
      call wfat_print_repeat(iunit, WF_RT_LINEWD, 0, '-')
      write(iunit, '(1x)')
      call wfat_print_rtcoln(iunit, nch, nqn)
      call wfat_print_rttitle(iunit, nch, orbid0, spinid, nqn, nid, mid,
     &     'orb')
      call wfat_print_rttitle(iunit, nch, orbid0, spinid, nqn, nid, mid,
     &     'qnt')
      write(iunit, '(a)', advance='no') '#'
      call wfat_print_repeat(iunit, WF_RT_LINEWD, 0, '-')
      write(iunit, '(1x)')


      WF_RT_PRINT_READY = .true.      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_rtcoln(iunit, nch, nqn)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer iunit, nch, nqn

      character*10 k_str
      logical isopen
      integer k, kk, ich, iqn
      

      inquire(unit=iunit, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: iunit must first be opened '//
     &     'before calling wfat_print_rtcoln', 1, UERR)

      write(iunit, '(a)', advance='no') '#'
      write(iunit, '(a15)', advance='no') 'col. #1'
      write(iunit, '(a15)', advance='no') 'col. #2'
      k = 2
      
            
      do ich = 1, nch
         write(iunit, WF_ORBGAP, advance='no')
         do iqn = 1, nqn
            write(iunit, WF_QNTGAP, advance='no') 

            k = k + 1
            write(k_str, '(i10)') k
            write(iunit, WF_RT_TITFMT, advance='no')
     &           'col. #'//trim(adjustl(k_str))            
         enddo
      enddo

      k = k + 1
      write(iunit, WF_ORBGAP, advance='no')
      write(k_str, '(i10)') k
      write(iunit, WF_RT_TITFMT)
     &     'col. #'//trim(adjustl(k_str))

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_print_rttitle(iunit, nch, orbid0, spinid,
     &           nqn, nid, mid, title)

      implicit none
#include "wfat.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer iunit
      integer nch, orbid0(nch), spinid(nch), nqn, nid(nqn), mid(nqn)
      character*3 title

      logical isopen
      character*6 io_str
      integer ich, iqn, io, is, n, m, i


      inquire(unit=iunit, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: iunit must be opened '//
     &     'before calling wfat_print_rttitle', 1, UERR)


      write(iunit, '(a)', advance='no') '#'
      call flush(iunit)

      
c     ==== beta space ====
      if (title == 'orb') then
         write(iunit, '(a15)', advance='no') 'beta (deg)'
      else
         write(iunit, '(a15)', advance='no') ' '
      endif
      call flush(iunit)
      

c     ==== gamma space ====
      if (title == 'orb') then
         write(iunit, '(a15)', advance='no') 'gamma (deg)'
      else
         write(iunit, '(a15)', advance='no') ' '
      endif
      call flush(iunit)
      

c     ==== titles ====
      do ich = 1, nch
         write(iunit, WF_ORBGAP, advance='no')
         call flush(iunit)
         do iqn = 1, nqn
            write(iunit, WF_QNTGAP, advance='no')
            call flush(iunit)

            if (title == 'orb') then
               io = orbid0(ich)
               is = spinid(ich)
               call wfat_print_orbtitle(iunit, is, io)
            else if (title == 'qnt') then
               n = nid(iqn)
               m = mid(iqn)
               call wfat_print_qnttitle(iunit, n, m)
            else
               call errquit('WFAT: Argument unit0 of '//
     &              'wfat_print_rttitle corresponds to an '//
     &              'unavailable option.', 2, UERR)
            endif

         enddo
      enddo

      write(iunit, WF_ORBGAP, advance='no')
      if (title == 'orb') then
         write(iunit, WF_RT_TITFMT) 'Total'
      else if (title == 'qnt') then
         write(iunit, '(" ")')
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_rate(bt, gm, nch, nqn, p_ionyd, 
     &                           ionyd, gap)

      implicit none
#include "wfat.fh"
      integer nch, nqn
      double precision bt, gm
      double precision p_ionyd(1:*), ionyd
      logical gap
      

      call wfat_print_rate0(WF_RT_UNIT, bt, gm, nch, nqn, p_ionyd, 
     &                      ionyd, gap)
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_print_rate0(iunit, bt, gm, nch, nqn, p_ionyd, 
     &                            ionyd, gap)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "geom.fh"
#include "stdio.fh"      
      integer iunit
      integer nch, nqn
      double precision bt, gm
      double precision p_ionyd(1:*), ionyd
      logical gap

      integer i, ich, iqn


      if (.not. WF_RT_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_rate must first '//
     &     'be called before calling wfat_print_rate.', 1, UERR)
      
         
      write(iunit, '(a)', advance='no') ' '
      write(iunit, '(es15.6)', advance='no') bt
      write(iunit, '(es15.6)', advance='no') gm
         

      i = 0
      do ich = 1, nch
         write(iunit, WF_ORBGAP, advance='no') 
         do iqn = 1, nqn
            write(iunit, WF_QNTGAP, advance='no') 

            i = i + 1
            write(iunit, WF_RT_VALFMT, advance='no') p_ionyd(i)
         enddo
      enddo

      write(iunit, WF_ORBGAP, advance='no')
      write(iunit, WF_RT_VALFMT) ionyd
      if (gap) write(iunit, '(" ")') 

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_endprint_rate()

      implicit none
#include "wfat.fh"

      call wfat_endprint_rate0(WF_RT_UNIT)
      
      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_endprint_rate0(iunit)

      implicit none
#include "errquit.fh"
#include "wfat.fh"
#include "stdio.fh"
      integer iunit

      character*255 fname


      if (.not. WF_RT_PRINT_READY) then
         call errquit('WFAT: wfat_initprint_rate must first '//
     &        'be called before calling wfat_endprint_rate.',
     &        1, UERR)
      else
         WF_RT_PRINT_READY=.false.
      endif
      fname = WF_RT_FNAME

      
      write(luout, '(/, 3x, "*** INFO ***")')      
      write(luout, '(5x, "Ionization rate has been printed in ", a)')
     &     ''''//trim(adjustl(fname))//'''.' 
      
      close(iunit)
      
      end
c=======================================================
      
