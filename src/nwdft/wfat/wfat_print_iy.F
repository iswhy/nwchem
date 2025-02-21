c=======================================================
c=======================================================
      subroutine wfat_initprint_yield()

      implicit none
#include "wfat.fh"      
      call wfat_initprint_yield0(WF_IY_UNIT)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_initprint_avbeta()

      implicit none
#include "wfat.fh"      
      call wfat_initprint_yield0(WF_ABTIY_UNIT)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_initprint_avgamma()

      implicit none
#include "wfat.fh"      
      call wfat_initprint_yield0(WF_AGMIY_UNIT)
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_initprint_yield0(iunit)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer iunit

      character*6 wd_str
      character*255 fname
      character*16 h
      integer dspn, dorb0, dnid, dmid, WF_IY_LINEWD, ang_wd


      call wfat_print_iy_ucheck(iunit)

      if ((iunit==WF_IY_UNIT .and. WF_IY_PRINT_READY) .or.
     &    (iunit==WF_ABTIY_UNIT .and. WF_ABTIY_PRINT_READY) .or.
     &    (iunit==WF_AGMIY_UNIT .and. WF_AGMIY_PRINT_READY))
     &     call errquit('WFAT: wfat_initprint_yield0 cannot be '//
     &     'called when the previous call to this routine has not '//
     &     'been terminated by calling wfat_endprint_yield.',
     &     1, UERR)
      

      if (iunit==WF_IY_UNIT) then
         call util_file_name('wfati', .false., .false., fname)
         WF_IY_FNAME = fname
         h = ' '
         ang_wd = 2
      else if (iunit==WF_ABTIY_UNIT) then
         call util_file_name('wfatab', .false., .false., fname)
         WF_ABTIY_FNAME = fname
         h = ' beta-averaged'
         ang_wd = 1
      else if (iunit==WF_AGMIY_UNIT) then
         call util_file_name('wfatag', .false., .false., fname)
         WF_AGMIY_FNAME = fname
         h = ' gamma-averaged'
         ang_wd = 1
      endif
      open(iunit, file=fname, action='write')
      
         
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Printing", a, " ion yield to ", a)') 
     &     trim(h),
     &     ''''//trim(adjustl(fname))//'''. Please wait ...'

      
c     ==== Column number and title ====
      WF_IY_LINEWD = ang_wd*15 + 6+14
      write(iunit, '(a)', advance='no') '#'
      call wfat_print_repeat(iunit, WF_IY_LINEWD, 0, '-')
      write(iunit, '(1x)')
      call wfat_print_iycoln(iunit)
      call wfat_print_iytitle(iunit)
      write(iunit, '(a)', advance='no') '#'
      call wfat_print_repeat(iunit, WF_IY_LINEWD, 0, '-')
      write(iunit, '(1x)')


      if (iunit==WF_IY_UNIT) then
         WF_IY_PRINT_READY = .true.
      else if (iunit==WF_ABTIY_UNIT) then
         WF_ABTIY_PRINT_READY = .true.
      else if (iunit==WF_AGMIY_UNIT) then
         WF_AGMIY_PRINT_READY = .true.
      endif
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_iycoln(iunit)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer iunit

      character*10 k_str
      logical isopen
      integer k
      

      call wfat_print_iy_ucheck(iunit)
      
      inquire(unit=iunit, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: iunit must first be opened '//
     &     'before calling wfat_print_iycoln', 1, UERR)

      write(iunit, '(a)', advance='no') '#'
      write(iunit, '(a15)', advance='no') 'col. #1'
      if (iunit==WF_IY_UNIT) then
         write(iunit, '(a15)', advance='no') 'col. #2'
         k = 2
      else
         k = 1
      endif
            
      k = k + 1
      write(iunit, WF_ORBGAP, advance='no')
      write(k_str, '(i10)') k
      write(iunit, WF_IY_TITFMT)
     &     'col. #'//trim(adjustl(k_str))

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_print_iytitle(iunit)

      implicit none
#include "wfat.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer iunit

      logical isopen


      call wfat_print_iy_ucheck(iunit)
      

      inquire(unit=iunit, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: iunit must be opened '//
     &     'before calling wfat_print_iytitle', 1, UERR)


      write(iunit, '(a)', advance='no') '#'
      call flush(iunit)

      
c     ==== beta space ====
      if (iunit /= WF_ABTIY_UNIT) then
         write(iunit, '(a15)', advance='no') 'beta (deg)'
         call flush(iunit)
      endif
      

c     ==== gamma space ====
      if (iunit /= WF_AGMIY_UNIT) then
         write(iunit, '(a15)', advance='no') 'gamma (deg)'
         call flush(iunit)
      endif
      

      write(iunit, WF_ORBGAP, advance='no')
      write(iunit, WF_IY_TITFMT) 'Total'
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_yield(bt, gm, ionyd, gap)

      implicit none
#include "wfat.fh"
      double precision bt, gm
      double precision ionyd
      logical gap
      call wfat_print_yield0(WF_IY_UNIT, bt, gm, ionyd, gap)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_avbeta(gm, ionyd)

      implicit none
#include "wfat.fh"
      double precision gm
      double precision ionyd
      call wfat_print_yield0(WF_ABTIY_UNIT, 0.0d0, gm, ionyd,
     &     .false.)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_avgamma(bt, ionyd)

      implicit none
#include "wfat.fh"
      double precision bt
      double precision ionyd
      call wfat_print_yield0(WF_AGMIY_UNIT, bt, 0.0d0, ionyd,
     &     .false.)
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_print_yield0(iunit, bt, gm, ionyd, gap)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "geom.fh"
#include "stdio.fh"      
      integer iunit
      double precision bt, gm
      double precision ionyd
      logical gap


      call wfat_print_iy_ucheck(iunit)

      if (iunit==WF_IY_UNIT) then
         if (.not. WF_IY_PRINT_READY)
     &        call errquit('WFAT: wfat_initprint_yield must first '//
     &        'be called before calling wfat_print_yield.', 1, UERR)
      else if (iunit==WF_ABTIY_UNIT) then
         if (.not. WF_ABTIY_PRINT_READY)
     &        call errquit('WFAT: wfat_initprint_avbeta must first '//
     &        'be called before calling wfat_print_avbeta.', 1, UERR)
      else if (iunit==WF_AGMIY_UNIT) then
         if (.not. WF_AGMIY_PRINT_READY)
     &        call errquit('WFAT: wfat_initprint_avgamma must first '//
     &        'be called before calling wfat_print_avgamma.', 1, UERR)
      endif
      
         
      write(iunit, '(a)', advance='no') ' '
      if (iunit /= WF_ABTIY_UNIT) 
     &     write(iunit, '(es15.6)', advance='no') bt
      if (iunit /= WF_AGMIY_UNIT) 
     &     write(iunit, '(es15.6)', advance='no') gm
         

      write(iunit, WF_ORBGAP, advance='no')
      write(iunit, WF_IY_VALFMT) ionyd
      if (gap) write(iunit, '(" ")') 

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_endprint_yield()

      implicit none
#include "wfat.fh"

      call wfat_endprint_yield0(WF_IY_UNIT)
      
      end
c=======================================================

      
c=======================================================      
c=======================================================
      subroutine wfat_endprint_avbeta()

      implicit none
#include "wfat.fh"

      call wfat_endprint_yield0(WF_ABTIY_UNIT)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_endprint_avgamma()

      implicit none
#include "wfat.fh"

      call wfat_endprint_yield0(WF_AGMIY_UNIT)
      
      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_endprint_yield0(iunit)

      implicit none
#include "errquit.fh"
#include "wfat.fh"
#include "stdio.fh"
      integer iunit

      character*255 fname
      character*16 h


      if (iunit==WF_IY_UNIT) then
         if (.not. WF_IY_PRINT_READY) then
            call errquit('WFAT: wfat_initprint_yield must first '//
     &           'be called before calling wfat_endprint_yield.',
     &           1, UERR)
         else
            WF_IY_PRINT_READY=.false.
         endif
         fname = WF_IY_FNAME
         h = ' '
      else if (iunit==WF_ABTIY_UNIT) then
         if (.not. WF_ABTIY_PRINT_READY) then
            call errquit('WFAT: wfat_initprint_avbeta must first '//
     &           'be called before calling wfat_endprint_avbeta.',
     &           1, UERR)
         else
            WF_ABTIY_PRINT_READY=.false.
         endif
         fname = WF_ABTIY_FNAME
         h = ' Beta-averaged'
      else if (iunit==WF_AGMIY_UNIT) then
         if (.not. WF_AGMIY_PRINT_READY) then
            call errquit('WFAT: wfat_initprint_avgamma must first '//
     &           'be called before calling wfat_endprint_avgamma.',
     &           1, UERR)
         else
            WF_AGMIY_PRINT_READY=.false.
         endif
         fname = WF_AGMIY_FNAME
         h = ' Gamma-averaged'
      endif

      
      write(luout, '(/, 3x, "*** INFO ***")')      
      write(luout, '(5x, a, " Ion yield has been printed in ", a)')
     &     trim(h), ''''//trim(adjustl(fname))//'''.' 
      
      close(iunit)
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_print_iy_ucheck(iunit)

      implicit none
#include "wfat.fh"
#include "errquit.fh"      
      integer iunit


      if (iunit /= WF_IY_UNIT .and. iunit /= WF_ABTIY_UNIT .and.
     &     iunit /= WF_AGMIY_UNIT) then
         call errquit('WFAT: The value of the argument iunit of '//
     &        'wfat_print_iy_ucheck is invalid.', 1, UERR)
      endif

      end
c=======================================================
