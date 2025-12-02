c=======================================================
c=======================================================
      subroutine wfat_initprint_tmd(bt, gm, nch, orbid0, spinid)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "stdio.fh"
      double precision bt, gm
      integer nch, orbid0(nch), spinid(nch)

      character*1 tt, tl
      character*6 wd_str
      integer dspn, dorb0, dnid, dmid, WF_TM_LINEWD
      integer k_len, val_len


      if (WF_TM_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_tmd cannot be '//
     &     'called when the previous call to this routine has not '//
     &     'been terminated by calling wfat_endprint_tmd.',
     &     1, UERR)

      tmdfile = tmdfile + 1
      if (tmdfile > tmdlim) then
         write(tl, '(i1)') tmdlim
         call errquit('WFAT: The number of TMD print files exceeded '//
     &        'its limit, which is '//tl//'.')
      endif
      write(tt, '(i1)') tmdfile
      call util_file_name('wfatt.'//tt, .false., .false., WF_TM_FNAME)
      open(WF_TM_UNIT, file=WF_TM_FNAME, action='write')

      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Printing TMD at beta = ", es15.5,'//
     &     '" and gamma = ", es15.5, " to ", a)') bt, gm,
     &     ''''//trim(adjustl(WF_TM_FNAME))//'''. Please wait ...'

      write(WF_TM_UNIT, '("# beta = ", es15.5, 5x, "gamma = ", es15.5)')
     &     bt, gm
      
c     ==== Column number and title ====
      k_len = 15
      val_len = 15
cOLD      orbgap_len = len(WF_ORBGAP)
      WF_TM_LINEWD =
     &     2*k_len +              ! total length of kx and ky
     &     (nch+1)*6 +            ! total length of space between orbitals
     &     nch*3*val_len +        ! total length of channel TMDs' real, imag, and abs.
     &     val_len                ! length of total TMD
      write(WF_TM_UNIT, '(a)', advance='no') '#'
      call wfat_print_repeat(WF_TM_UNIT, WF_TM_LINEWD, 0, '-')
      write(WF_TM_UNIT, '(1x)')
      call wfat_print_tmcoln(nch)
      call wfat_print_tmtitle(nch, orbid0, spinid, 'orb')
      call wfat_print_tmtitle(nch, orbid0, spinid, 'par')
      write(WF_TM_UNIT, '(a)', advance='no') '#'
      call wfat_print_repeat(WF_TM_UNIT, WF_TM_LINEWD, 0, '-')
      write(WF_TM_UNIT, '(1x)')
      
      WF_TM_PRINT_READY = .true.

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_tmcoln(nch)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer nch, nqn

      character*10 k_str
      logical isopen
      integer k, kk, ich
      

      inquire(unit=WF_TM_UNIT, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: WF_TM_UNIT must first be opened '//
     &     'before calling wfat_print_tmcoln', 1, UERR)

      write(WF_TM_UNIT, '(a)', advance='no') '#'
      write(WF_TM_UNIT, '(a15)', advance='no') 'col. #1'
      write(WF_TM_UNIT, '(a15)', advance='no') 'col. #2'

c     ==== Column numbers for ionization channel TMDs ====
      k = 2
      do ich = 1, nch
         write(WF_TM_UNIT, WF_ORBGAP, advance='no')
         do kk = 1, 3
            k = k + 1
            write(k_str, '(i10)') k
            write(WF_TM_UNIT, WF_TM_TITFMT, advance='no')
     &           'col. #'//trim(adjustl(k_str))
         enddo
      enddo

c     ==== Column number for the total TMD ====
      write(WF_TM_UNIT, WF_ORBGAP, advance='no')
      write(k_str, '(i10)') k + 1
      write(WF_TM_UNIT, WF_TM_TITFMT, advance='no')
     &           'col. #'//trim(adjustl(k_str))

      write(WF_TM_UNIT, '(" ")') 

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_print_tmtitle(nch, orbid0, spinid, title)

      implicit none
#include "wfat.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer nch, orbid0(nch), spinid(nch)
      character*3 title

      logical isopen
      character*6 io_str
      integer ich, io, is, i

      
      inquire(unit=WF_TM_UNIT, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: WF_TM_UNIT must be opened '//
     &     'before calling wfat_print_tmtitle', 1, UERR)


      write(WF_TM_UNIT, '(a)', advance='no') '#'
      call flush(WF_TM_UNIT)

c     ==== kx space ====
      if (title == 'orb') then
         write(WF_TM_UNIT, '(a15)', advance='no') 'kx (a.u.)'
      else
         write(WF_TM_UNIT, '(a15)', advance='no') ' '
      endif
      call flush(WF_TM_UNIT)

c     ==== ky space ====
      if (title == 'orb') then
         write(WF_TM_UNIT, '(a15)', advance='no') 'ky (a.u.)'
      else
         write(WF_TM_UNIT, '(a15)', advance='no') ' '
      endif
      call flush(WF_TM_UNIT)

c     ==== titles ====
      do ich = 1, nch
         write(WF_TM_UNIT, WF_ORBGAP, advance='no')
         call flush(WF_TM_UNIT)
         
         do i = 1, 3
            if (title == 'orb') then
               io = orbid0(ich)
               is = spinid(ich)
               call wfat_print_orbtitle(WF_TM_UNIT, is, io)
            else if (title == 'par') then
               call wfat_print_partitle(WF_TM_UNIT, i)
            else
               call errquit('WFAT: Argument unit0 of '//
     &              'wfat_print_tmtitle corresponds to an '//
     &              'unavailable option.', 2, UERR)
            endif
         enddo

      enddo


      write(WF_TM_UNIT, WF_ORBGAP, advance='no')
      if (title == 'orb') then
         write(WF_TM_UNIT, WF_TM_TITFMT, advance='no') 'Total'
      else if (title == 'par') then
         write(WF_TM_UNIT, WF_TM_TITFMT, advance='no') ' '
      else
         call errquit('WFAT: Argument unit0 of wfat_print_tmtitle '//
     &        'corresponds to an unavailable option.', 2, UERR)
      endif

      write(WF_TM_UNIT, '(" ")') 
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_print_tmd(kx, ky, nch, tmd, gap)

      implicit none
#include "wfat.fh"
#include "wfat_qm_param.fh"      
c#include "cdft.fh"
#include "errquit.fh"
#include "geom.fh"      
      integer nch
      double precision kx, ky
      double complex tmd(nch)
      logical gap

      integer ich

      
      if (.not. WF_TM_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_tmd must first '//
     &     'be called before calling wfat_print_tmd.', 1, UERR)

      
      write(WF_TM_UNIT, '(a)', advance='no') ' '
      write(WF_TM_UNIT, '(2es15.5)', advance='no') kx, ky
      
      do ich = 1, nch
         write(WF_TM_UNIT, WF_ORBGAP, advance='no')
         write(WF_TM_UNIT, WF_TM_VALFMT, advance='no')
     &         tmd(ich), abs(tmd(ich))
      enddo

      write(WF_TM_UNIT, WF_ORBGAP, advance='no')
      write(WF_TM_UNIT, '(es15.6)', advance='no')
     &         sum(abs(tmd) * abs(tmd))

      write(WF_TM_UNIT, '(" ")') 
      if (gap) write(WF_TM_UNIT, '(" ")') 
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_endprint_tmd()

      implicit none
#include "errquit.fh"
#include "wfat.fh"
#include "stdio.fh"
      

      if (.not. WF_TM_PRINT_READY) then
         call errquit('WFAT: wfat_initprint_tmd must first '//
     &        'be called before calling wfat_endprint_tmd.',
     &        1, UERR)
      else
         WF_TM_PRINT_READY=.false.
      endif

      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "TMD has been printed in ", a)')
     &     ''''//trim(adjustl(WF_TM_FNAME))//'''.'
      
      close(WF_TM_UNIT)
      
      end
c=======================================================
