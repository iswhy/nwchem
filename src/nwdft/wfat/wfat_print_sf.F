c=======================================================
c=======================================================
      subroutine wfat_initprint_structuref(nch, orbid0, spinid,
     &           nqn, nid, mid)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer nch, orbid0(nch), spinid(nch), nqn, nid(nqn), mid(nqn)

      character*6 wd_str
      integer dspn, dorb0, dnid, dmid, WF_SF_LINEWD


      if (WF_SF_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_structuref cannot be '//
     &     'called when the previous call to this routine has not '//
     &     'been terminated by calling wfat_endprint_structuref.',
     &     1, UERR)

      call util_file_name('wfats', .false., .false., WF_SF_FNAME)
      open(WF_SF_UNIT, file=WF_SF_FNAME, action='write')

      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Printing structure factor to ", a)')
     &     ''''//trim(adjustl(WF_SF_FNAME))//'''. Please wait ...'
      
c     ==== Column number and title ====
      WF_SF_LINEWD = 2*15 + nch*6 + nch*nqn*(3+3*15)
      write(WF_SF_UNIT, '(a)', advance='no') '#'
      call wfat_print_repeat(WF_SF_UNIT, WF_SF_LINEWD, 0, '-')
      write(WF_SF_UNIT, '(1x)')
      call wfat_print_sfcoln(WF_SF_UNIT, nch, nqn)
      call wfat_print_sftitle(WF_SF_UNIT, nch, orbid0, spinid, nqn,
     &     nid, mid, 'orb')
      call wfat_print_sftitle(WF_SF_UNIT, nch, orbid0, spinid, nqn,
     &     nid, mid, 'qnt')
      call wfat_print_sftitle(WF_SF_UNIT, nch, orbid0, spinid, nqn,
     &     nid, mid, 'par')
      write(WF_SF_UNIT, '(a)', advance='no') '#'
      call wfat_print_repeat(WF_SF_UNIT, WF_SF_LINEWD, 0, '-')
      write(WF_SF_UNIT, '(1x)')
      
      WF_SF_PRINT_READY = .true.
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_sfcoln(un, nch, nqn)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer un, nch, nqn

      character*10 k_str
      logical isopen
      integer k, kk, ich, iqn
      

      inquire(unit=un, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: un must first be opened '//
     &     'before calling wfat_print_sfcoln', 1, UERR)

      write(un, '(a)', advance='no') '#'
      write(un, '(a15)', advance='no') 'col. #1'
      write(un, '(a15)', advance='no') 'col. #2'
            
      k = 2
      do ich = 1, nch
         write(un, WF_ORBGAP, advance='no')
         do iqn = 1, nqn
            write(un, WF_QNTGAP, advance='no') 

            do kk = 1, 3
               k = k + 1
               write(k_str, '(i10)') k
               write(un, WF_SF_TITFMT, advance='no')
     &              'col. #'//trim(adjustl(k_str))
            enddo
            
         enddo
      enddo

      write(un, '(" ")') 

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_print_sftitle(un, nch, orbid0, spinid,
     &           nqn, nid, mid, title)

      implicit none
#include "wfat.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer un, nch, orbid0(nch), spinid(nch), nqn, nid(nqn), mid(nqn)
      character*3 title

      logical isopen
      character*6 io_str
      integer ich, iqn, io, is, n, m, i


      inquire(unit=un, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: un must be opened '//
     &     'before calling wfat_print_sftitle', 1, UERR)


      write(un, '(a)', advance='no') '#'
      call flush(un)

c     ==== beta space ====
      if (title == 'orb') then
         write(un, '(a15)', advance='no') 'beta (deg)'
      else
         write(un, '(a15)', advance='no') ' '
      endif
      call flush(un)

c     ==== gamma space ====      
      if (title == 'orb') then
         write(un, '(a15)', advance='no') 'gamma (deg)'
      else
         write(un, '(a15)', advance='no') ' '
      endif
      call flush(un)

c     ==== titles ====
      do ich = 1, nch
         write(un, WF_ORBGAP, advance='no')
         call flush(un)
         do iqn = 1, nqn
            write(un, WF_QNTGAP, advance='no')
            call flush(un)

            do i = 1, 3
               if (title == 'orb') then
                  io = orbid0(ich)
                  is = spinid(ich)
                  call wfat_print_orbtitle(un, is, io)
               else if (title == 'qnt') then
                  n = nid(iqn)
                  m = mid(iqn)
                  call wfat_print_qnttitle(un, n, m)
               else if (title == 'par') then
                  call wfat_print_partitle(un, i)
               else
                  call errquit('WFAT: Argument unit0 of '//
     &                 'wfat_print_sftitle corresponds to an '//
     &                 'unavailable option.', 2, UERR)
               endif
            enddo

         enddo
      enddo

      write(un, '(" ")') 
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_print_structuref(nat, bt, gm, nch, nqn, srfc, gap)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "geom.fh"      
      integer nat
      double precision bt, gm
      integer nch, nqn
      double complex srfc(1:*)
      logical gap

      integer i, ich, iqn
      double complex srfc_tot
      

      if (.not. WF_SF_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_structuref must first '//
     &     'be called before calling wfat_print_structuref.', 1, UERR)

      
      write(WF_SF_UNIT, '(a)', advance='no') ' '
      write(WF_SF_UNIT, '(2es15.5)', advance='no') bt, gm
      
      i = 1
      do ich = 1, nch
         write(WF_SF_UNIT, WF_ORBGAP, advance='no') 
         do iqn = 1, nqn
            write(WF_SF_UNIT, WF_QNTGAP, advance='no') 

            srfc_tot = sum(srfc(i : i+nat-1))
            write(WF_SF_UNIT, WF_SF_VALFMT, advance='no')
     &            srfc_tot, abs(srfc_tot)**2

            i = i + nat
         enddo
      enddo

      write(WF_SF_UNIT, '(" ")') 
      if (gap) write(WF_SF_UNIT, '(" ")') 

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_endprint_structuref()

      implicit none
#include "errquit.fh"
#include "wfat.fh"
#include "stdio.fh"
      

      if (.not. WF_SF_PRINT_READY) then
         call errquit('WFAT: wfat_initprint_structuref must first '//
     &        'be called before calling wfat_endprint_structuref.',
     &        1, UERR)
      else
         WF_SF_PRINT_READY=.false.
      endif

      write(luout, '(/, 3x, "*** INFO ***")')      
      write(luout, '(5x, "Structure factor has been printed in ", a)')
     &     ''''//trim(adjustl(WF_SF_FNAME))//'''.'
      
      close(WF_SF_UNIT)
      
      end
c=======================================================
