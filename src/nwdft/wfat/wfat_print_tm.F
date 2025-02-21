c=======================================================
c=======================================================
      subroutine wfat_initprint_tmd(nch, orbid0, spinid)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "stdio.fh"      
      integer nch, orbid0(nch), spinid(nch)

      character*6 wd_str
      integer dspn, dorb0, dnid, dmid, WF_TM_LINEWD
      integer k_len, val_len, orbgap_len


      if (WF_TM_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_tmd cannot be '//
     &     'called when the previous call to this routine has not '//
     &     'been terminated by calling wfat_endprint_tmd.',
     &     1, UERR)

      call util_file_name('wfatt', .false., .false., WF_TM_FNAME)
      open(WF_TM_UNIT, file=WF_TM_FNAME, action='write')

      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Printing TMD to ", a)')
     &     ''''//trim(adjustl(WF_TM_FNAME))//'''. Please wait ...'
      
c     ==== Column number and title ====
      k_len = 15
      val_len = 15
      orbgap_len = len(WF_ORBGAP)
      WF_TM_LINEWD = 2*k_len + nch*orbgap_len + nch*3*val_len
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
      write(WF_TM_UNIT, '(a15)', advance='no')
     &     'col. #1'
      write(WF_TM_UNIT, '(a15)', advance='no')
     &     'col. #2'
            
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
      double complex tmd(1:*)
      logical gap

      integer i, ich, natom
      double complex tmd_tot
      

      if (.not. WF_TM_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_tmd must first '//
     &     'be called before calling wfat_print_tmd.', 1, UERR)

      
c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom_hd, natom))
     &     call errquit('WFAT: Getting natom failed', 2, GEOM_ERR)

      write(WF_TM_UNIT, '(a)', advance='no') ' '
      write(WF_TM_UNIT, '(2es15.5)', advance='no') kx, ky
      
      i = 1
      do ich = 1, nch
         write(WF_TM_UNIT, WF_ORBGAP, advance='no') 

         tmd_tot = sum(tmd(i : i+natom-1))
         write(WF_TM_UNIT, WF_TM_VALFMT, advance='no')
     &         tmd_tot, abs(tmd_tot)**2

         i = i + natom
      enddo

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
      
