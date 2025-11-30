cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_initprint_tmd(nch, orbid0, spinid)
cOLD
cOLD      implicit none
cOLD#include "wfat.fh"
cOLD#include "errquit.fh"
cOLD#include "stdio.fh"      
cOLD      integer nch, orbid0(nch), spinid(nch)
cOLD
cOLD      character*6 wd_str
cOLD      integer dspn, dorb0, dnid, dmid, WF_TM_LINEWD
cOLD      integer k_len, val_len, orbgap_len
cOLD
cOLD
cOLD      if (WF_TM_PRINT_READY)
cOLD     &     call errquit('WFAT: wfat_initprint_tmd cannot be '//
cOLD     &     'called when the previous call to this routine has not '//
cOLD     &     'been terminated by calling wfat_endprint_tmd.',
cOLD     &     1, UERR)
cOLD
cOLD      call util_file_name('wfatt', .false., .false., WF_TM_FNAME)
cOLD      open(WF_TM_UNIT, file=WF_TM_FNAME, action='write')
cOLD
cOLD      write(luout, '(/, 3x, "*** INFO ***")')
cOLD      write(luout, '(5x, "Printing TMD to ", a)')
cOLD     &     ''''//trim(adjustl(WF_TM_FNAME))//'''. Please wait ...'
cOLD      
cOLDc     ==== Column number and title ====
cOLD      k_len = 15
cOLD      val_len = 15
cOLD      orbgap_len = len(WF_ORBGAP)
cOLD      WF_TM_LINEWD = 2*k_len + nch*orbgap_len + nch*3*val_len
cOLD      write(WF_TM_UNIT, '(a)', advance='no') '#'
cOLD      call wfat_print_repeat(WF_TM_UNIT, WF_TM_LINEWD, 0, '-')
cOLD      write(WF_TM_UNIT, '(1x)')
cOLD      call wfat_print_tmcoln(nch)
cOLD      call wfat_print_tmtitle(nch, orbid0, spinid, 'orb')
cOLD      call wfat_print_tmtitle(nch, orbid0, spinid, 'par')
cOLD      write(WF_TM_UNIT, '(a)', advance='no') '#'
cOLD      call wfat_print_repeat(WF_TM_UNIT, WF_TM_LINEWD, 0, '-')
cOLD      write(WF_TM_UNIT, '(1x)')
cOLD      
cOLD      WF_TM_PRINT_READY = .true.
cOLD      
cOLD      end
cOLDc=======================================================


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_print_tmcoln(nch)
cOLD
cOLD      implicit none
cOLD#include "wfat.fh"
cOLD#include "errquit.fh"
cOLD      integer nch, nqn
cOLD
cOLD      character*10 k_str
cOLD      logical isopen
cOLD      integer k, kk, ich
cOLD      
cOLD
cOLD      inquire(unit=WF_TM_UNIT, opened=isopen)
cOLD      if (.not. isopen)
cOLD     &     call errquit('WFAT: WF_TM_UNIT must first be opened '//
cOLD     &     'before calling wfat_print_tmcoln', 1, UERR)
cOLD
cOLD      write(WF_TM_UNIT, '(a)', advance='no') '#'
cOLD      write(WF_TM_UNIT, '(a15)', advance='no')
cOLD     &     'col. #1'
cOLD      write(WF_TM_UNIT, '(a15)', advance='no')
cOLD     &     'col. #2'
cOLD            
cOLD      k = 2
cOLD      do ich = 1, nch
cOLD         write(WF_TM_UNIT, WF_ORBGAP, advance='no')
cOLD         do kk = 1, 3
cOLD            k = k + 1
cOLD            write(k_str, '(i10)') k
cOLD            write(WF_TM_UNIT, WF_TM_TITFMT, advance='no')
cOLD     &           'col. #'//trim(adjustl(k_str))
cOLD         enddo
cOLD      enddo
cOLD
cOLD      write(WF_TM_UNIT, '(" ")') 
cOLD
cOLD      end
cOLDc=======================================================


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_print_tmtitle(nch, orbid0, spinid, title)
cOLD
cOLD      implicit none
cOLD#include "wfat.fh"
cOLD#include "stdio.fh"
cOLD#include "errquit.fh"
cOLD      integer nch, orbid0(nch), spinid(nch)
cOLD      character*3 title
cOLD
cOLD      logical isopen
cOLD      character*6 io_str
cOLD      integer ich, io, is, i
cOLD
cOLD
cOLD      inquire(unit=WF_TM_UNIT, opened=isopen)
cOLD      if (.not. isopen)
cOLD     &     call errquit('WFAT: WF_TM_UNIT must be opened '//
cOLD     &     'before calling wfat_print_tmtitle', 1, UERR)
cOLD
cOLD
cOLD      write(WF_TM_UNIT, '(a)', advance='no') '#'
cOLD      call flush(WF_TM_UNIT)
cOLD
cOLDc     ==== kx space ====
cOLD      if (title == 'orb') then
cOLD         write(WF_TM_UNIT, '(a15)', advance='no') 'kx (a.u.)'
cOLD      else
cOLD         write(WF_TM_UNIT, '(a15)', advance='no') ' '
cOLD      endif
cOLD      call flush(WF_TM_UNIT)
cOLD
cOLDc     ==== ky space ====      
cOLD      if (title == 'orb') then
cOLD         write(WF_TM_UNIT, '(a15)', advance='no') 'ky (a.u.)'
cOLD      else
cOLD         write(WF_TM_UNIT, '(a15)', advance='no') ' '
cOLD      endif
cOLD      call flush(WF_TM_UNIT)
cOLD
cOLDc     ==== titles ====
cOLD      do ich = 1, nch
cOLD         write(WF_TM_UNIT, WF_ORBGAP, advance='no')
cOLD         call flush(WF_TM_UNIT)
cOLD         
cOLD         do i = 1, 3
cOLD            if (title == 'orb') then
cOLD               io = orbid0(ich)
cOLD               is = spinid(ich)
cOLD               call wfat_print_orbtitle(WF_TM_UNIT, is, io)
cOLD            else if (title == 'par') then
cOLD               call wfat_print_partitle(WF_TM_UNIT, i)
cOLD            else
cOLD               call errquit('WFAT: Argument unit0 of '//
cOLD     &              'wfat_print_tmtitle corresponds to an '//
cOLD     &              'unavailable option.', 2, UERR)
cOLD            endif
cOLD         enddo
cOLD
cOLD      enddo
cOLD
cOLD      write(WF_TM_UNIT, '(" ")') 
cOLD      
cOLD      end
cOLDc=======================================================


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_print_tmd(kx, ky, nch, tmd, gap)
cOLD
cOLD      implicit none
cOLD#include "wfat.fh"
cOLD#include "wfat_qm_param.fh"      
cOLDc#include "cdft.fh"
cOLD#include "errquit.fh"
cOLD#include "geom.fh"      
cOLD      integer nch
cOLD      double precision kx, ky
cOLD      double complex tmd(1:*)
cOLD      logical gap
cOLD
cOLD      integer i, ich, natom
cOLD      double complex tmd_tot
cOLD      
cOLD
cOLD      if (.not. WF_TM_PRINT_READY)
cOLD     &     call errquit('WFAT: wfat_initprint_tmd must first '//
cOLD     &     'be called before calling wfat_print_tmd.', 1, UERR)
cOLD
cOLD      
cOLDc     ==== Get number of atoms ====
cOLD      if (.not. geom_ncent(geom_hd, natom))
cOLD     &     call errquit('WFAT: Getting natom failed', 2, GEOM_ERR)
cOLD
cOLD      write(WF_TM_UNIT, '(a)', advance='no') ' '
cOLD      write(WF_TM_UNIT, '(2es15.5)', advance='no') kx, ky
cOLD      
cOLD      i = 1
cOLD      do ich = 1, nch
cOLD         write(WF_TM_UNIT, WF_ORBGAP, advance='no') 
cOLD
cOLD         tmd_tot = sum(tmd(i : i+natom-1))
cOLD         write(WF_TM_UNIT, WF_TM_VALFMT, advance='no')
cOLD     &         tmd_tot, abs(tmd_tot)**2
cOLD
cOLD         i = i + natom
cOLD      enddo
cOLD
cOLD      write(WF_TM_UNIT, '(" ")') 
cOLD      if (gap) write(WF_TM_UNIT, '(" ")') 
cOLD
cOLD      end
cOLDc=======================================================      


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_endprint_tmd()
cOLD
cOLD      implicit none
cOLD#include "errquit.fh"
cOLD#include "wfat.fh"
cOLD#include "stdio.fh"
cOLD      
cOLD
cOLD      if (.not. WF_TM_PRINT_READY) then
cOLD         call errquit('WFAT: wfat_initprint_tmd must first '//
cOLD     &        'be called before calling wfat_endprint_tmd.',
cOLD     &        1, UERR)
cOLD      else
cOLD         WF_TM_PRINT_READY=.false.
cOLD      endif
cOLD
cOLD      write(luout, '(/, 3x, "*** INFO ***")')      
cOLD      write(luout, '(5x, "TMD has been printed in ", a)')
cOLD     &     ''''//trim(adjustl(WF_TM_FNAME))//'''.'
cOLD      
cOLD      close(WF_TM_UNIT)
cOLD      
cOLD      end
cOLDc=======================================================
      
