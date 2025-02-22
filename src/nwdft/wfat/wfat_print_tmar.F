c=======================================================
c=======================================================
      subroutine wfat_initprint_tmdar(nch, natom, orbid0, spinid)

      implicit none
#include "wfat.fh"
#include "stdio.fh"      
#include "errquit.fh"
      integer nch, natom, orbid0(nch), spinid(nch)

      character*6 wd_str
      integer dspn, dorb0, dnid, dmid, WF_TMAR_LINEWD
      integer k_len, val_len, orbgap_len, atmgap_len


      if (WF_TMAR_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_tmdar cannot be '//
     &     'called when the previous call to this routine has not '//
     &     'been terminated by calling wfat_endprint_tmdar.',
     &     1, UERR)

      call util_file_name('wfatt.ar', .false., .false., WF_TMAR_FNAME)
      open(WF_TMAR_UNIT, file=WF_TMAR_FNAME, action='write')

      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Printing atom-resolved TMD to ", a)')
     &     ''''//trim(adjustl(WF_TMAR_FNAME))//'''. Please wait ...'
      
c     ==== Column number and title ====
      k_len = 15
      val_len = 15
      atmgap_len = 3
      orbgap_len = 6
      WF_TMAR_LINEWD = 2*k_len + nch*orbgap_len +
     &                 nch*natom*(atmgap_len + 3*val_len)
      write(WF_TMAR_UNIT, '(a)', advance='no') '#'
      call wfat_print_repeat(WF_TMAR_UNIT, WF_TMAR_LINEWD, 0, '-')
      write(WF_TMAR_UNIT, '(1x)')
      call wfat_print_tmarcoln(nch, natom)
      call wfat_print_tmartitle(nch, natom, orbid0, spinid, 'orb')
      call wfat_print_tmartitle(nch, natom, orbid0, spinid, 'atm')
      call wfat_print_tmartitle(nch, natom, orbid0, spinid, 'par')
      write(WF_TMAR_UNIT, '(a)', advance='no') '#'
      call wfat_print_repeat(WF_TMAR_UNIT, WF_TMAR_LINEWD, 0, '-')
      write(WF_TMAR_UNIT, '(1x)')
      
      WF_TMAR_PRINT_READY = .true.
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_tmarcoln(nch, natom)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
      integer nch, natom

      character*10 k_str
      logical isopen
      integer k, kk, ich, iat
      

      inquire(unit=WF_TMAR_UNIT, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: WF_TMAR_UNIT must first be opened '//
     &     'before calling wfat_print_tmarcoln', 1, UERR)

      write(WF_TMAR_UNIT, '(a)', advance='no') '#'
      write(WF_TMAR_UNIT, '(a15)', advance='no') 'col. #1'
      write(WF_TMAR_UNIT, '(a15)', advance='no') 'col. #2'
            
      k = 2
      do ich = 1, nch
         write(WF_TMAR_UNIT, WF_ORBGAP, advance='no')
         do iat = 1, natom
            write(WF_TMAR_UNIT, WF_ATMGAP, advance='no') 

            do kk = 1, 3
               k = k + 1
               write(k_str, '(i10)') k
               write(WF_TMAR_UNIT, WF_TMAR_TITFMT, advance='no')
     &              'col. #'//trim(adjustl(k_str))
            enddo
            
         enddo
      enddo

      write(WF_TMAR_UNIT, '(" ")') 

      end
c=======================================================      


c=======================================================
c=======================================================
      subroutine wfat_print_tmartitle(nch, natom, orbid0, spinid,
     &           title)

      implicit none
#include "wfat.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer nch, orbid0(nch), spinid(nch), natom
      character*3 title

      logical isopen
      character*6 io_str
      integer ich, iat, io, is, i


      inquire(unit=WF_TMAR_UNIT, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: WF_TMAR_UNIT must be opened '//
     &     'before calling wfat_print_tmartitle', 1, UERR)


      write(WF_TMAR_UNIT, '(a)', advance='no') '#'
      call flush(WF_TMAR_UNIT)

c     ==== kx space ====
      if (title == 'orb') then
         write(WF_TMAR_UNIT, '(a15)', advance='no') 'kx (a.u.)'
      else
         write(WF_TMAR_UNIT, '(a15)', advance='no') ' '
      endif
      call flush(WF_TMAR_UNIT)

c     ==== ky space ====      
      if (title == 'orb') then
         write(WF_TMAR_UNIT, '(a15)', advance='no') 'ky (a.u.)'
      else
         write(WF_TMAR_UNIT, '(a15)', advance='no') ' '
      endif
      call flush(WF_TMAR_UNIT)

c     ==== titles ====
      do ich = 1, nch
         write(WF_TMAR_UNIT, WF_ORBGAP, advance='no')
         call flush(WF_TMAR_UNIT)
         do iat = 1, natom
            write(WF_TMAR_UNIT, WF_ATMGAP, advance='no')
            call flush(WF_TMAR_UNIT)

            do i = 1, 3
               if (title == 'orb') then
                  io = orbid0(ich)
                  is = spinid(ich)
                  call wfat_print_orbtitle(WF_TMAR_UNIT, is, io)
               else if (title == 'atm') then
                  call wfat_print_atmtitle(WF_TMAR_UNIT, iat)
               else if (title == 'par') then
                  call wfat_print_partitle(WF_TMAR_UNIT, i)
               else
                  call errquit('WFAT: Argument unit0 of '//
     &                 'wfat_print_tmartitle corresponds to an '//
     &                 'unavailable option.', 2, UERR)
               endif
            enddo

         enddo
      enddo

      write(WF_TMAR_UNIT, '(" ")') 
      
      end
c=======================================================      


c=======================================================
c=======================================================
      subroutine wfat_print_tmdar(kx, ky, nch, natom, tmd, gap)

      implicit none
#include "wfat.fh"
#include "wfat_qm_param.fh"      
c#include "cdft.fh"
#include "errquit.fh"
#include "geom.fh"      
      integer nch, natom
      double precision kx, ky
      double complex tmd(natom, nch)
      logical gap

      integer i, ich, iat
      

      if (.not. WF_TMAR_PRINT_READY)
     &     call errquit('WFAT: wfat_initprint_tmdar must first '//
     &     'be called before calling wfat_print_tmdar.', 1, UERR)

      
      write(WF_TMAR_UNIT, '(a)', advance='no') ' '
      write(WF_TMAR_UNIT, '(2es15.5)', advance='no') kx, ky
      
      i = 1
      do ich = 1, nch
         write(WF_TMAR_UNIT, WF_ORBGAP, advance='no') 
         do iat = 1, natom
            write(WF_TMAR_UNIT, WF_ATMGAP, advance='no') 

            write(WF_TMAR_UNIT, WF_TMAR_VALFMT, advance='no')
     &            tmd(iat, ich), abs(tmd(iat, ich))**2

            i = i + natom
         enddo
      enddo

      write(WF_TMAR_UNIT, '(" ")') 
      if (gap) write(WF_TMAR_UNIT, '(" ")') 

      end
c=======================================================      


c=======================================================
c=======================================================
      subroutine wfat_endprint_tmdar()

      implicit none
#include "errquit.fh"
#include "wfat.fh"
#include "stdio.fh"
      

      if (.not. WF_TMAR_PRINT_READY) then
         call errquit('WFAT: wfat_initprint_tmdar must first '//
     &        'be called before calling wfat_endprint_tmdar.',
     &        1, UERR)
      else
         WF_TMAR_PRINT_READY=.false.
      endif

      write(luout, '(/, 3x, "*** INFO ***")')      
      write(luout, '(5x, "Atom-resolved TMD has been printed in ", a)')
     &     ''''//trim(adjustl(WF_TMAR_FNAME))//'''.' 
      
      close(WF_TMAR_UNIT)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_tms(geom, nch, natom, tmds)

      implicit none
#include "stdio.fh"
#include "geom.fh"
#include "errquit.fh"
      integer geom
      integer nch
      integer natom
      double precision tmds(natom, nch)

      character*16 atomname
      character*10 i_str
      integer ia, i, titpos


      titpos = 6 + (9+nch*16)/2
      write(luout, '(2/)')
      call util_print_centered(luout, 'Integrated atom-resolved TMD',
     &     titpos, .true.)
      write(luout, '(" ")')


      call wfat_print_repeat(luout, 9+nch*16, 6, '-')
      write(luout, '(" ")')
      write(luout, '(6x, a4, a5)', advance='no') 'No.', 'Atom'
      call util_flush(luout)
      do i = 1, nch
         write(i_str, '(i10)') i
         write(luout, '(a16)', advance='no')
     &        'Ch. '//trim(adjustl(i_str))
         call util_flush(luout)
      enddo
      write(luout, '(" ")')
      call wfat_print_repeat(luout, 9+nch*16, 6, '-')
      write(luout, '(" ")')

      
      do ia = 1, natom
         if (.not. geom_cent_tag(geom, ia, atomname))
     &        call errquit('WFAT: Getting atomname failed',
     &        52, GEOM_ERR)
         write(luout, '(6x, i4, a5)', advance='no')
     &        ia, trim(adjustl(atomname))
         call util_flush(luout)
         
         do i = 1, nch
            write(luout, '(es16.6)', advance='no') tmds(ia, i)
            call util_flush(luout)
         enddo
         write(luout, '(" ")')
      enddo

      end
c=======================================================
