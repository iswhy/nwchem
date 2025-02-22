c=======================================================
c=======================================================
      subroutine wfat_print_orbtitle(unit0, is, io)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer unit0, is, io

      logical isopen
      character*1 spin0
      character*6 io_str
      character*12 form


      inquire(unit=unit0, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: Argument unit0 of wfat_print_orbtitle '//
     &     'corresponds to an unopened file.', 1, UERR)

cfordebug      write(luout, '("Unit inside wfat_print_orbtitle = ", i7)') unit0
      if ( unit0 /= WF_SF_UNIT .and.
     &     unit0 /= WF_RT_UNIT .and.
     &     unit0 /= WF_IY_UNIT .and.
     &     unit0 /= WF_ABTIY_UNIT .and.
     &     unit0 /= WF_AGMIY_UNIT .and. 
     &     unit0 /= WF_TM_UNIT .and.
     &     unit0 /= WF_TMAR_UNIT .and.
     &     unit0 /= WF_JN_UNIT .and.
     &     unit0 /= WF_KN_UNIT)
     &     call errquit('WFAT: Argument unit0 of wfat_print_orbtitle '//
     &     'corresponds to an unavailable option.', 2, UERR)


      if (io /= 0) then
         write(io_str, '(sp, i5)') io
      else
         write(io_str, '(" ")')
      endif

      if (unit0 == WF_SF_UNIT .or. unit0 == WF_JN_UNIT .or.
     &     unit0 == WF_KN_UNIT) then
         form = WF_SF_TITFMT
      else if (unit0 == WF_RT_UNIT .or. unit0 == WF_IY_UNIT .or. 
     &        unit0 == WF_ABTIY_UNIT .or. unit0 == WF_AGMIY_UNIT) then
         form = WF_IY_TITFMT
      else if (unit0 == WF_TM_UNIT) then
         form = WF_TM_TITFMT
      else if (unit0 == WF_TMAR_UNIT) then
         form = WF_TMAR_TITFMT
      endif
            
      if (is == 1) spin0 = 'a'
      if (is == 2) spin0 = 'b'
      write(unit0, form, advance='no')
     &     spin0//'/HOMO'//trim(adjustl(io_str))
      call flush(unit0)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_qnttitle(unit0, n, m)

      implicit none
#include "stdio.fh"      
#include "wfat.fh"
#include "errquit.fh"
      integer unit0, n, m

      logical isopen
      character*12 form
      character*6 n_str, m_str


      inquire(unit=unit0, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: Argument unit0 of wfat_print_qnttitle '//
     &     'corresponds to an unopened file.', 1, UERR)

      if ( unit0 /= WF_SF_UNIT .and.
     &     unit0 /= WF_RT_UNIT .and.
     &     unit0 /= WF_IY_UNIT .and.
     &     unit0 /= WF_ABTIY_UNIT .and.
     &     unit0 /= WF_AGMIY_UNIT .and.
     &     unit0 /= WF_JN_UNIT .and.
     &     unit0 /= WF_KN_UNIT)
     &     call errquit('WFAT: Argument unit0 of wfat_print_qnttitle '//
     &     'corresponds to an unavailable option.', 2, UERR)


      if (unit0 == WF_SF_UNIT .or. unit0 == WF_JN_UNIT .or.
     &     unit0 == WF_KN_UNIT) then
         form = WF_SF_TITFMT
      else if (unit0 == WF_RT_UNIT .or. unit0 == WF_IY_UNIT .or. 
     &        unit0 == WF_ABTIY_UNIT .or. unit0 == WF_AGMIY_UNIT) then
         form = WF_IY_TITFMT
      endif

      write(n_str, '(i5)') n
      write(m_str, '(sp, i5)') m
      
      write(unit0, form, advance='no')
     &     '('//trim(adjustl(n_str))//','//trim(adjustl(m_str))//')'
      call flush(unit0)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_atmtitle(unit0, iat)

      implicit none
#include "stdio.fh"      
#include "wfat.fh"
#include "errquit.fh"
      integer unit0, iat

      logical isopen
      character*12 form
      character*6 iat_str


      inquire(unit=unit0, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: Argument unit0 of wfat_print_atmtitle '//
     &     'corresponds to an unopened file.', 1, UERR)

      if ( unit0 /= WF_TMAR_UNIT )
     &     call errquit('WFAT: Argument unit0 of wfat_print_atmtitle '//
     &     'corresponds to an unavailable option.', 2, UERR)


      form = WF_TMAR_TITFMT


      write(iat_str, '(i5)') iat      
      write(unit0, form, advance='no') 'Atom #'//trim(adjustl(iat_str))
      call flush(unit0)


      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_partitle(unit0, i)

      implicit none
#include "stdio.fh"      
#include "wfat.fh"
#include "errquit.fh"
      integer unit0, i

      logical isopen
      character*12 form
      character*4 i_str


      inquire(unit=unit0, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: Argument unit0 of wfat_print_partitle '//
     &     'corresponds to an unopened file.', 1, UERR)

      if ( unit0 /= WF_SF_UNIT .and. 
     &     unit0 /= WF_TM_UNIT .and.
     &     unit0 /= WF_TMAR_UNIT .and.
     &     unit0 /= WF_JN_UNIT .and.
     &     unit0 /= WF_KN_UNIT)
     &     call errquit('WFAT: Argument unit0 of wfat_print_partitle '//
     &     'corresponds to an unavailable option.', 2, UERR)


      if (unit0 == WF_SF_UNIT .or. unit0 == WF_JN_UNIT .or.
     &     unit0 == WF_KN_UNIT) then
         form = WF_SF_TITFMT
      else if (unit0 == WF_TM_UNIT) then
         form = WF_TM_TITFMT
      else if (unit0 == WF_TMAR_UNIT) then
         form = WF_TMAR_TITFMT
      endif

      if (i == 1) then
         i_str = 'Real'
      else if (i == 2) then
         i_str = 'Imag'
      else if (i == 3) then
         i_str = 'Abs'
      endif
      
      write(unit0, form, advance='no') i_str
      call flush(unit0)

      end
c=======================================================


      
