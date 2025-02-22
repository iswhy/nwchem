c=======================================================
c=======================================================
      subroutine wfat_field_factor(fmagn, orberg_, zion,
     &                             nchan, mchan, fifc)
      
      implicit none
#include "stdio.fh"
      integer nchan, mchan
      double precision zion, fmagn, orberg_
      double precision fifc

      double precision kp

      if (fmagn > 1.0d-13) then
         kp = sqrt(2 * abs(orberg_))
         fifc = kp/2 *
     &        (4*kp**2/fmagn)**(2*zion/kp - 2*nchan - abs(mchan) - 1) *
     &        exp(-2*kp**3/3/fmagn)
      else
         fifc = 0.0d0
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_field_factor(rtdb, nch, nqn, nid, mid, zi, ef,
     &     erg, fifc)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      integer nch
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      double precision zi
      double precision ef
      double precision erg(nch)
      double precision fifc(nqn,nch)

      integer ich, iqn, ffc_type(nch)
      double precision ergf3(nch)

      if (.not. rtdb_get(rtdb, 'wfat:ffc_type', mt_int, nch, ffc_type))
     &     call errquit('WFAT: Cannot read ffc_type from RTDB.',
     &     30, UERR)
      
      do ich = 1, nch
         if (ffc_type(ich) == 3) then
            if (.not. rtdb_get(rtdb, 'wfat:erg_ffc', mt_dbl, nch,
     &           ergf3))
     &           call errquit('WFAT: Cannot read erg_ffc from RTDB.',
     &           30, UERR)
            exit
         endif
      enddo
         
      do ich = 1, nch
         do iqn = 1, nqn
            if (ffc_type(ich) == 1) then
               fifc(iqn, ich) = 1.0d0
            else if (ffc_type(ich) == 2) then
               call wfat_field_factor(ef, erg(ich), zi, nid(iqn),
     &              mid(iqn), fifc(iqn,ich))
            else if (ffc_type(ich) == 3) then
               call wfat_field_factor(ef, ergf3(ich), zi, nid(iqn), 
     &              mid(iqn), fifc(iqn,ich))
            else
               call errquit('WFAT: The input ffc_type of '//
     &              'wfat_get_field_factor corresponds to a '//
     &              'nonexisting option.', 33, UERR)
            endif
         enddo
      enddo

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_get_ergffc(rtdb, nch, erg, ergf)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      integer nch
      double precision erg(nch)
      double precision ergf(nch)

      integer ich, iqn, ffc_type(nch)
      double precision ergf3(nch)
      
      
      if (.not. rtdb_get(rtdb, 'wfat:ffc_type', mt_int, nch, ffc_type))
     &     call errquit('WFAT: Cannot read ffc_type from RTDB.',
     &     44, UERR)

      do ich = 1, nch
         if (ffc_type(ich) == 3) then
            if (.not. rtdb_get(rtdb, 'wfat:erg_ffc', mt_dbl, nch,
     &           ergf3))
     &           call errquit('WFAT: Cannot read erg_ffc from RTDB.',
     &           44, UERR)
            exit
         endif
      enddo

      do ich = 1, nch
         if (ffc_type(ich) == 1) then
            call errquit('WFAT: None of the elements of ffc_type '//
     &           'in wfat_get_ergffc can be equal to 1 if this '//
     &           'subroutine is to be called.', 60, UERR)
         else if (ffc_type(ich) == 2) then
            ergf(ich) = erg(ich)
         else if (ffc_type(ich) == 3) then
            ergf(ich) = ergf3(ich)
         else
            call errquit('WFAT: The input ffc_type of '//
     &           'wfat_get_ergffc corresponds to a nonexisting '//
     &           'option.', 35, UERR)
         endif
      enddo

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_print_field_factor(rtdb, ef, nch, nqn, erg, nid, 
     &     mid, zi, ff)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
      integer rtdb
      double precision ef
      integer nch
      integer nqn
      double precision erg(nch)
      integer nid(nqn)
      integer mid(nqn)
      double precision zi
      double precision ff(nqn,nch)

      character*5 i_str
      character*15 e_str
      character*23 ftps
      logical rtdb_stat
      integer i, j, n, m, ftp(nch)
      double precision ergf(nch)


      rtdb_stat = rtdb_parallel(.false.)
      
      if (.not. rtdb_get(rtdb, 'wfat:ffc_type', mt_int, nch, ftp))
     &     call errquit('WFAT: Cannot read ffc_type from RTDB.',
     &     29, UERR)
      if (.not. rtdb_get(rtdb, 'wfat:erg_ffc', mt_dbl, nch, ergf))
     &     call errquit('WFAT: Cannot read erg_ffc from RTDB.',
     &     29, UERR)

      write(luout, '(2/)')
      call util_print_centered(luout, 'Field factors', 20,
     &     .true.)
      write(luout, '("")')

      write(luout, '(3x, "Field strength = ", f10.6, " a.u.")') ef
      write(luout, '("")')
      
      do i = 1, nch
         write(i_str, '(i5)') i
         if (ftp(i) == 1) then
            ftps = 'unit'
            e_str = 'N/A'
         else if (ftp(i) == 2) then
            ftps = 'orbital energy'
            write(e_str, '(f15.6)') erg(i)
         else if (ftp(i) == 3) then
            ftps = 'modified orbital energy'
            write(e_str, '(f15.6)') ergf(i)
         endif
         
         write(luout, '(3x, "Channel = ", a, ",", 2x, "energy = ", '//
     &        'a, ",", 2x, "type = ", a)')
     &        trim(adjustl(i_str)), trim(adjustl(e_str)),
     &        trim(adjustl(ftps))

         write(luout, '(5x, 2a4, a21)') 'n', 'm', 'field factor'
         call wfat_print_repeat(luout, 29, 5, '-')
         write(luout, '("")')
         do j = 1, nqn
            n = nid(j)
            m = mid(j)
            write(luout, '(5x, 2i4, es21.12)') n, m, ff(j,i)
         enddo
         write(luout, '("")')
      enddo


      rtdb_stat = rtdb_parallel(.true.)
      end
c=======================================================
      
