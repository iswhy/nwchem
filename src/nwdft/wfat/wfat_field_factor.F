c=======================================================
c=======================================================
      subroutine wfat_field_factor(fmagn, orberg_, zion, nchan, mchan, 
     &     fifc)
      
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
     &     erg, ffc)

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
      double precision ffc(nqn,nch)

      logical ffc_l(nqn,nch)
      integer ic, iq
      double precision mffc(nqn,nch)

      if (.not. rtdb_get(rtdb, 'wfat:use_manffc', mt_log, nqn*nch,
     &     ffc_l)) call errquit('WFAT: Cannot read use_manffc from '//
     &     'RTDB.', 30, UERR)
      if (.not. rtdb_get(rtdb, 'wfat:manffc', mt_dbl, nqn*nch, mffc)) 
     &     call errquit('WFAT: Cannot read manffc from RTDB.', 30, UERR)

      do ic = 1, nch
         do iq = 1, nqn
            if (ffc_l(iq,ic)) then
               ffc(iq,ic) = mffc(iq,ic)
            else
               call wfat_field_factor(ef, erg(ic), zi, nid(iq), mid(iq),
     &              ffc(iq,ic))
            endif
         enddo
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
      character*4 ftps
      logical ffc_l(nqn,nch)
      integer i, j

      
      if (.not. rtdb_get(rtdb, 'wfat:use_manffc', mt_log, nqn*nch,
     &     ffc_l)) call errquit('WFAT: Cannot read use_manffc from '//
     &     'RTDB.', 29, UERR)

      write(luout, '(2/)')
      call util_print_centered(luout, 'Field factors', 20, .true.)
      write(luout, '("")')

      write(luout, '(3x, "Field strength = ", f10.6, " a.u.")') ef
      write(luout, '("")')
      
      do i = 1, nch
         write(i_str, '(i5)') i
         write(e_str, '(f15.6)') erg(i)         
         write(luout, '(3x, "Channel = ", a, ",", 2x, "energy = ", '//
     &        'a)') trim(adjustl(i_str)), trim(adjustl(e_str))

         write(luout, '(5x, 2a4, a21, 2x, a4)')
     &        'n', 'm', 'field factor', 'type'
         call wfat_print_repeat(luout, 35, 5, '-')
         write(luout, '("")')
         
         do j = 1, nqn
            if (ffc_l(j,i)) then
               ftps = 'user'
            else
               ftps = 'auto'
            endif
            write(luout, '(5x, 2i4, es21.12, 2x, a4)')
     &           nid(j), mid(j), ff(j,i), ftps
         enddo
         write(luout, '("")')
      enddo

      end
c=======================================================
      
