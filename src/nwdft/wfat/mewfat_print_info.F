c=======================================================
c=======================================================
      subroutine mewfat_print_info(rtdb)

      implicit none
#include "global.fh"      
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mewfat_qm_param.fh"
#include "geom.fh"
#include "bas.fh"
      integer rtdb

      logical rtdb_stat
      integer ndid
      

      ndid = ga_nodeid()
      if (ndid /= 0) return

      rtdb_stat = rtdb_parallel(.false.)

      
c     ==== Printing info ====
      call wfat_print_printsw(rtdb)      
      

c     ==== Orientation angles ====
      call wfat_print_angles(rtdb)
      
            
c     ==== Remaining WFAT info ====
      call wfat_print_others(rtdb)

      
      rtdb_stat = rtdb_parallel(.true.)
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine mewfat_print_qmsys()

      implicit none
#include "global.fh"
#include "mewfat_qm_param.fh"
#include "stdio.fh"
      integer ndid
      character*5 ss


      ndid = ga_nodeid()

      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Previous SCF info', 
     &        20, .true.)
         write(luout, '("")')
      endif
      
      if (ndid == 0)
     &     write(luout, 900) 'initial species (neutral)',
     &     trim(movecs0), trim(title0), trim(scftype0), trim(basname0),
     &     nbf0, ns_ion0, nmo0, nel0, multip0, charge0

      if (ndid == 0)
     &     write(luout, 900) 'final species (cation)',
     &     trim(movecs1), trim(title1), trim(scftype1), trim(basname1),
     &     nbf1, ns_ion1, nmo1, nel1, multip1, charge1

      if (ndid == 0) then
         if (spinx == 1) then
            ss = 'alpha'
         else if (spinx == 2) then
            ss = 'beta'
         endif
         write(luout, '(3x, "Ionized spin channel : ", a)') ss
      endif
      
 900  format(
     &     3x, 'Information about the ', a, ' :', /,
     &     5x, '* MO file = ', a, /,
     &     5x, '* Job title = ', a, /,
     &     5x, '* SCF type = ', a, /,
     &     5x, '* Basis set name = ', a, /,
     &     5x, '* No. of bases = ', i5, /,
     &     5x, '* No. of spin channels = ', i3, /,
     &     5x, '* No. of MO''s in each channel = ', 2i5, /,
     &     5x, '* No. electrons in each channel = ', 2i4, /,
     &     5x, '* Spin multiplicity = ', i3, /,
     &     5x, '* Charge state = ', f7.3, /)
      
      end
c=======================================================
      
