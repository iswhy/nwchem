c===============================================================
c===============================================================
      subroutine wfat_2a1e_init(rtdb)

      implicit none
#include "stdio.fh"
#include "geom.fh"
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "wfat.fh"
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
      integer rtdb

      logical rtdb_stat
      integer nel_2a1e
      double precision zn


      call ga_sync()
      
c     ==== Printing units ====
      WF_SF_PRINT_READY = .false.
      WF_IY_PRINT_READY = .false.
      WF_SF1_PRINT_READY(1:2) = .false.
      WF_SF1_UNIT = (/ WF_JN_UNIT, WF_KN_UNIT /)
      

c     ==== Check existence, format, and ====
c     ====     necessities of inputs    ====
      if (ga_nodeid() == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         call wfat_process_input(rtdb)
         rtdb_stat = rtdb_parallel(.true.)
      endif

      
      call wfat_2a1e_geom(rtdb, geom_hd)
      if (.not. geom_ncent(geom_hd, natm))
     &     call errquit('WFAT: Getting natm failed', 89, GEOM_ERR)


c     ==== Determine zion ====
      if (.not. geom_nuc_charge(geom_hd, zn))
     &     call errquit('WFAT: Getting zion failed', 52, GEOM_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_nel', mt_int, 1, nel_2a1e))
     &     call errquit('WFAT: Failed reading ipol from RTDB.', 1,
     &     RTDB_ERR)
      zion = zn - nel_2a1e + 1

      
c     ==== Miscellaneous ====
      call wfat_reset_lfqlib()  ! Reset Gauss-Laguerre nodes and weights for the calculations of the various U(a,b,z).

      
c     ==== Printing important information ====
      call wfat_print_info(rtdb)
      
      call ga_sync()
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_geom(rtdb, geom)

      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "geom.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input]
      integer geom              ! [output]

      character*16 tag(2)
      integer nc
      double precision ch(2), cd(3,2), r
      

c     ==== Create the geometry object for 2a1e ====
      if (.not. geom_create(geom, 'wfat_geom'))
     &     call errquit('WFAT: Error in creating geometry object.', 71,
     &     GEOM_ERR)

      
c     ==== Read some quantities from RDB ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_nc', mt_int, 1, nc))
     &     call errquit('WFAT: Failed reading 2a1e_nc from RTDB.', 1,
     &     RTDB_ERR)
      if (nc <= 0 .or. nc >=3)
     &     call errquit('WFAT: The number of 2a1e centers must be '//
     &     'either 1 or 2.', 93, UERR)
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_q', mt_dbl, 2, ch))
     &     call errquit('WFAT: Failed reading 2a1e_q from RTDB.',
     &     1, RTDB_ERR)
      if (nc == 2) then
         if (.not. rtdb_get(rtdb, 'wfat:internuc', mt_dbl, 1, r))
     &        call errquit('WFAT: Failed reading internuc from RTDB.',
     &        1, RTDB_ERR)
         cd = 0.0d0
         cd(3,1) = ch(2)/sum(ch) * r
         cd(3,2) = -ch(1)/sum(ch) * r
      else if (nc == 1) then
         cd = 0.0d0
      endif

      
c     ==== Set the attributes of 2a1e geom using ====
c     ====       the above RTDB quantities       ====
      tag = (/ 'a1', 'a2' /)
      if (.not. geom_cart_set(geom, nc, tag(1:nc), cd(:,1:nc),            ! 1)
     &     ch(1:nc)))
     &     call errquit('WFAT: Error in initializing geometry object.',
     &     71, GEOM_ERR)

      
c     ==== Print geometries ====
      if (ga_nodeid() == 0) then
         if (.not. geom_print(geom))
     &        call errquit('WFAT: Cannot print geometry information '//
     &        'for 2a1e system', 94, UERR)
      endif
      
      end
c     NOTES:
c     1) The coordinates of the two nuclei must be assumed to be defined relative
c        to the basis center (since the 2a1e system employs a single center expansion)
c        throughout the code. This means that in any calculations that involve both
c        the coordinates of the nuclei and electron, e.g. in the code that evaluates 
c        the single center basis at a certain point, the two coordinates should be used
c        with the same reference. If Z is the z-coordinate of one nuclei and z is the
c        z-coordinate for the electron, then if Z-A (with A some reference coordinate)
c        appears, z-A must apply (instead of z-B) wherever the electronic coordinate is
c        needed in any given code.
c===============================================================
      
      
