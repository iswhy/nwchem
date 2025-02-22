c========================================================
c========================================================
      subroutine geom_cent2cent(geom, icent1, v1, icent2, v2, tt)

      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [Input] the geometry handle
      integer icent1            !< [Input] the center1 rank
      integer icent2            !< [Input] the center2 rank
      double precision v1(3)    !< [Input] the coordinate w.r.t icent1
      double precision v2(3)    !< [Output] the coordinate w.r.t icent2
      logical tt                !< [Output] true if successful, false otherwise
      double precision at1(3), at2(3)
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c

      if (icent1 /= 0 .or. icent2 /=0) then
         tt = geom_check_handle(geom, 'geom_cent2cent')
         if (.not. tt) return
      endif

      if (icent1 /= 0) then
         tt = geom_check_cent(geom, 'geom_cent2cent - icent1', icent1)
         if (.not. tt) return
         at1(1:3) = coords(1:3,icent1,geom)
      else
         at1(1:3) = 0.0d0
      endif

      if (icent2 /= 0) then
         tt = geom_check_cent(geom, 'geom_cent2cent - icent2', icent2)
         if (.not. tt) return
         at2(1:3) = coords(1:3,icent2,geom)
      else
         at2(1:3) = 0.0d0
      endif
      
      v2 = v1 + at1 - at2
      tt = .true.
c
      end
c========================================================


c========================================================
c========================================================
      subroutine geom_cent2arb(geom, icent1, v1, cent2, v2, tt)

      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
      integer icent1
      double precision v1(3)
      double precision cent2(3)
      double precision v2(3)
      logical tt

      double precision at1(3)
      
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      
      if (icent1 /= 0) then
         tt = geom_check_handle(geom, 'geom_cent2cent')
         if (.not. tt) return
         tt = geom_check_cent(geom, 'geom_cent2arb - icent1', icent1)
         if (.not. tt) return
         at1(1:3) = coords(1:3, icent1, geom)
      else
         at1(1:3) = 0.0d0
      endif
      
      v2 = v1 + at1 - cent2

      end
c========================================================


c========================================================
c========================================================
      subroutine geom_arb2cent(geom, cent1, v1, icent2, v2, tt)

      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
      double precision cent1(3)
      double precision v1(3)
      integer icent2
      double precision v2(3)
      logical tt

      double precision at2(3)
      
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      
      if (icent2 /= 0) then
         tt = geom_check_handle(geom, 'geom_cent2cent')
         if (.not. tt) return
         tt = geom_check_cent(geom, 'geom_cent2arb - icent2', icent2)
         if (.not. tt) return
         at2(1:3) = coords(1:3, icent2, geom)
      else
         at2(1:3) = 0.0d0
      endif
      
      v2 = v1 + cent1 - at2

      end
c========================================================


c========================================================
c========================================================
      subroutine wfat_channel_geom0(geom, natom, tag, cd, cg, geomo)

      implicit none
#include "errquit.fh"
#include "geom.fh"
      integer geom
      integer natom
      character*16 tag
      double precision cd(3)
      double precision cg
      integer geomo

      character*16 tag_(natom+1)
      integer ic
      double precision cd_(3,natom+1)
      double precision cg_(natom+1)


      if (.not. geom_check_handle(geomo, 'geomo does not exist'))
     &     call errquit('WFAT: The input geomo must be an existing'//
     &     'geometry handle before input to wfat_channel_geom0.', 33,
     &     GEOM_ERR)
      if (.not. geom_cart_set(geomo, natom+1, tag_, cd_, cg_))
     &     call errquit('WFAT: Error in initializing geometry object.',
     &     65, GEOM_ERR)
      

      if (.not. geom_cent_set(geomo, 1, tag, cd, cg))
     &     call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)

      
      do ic = 2, natom+1
         if (.not. geom_cent_get(geom, ic-1, tag_(1), cd_(:,1), cg_(1)))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
         if (.not. geom_cent_set(geomo, ic, tag_(1), cd_(:,1), cg_(1)))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
      enddo

      end
c========================================================


c========================================================
c========================================================
      subroutine wfat_cgeom_create(geom, ozn, nch, chpos, geomc)

      implicit none
#include "geom.fh"
#include "errquit.fh"
      integer geom
      integer ozn
      integer nch
      double precision chpos(3,nch)
      integer geomc(nch)

      character*16 ctag
      integer i, natom

      if (.not. geom_ncent(geom, natom))
     &     call errquit('WFAT: Cannot get natom.', 30, GEOM_ERR)
      
      do i = 1, nch
         if (.not. geom_create(geomc(i), 'geomc'))
     &        call errquit('WFAT: Error in creating geomc.', 30,
     &        GEOM_ERR)
         ctag = 'chcent'
         call wfat_channel_geom0(geom, natom, ctag, chpos(:,i), 
     &        dble(ozn), geomc(i))
      enddo

      end
c========================================================


c========================================================
c========================================================
      subroutine wfat_cgeom_destroy(nch, geomc)

      implicit none
#include "geom.fh"
#include "errquit.fh"
      integer nch
      integer geomc(nch)

      integer i
      
      do i = 1, nch
         if (.not. geom_destroy(geomc(i)))
     &        call errquit('WFAT: Cannot destroy geomc.', 30,
     &        GEOM_ERR)
      enddo
      
      end
c========================================================
      
      
      
