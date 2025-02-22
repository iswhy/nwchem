c===============================================================
c===============================================================
      subroutine wfat_print_dens_integxy(rtdb, nvec, vec)

      implicit none
#include "wfat.fh"
#include "wfat_qm_param.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"      
      integer rtdb
      integer nvec
      double precision vec(nbas0, nvec)

      character*255 fname
      integer nx, ny, nz
      integer l_z, k_z, l_rh, k_rh, i_rh, g_rh
      integer iz, i, ndid
      double precision xlims(2), ylims(2), zlims(2), bohr2ang
      double precision x_low, x_upp, y_low, y_upp, z_low, z_upp, dz

      parameter (bohr2ang=0.52917721065d0,
     &     x_low=-50.0d0, x_upp=50.0d0, y_low=-50.0d0, y_upp=50.0d0,
     &     z_low=-25.0d0, z_upp=25.0d0)            ! Default values *_upp and *_low are in Bohr.
      

      call ga_sync()
      ndid = ga_nodeid()
      

c     ==== Z range for printing ====
      if (.not. rtdb_get(rtdb, 'wfat:res_mode:nz_densprint', mt_int,
     &     1, nz))
     &     nz = 500
      if (rtdb_get(rtdb, 'wfat:res_mode:zlims_densprint', mt_dbl,
     &     2, zlims)) then
         zlims = zlims / bohr2ang         ! Convert to Bohr.
      else
         zlims = (/ z_low, z_upp /)       ! Default values already in Bohr.
      endif
      
      if (.not. ma_alloc_get(mt_dbl, nz, 'z_point', l_z, k_z))
     &     call errquit('WFAT: Cannot allocate z_point', 1, MA_ERR)
      dz = (zlims(2) - zlims(1)) / (nz - 1)
      do iz = 1, nz
         dbl_mb(k_z + iz - 1) = zlims(1) + (iz-1)*dz
      enddo


c     ==== X range for integration ====
      if (.not. rtdb_get(rtdb, 'wfat:res_mode:nx_densint', mt_int,
     &     1, nx))
     &     nx = 1000
      if (rtdb_get(rtdb, 'wfat:res_mode:xlims_densint', mt_dbl,
     &     2, xlims)) then
         xlims = xlims / bohr2ang         ! Convert to Bohr.
      else
         xlims = (/ x_low, x_upp /)       ! Default values already in Bohr.
      endif


c     ==== Y range for integration ====      
      if (.not. rtdb_get(rtdb, 'wfat:res_mode:ny_densint', mt_int,
     &     1, ny))
     &     ny = 1000
      if (rtdb_get(rtdb, 'wfat:res_mode:ylims_densint', mt_dbl,
     &     2, ylims)) then
         ylims = ylims / bohr2ang         ! Convert to Bohr.
      else
         ylims = (/ y_low, y_upp /)       ! Default values already in Bohr.
      endif


c     ==== Calculate xy-integrated density ====
      if (.not. ga_create(mt_dbl, nz, 1, 'dens_integxy', 0, 0, g_rh))
     &     call errquit('WFAT: Cannot create dens_integxy', 1, GA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nz*nvec, 'dens_integxy',
     &     l_rh, k_rh))
     &     call errquit('WFAT: Cannot allocate dens_integxy', 1, MA_ERR)

      do i = 1, nvec
         call ga_zero(g_rh)
cforcheck         call wfat_test_numovl(geom_hd, bas_hd, nz, dbl_mb(k_z),
cforcheck     &        nx, xlims, ny, ylims, g_rh, oskel)
         call wfat_vecdens_numinteg_xy(geom_hd, bas_hd, nz, dbl_mb(k_z),
     &        nx, xlims, ny, ylims, vec(1:nbas0,i), g_rh)

         i_rh = k_rh + (i-1)*nz
         call ga_get(g_rh, 1, nz, 1, 1, dbl_mb(i_rh), nz)

         if (ndid==0)
     &        write(luout, '(5x, i4, " All space integral = ", es14.6)')
     &        i, sum(dbl_mb(i_rh:i_rh+nz-1))*dz
      enddo
      call wfat_dtranspose(nz, nvec, dbl_mb(k_rh), dbl_mb(k_rh))


c     ==== Print xy-integrated density ====
      if (ndid == 0) then
         call util_file_name('wrhz', .false., .false., fname)
         open(WF_DENSZ_UNIT, file=fname, action='write')

c        ==== Column titles ====
         write(WF_DENSZ_UNIT, '("#", a14, 3x)', advance='no')
     &        'z (angstrom)'
         call util_flush(WF_DENSZ_UNIT)
         do i = 1, nvec
            write(WF_DENSZ_UNIT, '(i14, 1x)', advance='no') i
            call util_flush(WF_DENSZ_UNIT)
         enddo
         write(WF_DENSZ_UNIT, '(" ")')
         
c        ==== Values ====
         do iz = 1, nz
c            write(WF_DENSZ_UNIT, '(1x, es14.6, 3x)', advance='no')
c     &           dbl_mb(k_z + iz - 1) * bohr2ang
            write(WF_DENSZ_UNIT, '(1x, i6, 3x, es14.6, 3x)',
     &           advance='no')
     &           iz, dbl_mb(k_z + iz - 1) * bohr2ang            
            call util_flush(WF_DENSZ_UNIT)
            
            do i = 1, nvec
               i_rh = k_rh + (iz-1)*nvec + i - 1
               write(WF_DENSZ_UNIT, '(es14.6, 1x)', advance='no')
     &              dbl_mb(i_rh)
               call util_flush(WF_DENSZ_UNIT)
            enddo
            
            write(WF_DENSZ_UNIT, '(" ")')
         enddo
         close(WF_DENSZ_UNIT)
      endif
      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_z))
     &     call errquit('WFAT: Cannot deallocate z_point', 1, MA_ERR)
      if (.not. ma_free_heap(l_rh))
     &     call errquit('WFAT: Cannot deallocate dens_integxy',
     &     1, MA_ERR)
      if (.not. ga_destroy(g_rh))
     &     call errquit('WFAT: Cannot destroy dens_integxy', 1, GA_ERR)

      call ga_sync()

      
      end
c===============================================================            
      
