c===============================================================
c===============================================================
      subroutine wfat_sumyield_lebedev(rtdb, geom, nlas, ftlas, alas, 
     &           wlas, plas, tlas, dlas, nt, tmax, extf, nch, 
     &           erg, chmmt, nqn, nid, mid, lmax, zi, intgl)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "wfat_field_param.fh"
      integer rtdb
      integer geom
      integer nlas
      character(len=FLD_TP_LEN) ftlas(nlas)
      double precision alas(nlas)   ! Peak intensity, TW/cm^2                              
      double precision wlas(nlas)   ! Wavelength, nm                                       
      double precision plas(nlas)   ! Carrier-envelope phase, rad                          
      double precision tlas(nlas)   ! Delay (temporal position of the envelope peak), fs   
      double precision dlas(nlas)   ! Field FWHM, fs                                       
      double precision nt
      double precision tmax
      logical extf
      integer nch                     ! [input], number of ionized channels.
      double precision erg(nch)       ! [input]
      double precision chmmt(3, nch)  ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax                    ! [input]
      double precision zi             ! [input]
      double complex intgl(1:*)       ! [input]
      
      integer NTRIAL, nmax, nsf, nat, i, j
      integer l_xx, l_yy, l_zz, l_w, l_sf, l_sfn, l_egf
      integer k_xx, k_yy, k_zz, k_w, k_sf, k_sfn, k_egf
      parameter (NTRIAL = 10)
      integer nang(NTRIAL)
      double precision bt, gm, deg2rad, pi_const, yield, sumbtgm(NTRIAL)
      parameter (deg2rad = 0.01745329251994d0)
      parameter (pi_const = 3.14159265359d0)


      nang = (/ 6, 14, 26, 38, 50, 74, 86, 146, 194, 266 /)
      nmax = maxval(nang)

      
c     ==== Allocation of Lebedev angular grids and weights ====
      if (.not. MA_alloc_get(MT_DBL, nmax, 'angx', l_xx, k_xx))
     &     call errquit('WFAT: angx allocation failed', 81, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nmax, 'angy', l_yy, k_yy))
     &     call errquit('WFAT: angy allocation failed', 82, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nmax, 'angz', l_zz, k_zz))
     &     call errquit('WFAT: angz allocation failed', 83, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nmax, 'angwt', l_w, k_w))
     &     call errquit('WFAT: angwt allocation failed', 84, MA_ERR)


c     ==== Allocation of structure factors ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 21, GEOM_ERR)
      nsf = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfc', l_sf, k_sf))
     &     call errquit('WFAT: srfc allocation failed', 22, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfcn', l_sfn, k_sfn))
     &     call errquit('WFAT: srfc allocation failed', 23, MA_ERR)

      
c     ==== Modified energy for field factors ====
      if (.not. MA_alloc_get(MT_DBL, nch, 'ergf', l_egf, k_egf))
     &     call errquit('WFAT: ergf allocation failed', 22, MA_ERR)
      call wfat_get_ergffc(rtdb, nch, erg, dbl_mb(k_egf))
      

c     ==== Loop over the sample sets of Lebedev grids ====
      do i = 1, NTRIAL
         call wfat_grid_lebedev(nang(i), dbl_mb(k_xx), dbl_mb(k_yy), 
     &        dbl_mb(k_zz), dbl_mb(k_w))

c        ==== Do the quadrature summation ====
         sumbtgm(i) = 0.0d0
         do j = 1, nang(i)
            bt = acos(dbl_mb(k_zz+j-1)) / deg2rad
            gm = atan(dbl_mb(k_yy+j-1), dbl_mb(k_xx+j-1)) / deg2rad

c           ==== Calculate structure factor ====
            call wfat_get_srfc0(geom, bt, gm, nch, erg, nqn, mid, 
     &           lmax, chmmt, intgl, .false., dcpl_mb(k_sf))
            call wfat_get_srfc0(geom, 180.0d0-bt, 180.0d0+gm, nch, erg,
     &           nqn, mid, lmax, chmmt, intgl, .false., dcpl_mb(k_sfn))

c           ==== Calculate ionization yield ====
            call wfat_calc_yield(nat, zi, nlas, ftlas, alas, wlas, 
     &           plas, tlas, dlas, nt, tmax, extf, nch, erg, nqn,
     &           nid, mid, dcpl_mb(k_sf), dcpl_mb(k_sfn), yield)

c           ==== Integrate using the quadrature formula ====
            sumbtgm(i) = sumbtgm(i) +
     &                   yield * 4 * pi_const * dbl_mb(k_w+j-1)
         enddo

         write(luout, '(4x, i6, es14.6)') nang(i), sumbtgm(i)

         dbl_mb(k_xx:k_xx+nmax-1) = 0.0d0
         dbl_mb(k_yy:k_yy+nmax-1) = 0.0d0
         dbl_mb(k_zz:k_zz+nmax-1) = 0.0d0
         dbl_mb(k_w:k_w+nmax-1) = 0.0d0
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_egf))
     &     call errquit('WFAT: ergf deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_xx))
     &     call errquit('WFAT: angx deallocation failed', 5, MA_ERR)
      if (.not. MA_free_heap(l_yy))
     &     call errquit('WFAT: angy deallocation failed', 5, MA_ERR)
      if (.not. MA_free_heap(l_zz))
     &     call errquit('WFAT: angz deallocation failed', 5, MA_ERR)
      if (.not. MA_free_heap(l_w))
     &     call errquit('WFAT: angwt deallocation failed', 5, MA_ERR)
      if (.not. MA_free_heap(l_sf))
     &     call errquit('WFAT: srfc deallocation failed', 5, MA_ERR)
      if (.not. MA_free_heap(l_sfn))
     &     call errquit('WFAT: srfcn deallocation failed', 5, MA_ERR)

      end
c===============================================================
