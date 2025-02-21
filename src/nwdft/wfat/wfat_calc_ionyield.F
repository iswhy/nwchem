c=======================================================
c=======================================================
      subroutine wfat_calc_rate(nat, nch, nqn, srfc, ffc, p_ionyd, 
     &     ionyd)

      implicit none
      integer nat                          ! [input]
      integer nch                          ! [input]
      integer nqn                          ! [input]
      double complex srfc(1:*)             ! [input]
      double precision ffc(nqn, nch)       ! [input]
      double precision p_ionyd(nqn, nch)   ! [output]
      double precision ionyd               ! [output]
      
      integer ich, iqn, i
      double complex srfc_tot      

      i = 1
      do ich = 1, nch
         do iqn = 1, nqn
            srfc_tot = sum(srfc(i : i+nat-1))
            p_ionyd(iqn, ich) = ffc(iqn, ich) * abs(srfc_tot)**2

            i = i + nat
         enddo
      enddo
      ionyd = sum(p_ionyd)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_yield(nat, zi, nlas, ftlas, alas, wlas, 
     &                 plas, tlas, dlas, nt, tmax, extf, nch, ergf, nqn,
     &                 nid, mid, srfc, srfcn, yield)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "wfat_field_param.fh"
      integer nat
      double precision zi
      integer nlas
      character(len=FLD_TP_LEN) ftlas(nlas)       ! Function type of envelope, 'gaussian' or 'hann'         
      double precision alas(nlas)   ! Peak electric field, in a.u. of e. field.                              
      double precision wlas(nlas)   ! angular frequency, in a.u. (Hartree).                                      
      double precision plas(nlas)   ! Carrier-envelope phase, rad                          
      double precision tlas(nlas)   ! Delay (temporal position of the envelope peak), a.u. of time
      double precision dlas(nlas)   ! Field FWHM, a.u. of time                                       
      integer nt                    
      double precision tmax         ! Time window width, a.u. of time
      logical extf
      integer nch
      double precision ergf(nch)
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      double complex srfc(1:*)
      double complex srfcn(1:*)
      double precision yield

      integer i, j, k, npi, l_pi, k_pi, nfc, l_fc, k_fc, i_fc, ich, iqn
      double precision t0, dt, t, las, rate
      parameter (t0=0.0d0)

      
      npi = nqn * nch
      if (.not. ma_alloc_get(mt_dbl, npi, 'partial_ion', l_pi, k_pi))
     &     call errquit('WFAT: Cannot allocate partial_ion', 1, MA_ERR)
      nfc = nqn * nch
      if (.not. MA_alloc_get(MT_DBL, nfc, 'ffc', l_fc, k_fc))
     &     call errquit('WFAT: ffc allocation failed', 2, MA_ERR)
      
      
      dt = (tmax-t0)/(nt-1)

      yield = 0.0d0
      do i = 1, nt
         t = t0 + (i-1)*dt

c        ==== Evaluate laser field ====
         if (extf) then
            call wfat_eval_extlaser(i, las)
         else
            call wfat_eval_laser(nlas, ftlas, alas, wlas, plas, tlas, 
     &           dlas, t, las)
         endif

c        ==== Calculate field factors ====
         do ich = 1, nch
            do iqn = 1, nqn
               i_fc = k_fc + (ich-1)*nqn + (iqn-1)
               if (las >= 1.0d-13) then
                  call wfat_field_factor(las, ergf(ich), zi, nid(iqn),
     &                 mid(iqn), dbl_mb(i_fc))
               else if (las <= -1.0d-13) then
                  call wfat_field_factor(abs(las), ergf(ich), zi, 
     &                 nid(iqn), mid(iqn), dbl_mb(i_fc))
               else
                  dbl_mb(i_fc) = 0.0d0
               endif
            enddo
         enddo

c        ==== Calculate rates ====
         dbl_mb(k_pi:k_pi+npi-1) = 0.0d0
         if (las >= 1.0d-13) then
            call wfat_calc_rate(nat, nch, nqn, srfc, dbl_mb(k_fc),
     &           dbl_mb(k_pi), rate)
         else if (las <= -1.0d-13) then
            call wfat_calc_rate(nat, nch, nqn, srfcn, dbl_mb(k_fc),
     &           dbl_mb(k_pi), rate)
         else
            dbl_mb(k_pi : k_pi+npi-1) = 0.0d0
            rate = 0.0d0
         endif

c        ==== Calculate yield ====
         yield = yield + rate * dt
      enddo

      yield = 1.0d0 - exp(-yield)

      if (.not. ma_free_heap(l_pi))
     &     call errquit('WFAT: Cannot deallocate partial_ion',
     &     1, MA_ERR)
      if (.not. MA_free_heap(l_fc))
     &     call errquit('WFAT: ffc deallocation failed',
     &     3, MA_ERR)
      
      end
c=======================================================
