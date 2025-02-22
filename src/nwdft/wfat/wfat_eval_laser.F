c=======================================================
c=======================================================
      subroutine wfat_eval_laser(nl, ftl, al, wl, pl, tl, dl, t, ef)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "wfat_field_param.fh"
      integer nl
      character(len=FLD_TP_LEN) ftl(nl)         ! Function type of envelope, 'gaussian' or 'hann'
      double precision al(nl)     ! Peak electric field, in a.u. of e. field.
      double precision wl(nl)     ! angular frequency, in a.u. (Hartree).
      double precision pl(nl)     ! Carrier-envelope phase, rad
      double precision tl(nl)     ! Delay (temporal position of the envelope peak), a.u. of time
      double precision dl(nl)     ! Field FWHM, a.u. of time
      double precision t          ! Evaluation time, fs
      double precision ef         ! Electric field, a.u. of field

      integer i
      double precision al0, wl0, tl0, dl0, tt, ef0
      double precision pi_const, lspeed, efau2wcm, fs2au, nm2bohr
      parameter (pi_const=3.14159265359d0, lspeed=137.035999084d0)
      parameter (efau2wcm=3.5094448d16, fs2au=41.3413733352d0)
      parameter (nm2bohr=18.897261255d0)


      ef = 0.0d0
      do i = 1, nl
         al0 = al(i)
         wl0 = wl(i)
         tl0 = tl(i)
         dl0 = dl(i)
         
         if (trim(adjustl(ftl(i))) == 'gaussian') then
            tt = t - tl0
            ef0 = al0 * exp(-4.0d0*log(2.0d0)*tt**2/dl0**2) *
     &            cos(wl0*tt + pl(i))
         else if (trim(adjustl(ftl(i))) == 'hann') then
            call wfat_eval_hann(al0, wl0, pl, tl0, dl0, t, ef0)
         else
            call errquit('WFAT: Unknown value for an element of '//
     &           '''ftl'', which is '//ftl(i)//'.', 1, UERR)
         endif
         ef = ef + ef0
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval_hann(al, wl, pl, tl, dl, t, ef)

      implicit none
c     envelope = A * sin^2(omg_env * (t-t0))
c     carrier = cos(wl*(t-tl) + pl)
      double precision al     ! Peak electric field, in a.u. of e. field.
      double precision wl     ! angular frequency, in a.u. (Hartree).
      double precision pl     ! Carrier-envelope phase, rad
      double precision tl     ! Delay (temporal position of the envelope peak), a.u. of time
      double precision dl     ! [input], FWHM of intensity, a.u. of time.
      double precision t      ! Evaluation time, fs
      double precision ef     ! Electric field, a.u. of field

      double precision tl_, omg_env, t0, lim_l, lim_r, pi_const, envl
      parameter(pi_const=3.14159265359d0)

      
      tl_ = tl
      omg_env = (pi_const - acos(1.0d0 - sqrt(2.0d0))) / dl
      t0 = tl_ - pi_const/2.0d0/omg_env
      
      lim_l = t0
      lim_r = tl_ + pi_const/2.0d0/omg_env
    
      if (t > lim_l .and. t < lim_r) then
         envl = al * (sin(omg_env * (t - t0)))**2
      else
         envl = 0.0d0
      end if

      ef = envl * cos(wl*(t-tl) + pl)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_laser(nl, ftl, al, wl, pl, 
     &                            tl, dl, nt, tmax, extf)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "wfat.fh"
#include "wfat_field_param.fh"
#include "stdio.fh"
      integer nl
      character(len=FLD_TP_LEN) ftl(nl)       ! Function type of envelope, 'gaussian' or 'hann'         
      double precision al(nl)   ! Peak electric field, in a.u.                              
      double precision wl(nl)   ! Angular frequency in a.u. (Hartree)                         
      double precision pl(nl)   ! Carrier-envelope phase, rad                          
      double precision tl(nl)   ! Delay (temporal position of the envelope peak), in a.u. of time   
      double precision dl(nl)   ! Field FWHM, in a.u. of time                                       
      integer nt                    
      double precision tmax     ! Max simulation time, in a.u. of time
      logical extf
      
      character*255 fname
      integer i
      double precision t, t0, dt, las
      parameter (t0=0.0d0)

      double precision wfat_convert_time
      external wfat_convert_time
      

      dt = (tmax-t0) / (nt-1)

      call util_file_name('wlas', .false., .false., fname)
      open(WF_LAS_UNIT, file=fname, action='write')
      do i = 1, nt
         t = t0 + (i-1)*dt
         if (extf) then
            call wfat_eval_extlaser(i, las)
         else
            call wfat_eval_laser(nl, ftl, al, wl, pl, tl, dl, t, las)
         endif

         write(WF_LAS_UNIT, '(i4, 2x, es14.6, 2x es14.6)') i,
     &        wfat_convert_time('au',t,'fs'), las
      enddo
      close(WF_LAS_UNIT)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval_extlaser(i, las)

      implicit none
#include "wfat_input.fh"
      integer i
      double precision las


      las = fld_f_ext(i)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_laser_info(nl, al, las_env, wl, pl, tl, dl,
     &     ntime, tmax)

      implicit none
#include "stdio.fh"
#include "wfat_field_param.fh"
      integer nl
      double precision al(nl)
      character(len=FLD_TP_LEN) las_env(max_nfield)
      double precision wl(nl)
      double precision pl(nl)
      double precision tl(nl)
      double precision dl(nl)
      integer ntime
      double precision tmax

      integer i
      double precision wfat_convert_freq, wfat_convert_time,
     &     wfat_convert_efield, wfat_convert_angle
      external wfat_convert_freq, wfat_convert_time,
     &     wfat_convert_efield, wfat_convert_angle


c     ==== tmax ====
      write(luout, '(3x, a5, " = ", es10.3, " a.u. = ", f7.2, " fs")')
     &     'tmax ', tmax, wfat_convert_time('au', tmax, 'fs')

c     ==== ntime ====
      write(luout, '(3x, a5, " = ", i7)') 'ntime', ntime
      
      do i = 1, nl
         write(luout, '(" ")')
         write(luout, '(3x, "Field no. ", i2)') i
c        ==== type ====
         write(luout, '(5x, "* ", a10, " = ", a)')
     &        'type      ', las_env(i) !byte_mb(k_ftl+(i-1)*FLD_TP_LEN)

         if (las_env(i)=='hann' .or. las_env(i)=='gaussian' .or.
     &        las_env(i)=='dc') then
c           ==== amplitude ====
            write(luout, '(5x, "* ", a10, " = ", f10.3, " a.u.")')
     &           'amplitude ', al(i)
c           ==== intensity ====
            write(luout, '(5x, "* ", a10, " = ", f10.3, " TW/cm^2")')
     &           'intensity ', wfat_convert_efield('au', al(i),
     &           'tw/cm2')
         endif

         if (las_env(i)=='hann' .or. las_env(i)=='gaussian') then
c           ==== frequency ====
            write(luout, '(5x, "* ", a10, " = ", f10.3, " a.u. = ", '//
     &           'f10.3, " THz = ", f10.3, " eV")')
     &           'frequency ',
     &           wl(i),
     &           wfat_convert_freq('au', wl(i), 'thz'),
     &           wfat_convert_freq('au', wl(i), 'ev')

c           ==== wavelength ====
            write(luout, '(5x, "* ", a10, " = ", f10.3, " nm = ", '//
     &           'es10.3, " Bohr")')
     &           'wavelength',
     &           wfat_convert_freq('au', wl(i), 'nm'),
     &           wfat_convert_freq('au', wl(i), 'bohr')
            
c           ==== CEP ====
            write(luout, '(5x, "* ", a10, " = ", f10.3, " rad = ", '//
     &           'f10.3, " degree")')
     &           'CEP       ',
     &           pl(i), wfat_convert_angle('rad', pl(i), 'deg')

c           ==== center ====
            write(luout, '(5x, "* ", a10, " = ", es10.3, " a.u. = ", '//
     &           'f10.3, " fs")')
     &           'center    ',
     &           tl(i), wfat_convert_time('au', tl(i), 'fs')
            
c           ==== width ====
            write(luout, '(5x, "* ", a10, " = ", es10.3, " a.u. = ", '//
     &           'f10.3, " fs")')
     &           'width     ',
     &           dl(i), wfat_convert_time('au', dl(i), 'fs')
         endif
      enddo
         
      end
c=======================================================
      
      
