c=================================================
c=================================================
      subroutine mewfat_store_intgl(rtdb, nsh, ncent, nqn, nch, g_il)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "wfat.fh"
#include "global.fh"
#include "errquit.fh"
      integer rtdb                    ! [input]
      integer nsh                     ! [input]
      integer ncent                   ! [input]
      integer nqn                     ! [input]
      integer nch                     ! [input]
      integer g_il                    ! [input]

      character*4096 fname
      logical rtdb_stat, isthere
      integer ndid
      

      ndid = ga_nodeid()
      if (ndid /= 0)
     &     call errquit('MEWFAT: The subroutine ''mewfat_store_'//
     &     'intgl'' cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)


cOLD      if (.not. rtdb_put(rtdb, 'mewfat:gserg0', mt_dbl, nch, gserg0))
cOLD     &     call errquit('MEWFAT: Failed writing gserg0 to RTDB.',
cOLD     &     1, RTDB_ERR)
cOLD      if (.not. rtdb_put(rtdb, 'mewfat:gserg1', mt_dbl, nch, gserg1))
cOLD     &     call errquit('MEWFAT: Failed writing gserg1 to RTDB.',
cOLD     &     1, RTDB_ERR)
cOLD      if (.not. rtdb_put(rtdb, 'mewfat:ionerg', mt_dbl, nch, ionerg))
cOLD     &     call errquit('MEWFAT: Failed writing erg to RTDB.',
cOLD     &     1, RTDB_ERR)
cOLD      if (.not. rtdb_put(rtdb, 'mewfat:edipm0', mt_dbl, 3*nch, edipm0))
cOLD     &     call errquit('MEWFAT: Failed writing edipm0 to RTDB.',
cOLD     &     1, RTDB_ERR)
cOLD      if (.not. rtdb_put(rtdb, 'mewfat:edipm1', mt_dbl, 3*nch, edipm1))
cOLD     &     call errquit('MEWFAT: Failed writing edipm1 to RTDB.',
cOLD     &     1, RTDB_ERR)
cOLD      if (.not. rtdb_put(rtdb, 'mewfat:chcent', mt_dbl, 3*nch, chcent))
cOLD     &     call errquit('MEWFAT: Failed writing chcent to RTDB.',
cOLD     &     1, RTDB_ERR)
cOLD      if (.not. rtdb_put(rtdb, 'mewfat:chmmt', mt_dbl, 3*nch, chmmt))
cOLD     &     call errquit('MEWFAT: Failed writing chmmt to RTDB.',
cOLD     &     1, RTDB_ERR)
cOLD      if (.not. rtdb_put(rtdb, 'mewfat:iondip', mt_dbl, 3*nch, iondip))
cOLD     &     call errquit('MEWFAT: Failed writing iondip to RTDB.',
cOLD     &     1, RTDB_ERR)
      

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_sname', 1, fname))
     &     call errquit('MEWFAT: Cannot read wfat:intgl_sname from '//
     &     'RTDB.', 1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (isthere) then
         write(luout, '(/, 3x, "*** ATTENTION ***")')
         write(luout, '(5x, "A MEWFAT integral file named '''//
     &        trim(adjustl(fname))//''' has already existed. This '//
     &        'file will be overwritten.")')
      endif

      
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Writing MEWFAT integrals to '//
     &     trim(adjustl(fname))//' ... ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionization channels = ", i4)') nch


      open(WF_INTGL_UNIT, file=fname, action='write',
     &     form='unformatted')
      call wfat_store_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, nch, g_il)
      close(WF_INTGL_UNIT)
      
      rtdb_stat = rtdb_parallel(.true.)

      end
c=================================================

      
c=================================================
c=================================================
      subroutine mewfat_load_intgl(rtdb, nsh, ncent, nqn, nch, g_il)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "wfat.fh"
#include "global.fh"
#include "errquit.fh"
      integer rtdb                    ! [input]
      integer nsh                     ! [input]
      integer ncent                   ! [input]
      integer nqn                     ! [input]
      integer nch                     ! [input]
      integer g_il                    ! [output]

      character*4096 fname
      logical rtdb_stat, isthere
      integer ndid

      
      ndid = ga_nodeid()
      if (ndid /= 0)
     &     call errquit('MEWFAT: The subroutine ''mewfat_load_intgl'''//
     &     ' cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)
      

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_lname', 9, fname))
     &     call errquit('MEWFAT: Cannot read wfat:intgl_lname from '//
     &     'RTDB.', 1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (.not. isthere) then
         call errquit('MEWFAT: Cannot find the MEWFAT integral file '//
     &        ''''//trim(adjustl(fname))//'''.', 35, UERR)
      endif


      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Reading WFAT integrals from '//
     &     trim(adjustl(fname))//' ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionization channels = ", i4)') nch
      

      open(WF_INTGL_UNIT, file=fname, action='read', form='unformatted')
      call wfat_load_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, nch, g_il)
      close(WF_INTGL_UNIT)

      rtdb_stat = rtdb_parallel(.true.)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_prep_channelprop(nch, l_eg0, k_eg0, l_eg1, 
     &     k_eg1, l_ed0, k_ed0, l_ed1, k_ed1)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer nch           ! [input]
      integer l_eg0, k_eg0  ! [output]
      integer l_eg1, k_eg1  ! [output]
      integer l_ed0, k_ed0  ! [output]
      integer l_ed1, k_ed1  ! [output]

      if (.not. ma_alloc_get(mt_dbl, nch, 'gserg0', l_eg0, k_eg0))
     &     call errquit('MEWFAT: Cannot allocate gserg0', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nch, 'gserg1', l_eg1, k_eg1))
     &     call errquit('MEWFAT: Cannot allocate gserg1', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nch*3, 'edipm0', l_ed0, k_ed0))
     &     call errquit('MEWFAT: Cannot allocate edipm0', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nch*3, 'edipm1', l_ed1, k_ed1))
     &     call errquit('MEWFAT: Cannot allocate edipm1', 1, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_erase_channelprop(l_eg0, l_eg1, l_ed0, l_ed1)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer l_eg0, l_eg1, l_ed0, l_ed1

      if (.not. ma_free_heap(l_eg0))
     &     call errquit('MEWFAT: Cannot deallocate gserg0', 1, MA_ERR)
      if (.not. ma_free_heap(l_eg1))
     &     call errquit('MEWFAT: Cannot deallocate gserg1', 1, MA_ERR)
      if (.not. ma_free_heap(l_ed0))
     &     call errquit('MEWFAT: Cannot deallocate edipm0', 1, MA_ERR)
      if (.not. ma_free_heap(l_ed1))
     &     call errquit('MEWFAT: Cannot deallocate edipm1', 1, MA_ERR)

      end
c=================================================
      
