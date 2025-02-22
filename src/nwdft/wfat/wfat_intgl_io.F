c=======================================================
c=======================================================
      subroutine wfat_store_intgl(rtdb, nsh, ncent, nqn, nch,
     &           erg, chpos, chmmt, iondip, g_il)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "wfat.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb                    ! [input]
      integer nsh                     ! [input]
      integer ncent                   ! [input]
      integer nqn                     ! [input]
      integer nch                     ! [input]
      double precision erg(nch)       ! [input]
      double precision chpos(3,nch)   ! [input]
      double precision chmmt(3,nch)   ! [input]
      double precision iondip(3,nch)  ! [input]
      integer g_il                    ! [input]

      character*4096 fname
      logical rtdb_stat, isthere
      integer ndid, l_il_, k_il_
      integer i1, i2, j, ich, iqn, icent, i_lm
      integer nil_
      

      ndid = ga_nodeid()
      if (ndid /= 0)
     &     call errquit('WFAT: The subroutine ''wfat_store_intgl'' '//
     &     'cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)


      if (.not. rtdb_put(rtdb, 'wfat:erg', mt_dbl, nch, erg))
     &     call errquit('WFAT: Failed writing erg to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:chpos', mt_dbl, 3*nch, chpos))
     &     call errquit('WFAT: Failed writing chpos to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:chmmt', mt_dbl, 3*nch, chmmt))
     &     call errquit('WFAT: Failed writing chmmt to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:iondip', mt_dbl, 3*nch, iondip))
     &     call errquit('WFAT: Failed writing iondip to RTDB.',
     &     1, RTDB_ERR)
      

      nil_ = nsh * ncent * nqn
      if (.not. ma_alloc_get(mt_dcpl, nil_, 'intgl_', l_il_, k_il_))
     &     call errquit('WFAT: Cannot allocate intgl_.', 1, MA_ERR)
      

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('WFAT: Cannot read wfat:intgl_name from RTDB.',
     &     1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (isthere) then
         write(luout, '(/, 3x, "*** ATTENTION ***")')
         write(luout, '(5x, "A WFAT integral file named '''//
     &        trim(adjustl(fname))//''' has already existed. This '//
     &        'file will be overwritten.")')
      endif
      open(WF_INTGL_UNIT, file=fname, action='write',
     &     form='unformatted')

      
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Writing WFAT integrals to '//
     &     trim(adjustl(fname))//' ... ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionized channels = ", i4)') nch
      

      i1 = -nil_ + 1
      i2 = i1 + nil_ - 1
      do ich = 1, nch
         i1 = i1 + nil_
         i2 = i2 + nil_
         call ga_get(g_il, i1, i2, 1, 1, dcpl_mb(k_il_), nil_)

         j = 0
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  j = j + 1
                  write(WF_INTGL_UNIT) dcpl_mb(k_il_+j-1)
               enddo
            enddo
         enddo
      enddo

      close(WF_INTGL_UNIT)

      if (.not. ma_free_heap(l_il_))
     &     call errquit('WFAT: Cannot deallocate intgl_.', 1, MA_ERR)
      rtdb_stat = rtdb_parallel(.true.)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_load_intgl(rtdb, nsh, ncent, nqn, nch,
     &           erg, chpos, chmmt, iondip, g_il)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "wfat.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb                    ! [input]
      integer nsh                     ! [input]
      integer ncent                   ! [input]
      integer nqn                     ! [input]
      integer nch                     ! [input]
      double precision erg(nch)       ! [output]
      double precision chpos(3,nch)   ! [output]
      double precision chmmt(3,nch)   ! [output]
      double precision iondip(3,nch)  ! [output]
      integer g_il                    ! [output]

      character*4096 fname
      logical rtdb_stat, isthere
      integer ndid, l_il_, k_il_
      integer i1, i2, j, ich, iqn, icent, i_lm
      integer nil_

      
      ndid = ga_nodeid()
      if (ndid /= 0)
     &     call errquit('WFAT: The subroutine ''wfat_load_intgl'' '//
     &     'cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)


      if (.not. rtdb_get(rtdb, 'wfat:erg', mt_dbl, nch, erg))
     &     call errquit('WFAT: Failed reading erg from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:chpos', mt_dbl, 3*nch, chpos))
     &     call errquit('WFAT: Failed reading chpos from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:chmmt', mt_dbl, 3*nch, chmmt))
     &     call errquit('WFAT: Failed reading chmmt from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:iondip', mt_dbl, 3*nch, iondip))
     &     call errquit('WFAT: Failed reading iondip from RTDB.',
     &     1, RTDB_ERR)
      

      nil_ = nsh * ncent * nqn
      if (.not. ma_alloc_get(mt_dcpl, nil_, 'intgl_', l_il_, k_il_))
     &     call errquit('WFAT: Cannot allocate intgl_.', 1, MA_ERR)
      

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('WFAT: Cannot read wfat:intgl_name from RTDB.',
     &     1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (.not. isthere) then
         call errquit('WFAT: Cannot find the WFAT integral file '''//
     &        trim(adjustl(fname))//'''. This file is needed when '//
     &        'using ''restart'' as the startup mode. If it is not '//
     &        'specified in the input file, make sure it is located '//
     &        'in the current directory and is named according to '//
     &        '''<file_prefix>.owigl'' where file_prefix is the file '//
     &        'prefix of NWChem files for the present calculation. '//
     &        'Otherwise, it can be given through the input '//
     &        'file using the ''intgl'' input parameter '//
     &        'following the format ''intgl "<string value>"'' '//
     &        '(do not skip the double quotes).',
     &        0, UERR)
      endif
      open(WF_INTGL_UNIT, file=fname, action='read',
     &     form='unformatted')


      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Reading WFAT integrals from '//
     &     trim(adjustl(fname))//' ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionized channels = ", i4)') nch
      

      i1 = -nil_ + 1
      i2 = i1 + nil_ - 1
      do ich = 1, nch
         j = 0
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  j = j + 1
                  read(WF_INTGL_UNIT) dcpl_mb(k_il_+j-1)
               enddo
            enddo
         enddo

         i1 = i1 + nil_
         i2 = i2 + nil_
         call ga_put(g_il, i1, i2, 1, 1, dcpl_mb(k_il_), nil_)
      enddo

      close(WF_INTGL_UNIT)

      if (.not. ma_free_heap(l_il_))
     &     call errquit('WFAT: Cannot deallocate intgl_.', 1, MA_ERR)

      rtdb_stat = rtdb_parallel(.true.)

      end
c=======================================================
      
