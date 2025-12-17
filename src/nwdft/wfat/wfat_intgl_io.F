c=======================================================
c=======================================================
      subroutine wfat_store_intgl(rtdb, nsh, ncent, nqn, nch, wfat1,
     &           g_il, g_j, g_k)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "wfat.fh"
#include "global.fh"
#include "errquit.fh"
      integer rtdb     ! [input]
      integer nsh      ! [input]
      integer ncent    ! [input]
      integer nqn      ! [input]
      integer nch      ! [input]
      logical wfat1    ! [input]
      integer g_il     ! [input]
      integer g_j      ! [input]
      integer g_k(2)   ! [input]
      
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


      if (.not. rtdb_cget(rtdb, 'wfat:intgl_sname', 1, fname))
     &     call errquit('WFAT: Cannot read intgl_sname from RTDB.',
     &     1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (isthere) then
         write(luout, '(/, 3x, "*** ATTENTION ***")')
         write(luout, '(5x, "A WFAT integral file named '''//
     &        trim(adjustl(fname))//''' has already existed. This '//
     &        'file will be overwritten.")')
      endif
      
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Writing WFAT integrals to '//
     &     trim(adjustl(fname))//' ... ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionized channels = ", i4)') nch
            

      open(WF_INTGL_UNIT, file=fname, action='write',
     &     form='unformatted')
      call wfat_store_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, nch, g_il)
      if (wfat1) then
         call wfat_store_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, 3*nch, 
     &        g_j)
         call wfat_store_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, nch,
     &        g_k(1))
         call wfat_store_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, nch, 
     &        g_k(2))
      endif
      close(WF_INTGL_UNIT)

      rtdb_stat = rtdb_parallel(.true.)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_load_intgl(rtdb, nsh, ncent, nqn, nch, wfat1,
     &           g_il, g_j, g_k)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "wfat.fh"
#include "global.fh"
#include "errquit.fh"
      integer rtdb     ! [input]
      integer nsh      ! [input]
      integer ncent    ! [input]
      integer nqn      ! [input]
      integer nch      ! [input]
      logical wfat1    ! [input]
      integer g_il     ! [output]
      integer g_j      ! [output]
      integer g_k(2)   ! [output]
      
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
            

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_lname', 17, fname))
     &     call errquit('WFAT: Cannot read intgl_lname from RTDB.',
     &     1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (.not. isthere) then
         call errquit('WFAT: Cannot find the WFAT integral file '''//
     &        trim(adjustl(fname))//'''.', 35, UERR)
      endif
      

      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Reading WFAT integrals from '//
     &     trim(adjustl(fname))//' ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionized channels = ", i4)') nch
      

      open(WF_INTGL_UNIT, file=fname, action='read',
     &     form='unformatted')
      call wfat_load_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, nch, g_il)
      if (wfat1) then
         call wfat_load_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, 3*nch, 
     &        g_j)
         call wfat_load_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, nch,
     &        g_k(1))
         call wfat_load_intgl0(WF_INTGL_UNIT, nsh, ncent, nqn, nch, 
     &        g_k(2))
      endif
      close(WF_INTGL_UNIT)

      rtdb_stat = rtdb_parallel(.true.)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_store_intgl0(u, nsh, ncent, nqn, nch, g_i)

      implicit none
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer u
      integer nsh          ! [input]
      integer ncent        ! [input]
      integer nqn          ! [input]
      integer nch          ! [input]
      integer g_i          ! [input]

      integer l_i, k_i, ni, i1, i2, ich, iqn, icent, i_lm, j
      

      ni = nsh * ncent * nqn
      if (.not. ma_alloc_get(mt_dcpl, ni, 'intgl', l_i, k_i))
     &     call errquit('WFAT: Cannot allocate intgl.', 1, MA_ERR)

      i1 = -ni + 1
      i2 = i1 + ni - 1
      do ich = 1, nch
         i1 = i1 + ni
         i2 = i2 + ni
         call ga_get(g_i, i1, i2, 1, 1, dcpl_mb(k_i), ni)

         j = 0
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  j = j + 1
                  write(u) dcpl_mb(k_i+j-1)
               enddo
            enddo
         enddo

      enddo

      if (.not. ma_free_heap(l_i))
     &     call errquit('WFAT: Cannot deallocate intgl.', 1, MA_ERR)
      
      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_load_intgl0(u, nsh, ncent, nqn, nch, g_i)

      implicit none
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer u
      integer nsh          ! [input]
      integer ncent        ! [input]
      integer nqn          ! [input]
      integer nch          ! [input]
      integer g_i          ! [output]

      integer l_i, k_i, ni, i1, i2, ich, iqn, icent, i_lm, j


      ni = nsh * ncent * nqn
      if (.not. ma_alloc_get(mt_dcpl, ni, 'intgl', l_i, k_i))
     &     call errquit('WFAT: Cannot allocate intgl.', 3, MA_ERR)
      
      i1 = -ni + 1
      i2 = i1 + ni - 1
      do ich = 1, nch
         j = 0
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  j = j + 1
                  read(u) dcpl_mb(k_i+j-1)
               enddo
            enddo
         enddo

         i1 = i1 + ni
         i2 = i2 + ni
         call ga_put(g_i, i1, i2, 1, 1, dcpl_mb(k_i), ni)
      enddo

      if (.not. ma_free_heap(l_i))
     &     call errquit('WFAT: Cannot deallocate intgl.',3, MA_ERR)

      end
c=======================================================
      
