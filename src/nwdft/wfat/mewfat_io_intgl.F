cOLDc=================================================
cOLDc=================================================
cOLD      subroutine mewfat_store_intgl_OLD(rtdb, nsh, ncent, nqn, nid, mid,
cOLD     &     lmax_, gqual, smh, rmode, ismask, mask_param, ionerg, chpos,
cOLD     &     chmmt, dipm1, g_il)
cOLD
cOLD      implicit none
cOLD#include "stdio.fh"
cOLD#include "global.fh"
cOLD#include "mafdecls.fh"
cOLD#include "rtdb.fh"
cOLD#include "errquit.fh"
cOLD#include "wfat.fh"
cOLD#include "mewfat_qm_param.fh"
cOLD      integer rtdb                ! [input],
cOLD      integer nsh                 ! [input],
cOLD      integer ncent               ! [input],
cOLD      integer nqn                 ! [input],
cOLD      integer nid(nqn)            ! [input],
cOLD      integer mid(nqn)            ! [input],
cOLD      integer lmax_               ! [input],
cOLD      integer gqual               ! [input],
cOLD      double precision smh        ! [input],
cOLD      logical rmode               ! [input],
cOLD      logical ismask                  ! [input], 
cOLD      double precision mask_param(2)  ! [input]. Only used when ismask=.true.. Must always be supplied regardless.
cOLD      double precision ionerg     ! [output],
cOLD      double precision chpos(3)   ! [output],
cOLD      double precision chmmt(3)   ! [output],
cOLD      double precision dipm1(3)   ! [output],
cOLD      integer g_il                ! [output],
cOLD
cOLD      logical int_nosave
cOLD      integer ndid, cx, cy, nocmo0, nocmo1, npcf, nrcf, nqcf, nscf
cOLD      integer l_mo0, l_mo1, l_pcf, l_rcf, l_qcf, l_scf, l_eg0, l_eg1
cOLD      integer l_ed0, l_ed1
cOLD      integer k_mo0, k_mo1, k_pcf, k_rcf, k_qcf, k_scf, k_eg0, k_eg1
cOLD      integer k_ed0, k_ed1
cOLDcOLD      double precision gserg0, gserg1, edipm0(3), edipm1(3)
cOLD
cOLD      
cOLD      end
cOLDc=================================================

      
c=================================================
c=================================================
      subroutine mewfat_store_intgl(rtdb, nsh, ncent, nqn, nch,
     &     gserg0, gserg1, ionerg, edipm0, edipm1, chcent, chmmt,
     &     iondip, g_il)

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
      double precision gserg0(nch)    ! [input]
      double precision gserg1(nch)    ! [input]
      double precision ionerg(nch)    ! [input]
      double precision edipm0(3,nch)  ! [input]
      double precision edipm1(3,nch)  ! [input]
      double precision chcent(3,nch)  ! [input]
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
     &     call errquit('MEWFAT: The subroutine ''mewfat_store_'//
     &     'intgl'' cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)


      if (.not. rtdb_put(rtdb, 'mewfat:gserg0', mt_dbl, nch, gserg0))
     &     call errquit('MEWFAT: Failed writing gserg0 to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'mewfat:gserg1', mt_dbl, nch, gserg1))
     &     call errquit('MEWFAT: Failed writing gserg1 to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'mewfat:ionerg', mt_dbl, nch, ionerg))
     &     call errquit('MEWFAT: Failed writing erg to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'mewfat:edipm0', mt_dbl, 3*nch, edipm0))
     &     call errquit('MEWFAT: Failed writing edipm0 to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'mewfat:edipm1', mt_dbl, 3*nch, edipm1))
     &     call errquit('MEWFAT: Failed writing edipm1 to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'mewfat:chcent', mt_dbl, 3*nch, chcent))
     &     call errquit('MEWFAT: Failed writing chcent to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'mewfat:chmmt', mt_dbl, 3*nch, chmmt))
     &     call errquit('MEWFAT: Failed writing chmmt to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'mewfat:iondip', mt_dbl, 3*nch, iondip))
     &     call errquit('MEWFAT: Failed writing iondip to RTDB.',
     &     1, RTDB_ERR)
      

      nil_ = nsh * ncent * nqn
      if (.not. ma_alloc_get(mt_dcpl, nil_, 'intgl_', l_il_, k_il_))
     &     call errquit('MEWFAT: Cannot allocate intgl_.', 1, MA_ERR)
      

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('MEWFAT: Cannot read wfat:intgl_name from '//
     &     'RTDB.', 1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (isthere) then
         write(luout, '(/, 3x, "*** ATTENTION ***")')
         write(luout, '(5x, "A MEWFAT integral file named '''//
     &        trim(adjustl(fname))//''' has already existed. This '//
     &        'file will be overwritten.")')
      endif
      open(WF_INTGL_UNIT, file=fname, action='write',
     &     form='unformatted')

      
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Writing MEWFAT integrals to '//
     &     trim(adjustl(fname))//' ... ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionization channels = ", i4)') nch
      

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
     &     call errquit('MEWFAT: Cannot deallocate intgl_.', 1, MA_ERR)
      rtdb_stat = rtdb_parallel(.true.)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_load_intgl(rtdb, nsh, ncent, nqn, nch,
     &     sx, ionerg, chcent, chmmt, iondip, g_il)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb                    ! [input]
      integer nsh                     ! [input]
      integer ncent                   ! [input]
      integer nqn                     ! [input]
      integer nch                     ! [input]
      integer sx                      ! [input]
      double precision ionerg(nch)    ! [output]
      double precision chcent(3,nch)  ! [output]
      double precision chmmt(3,nch)   ! [output]
      double precision iondip(3,nch)  ! [output]
      integer g_il                    ! [output]

      logical rtdb_stat
      integer l_eg0, l_eg1, l_ed0, l_ed1
      integer k_eg0, k_eg1, k_ed0, k_ed1
      integer ndid, state0, state1


      ndid = ga_nodeid()
      if (ndid /= 0)
     &     call errquit('MEWFAT: The subroutine ''mewfat_load_intgl'''//
     &     ' cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)
      

      call mewfat_prep_channelprop(nch, l_eg0, k_eg0, l_eg1, k_eg1,
     &     l_ed0, k_ed0, l_ed1, k_ed1)

      
      call mewfat_load_intgl0(rtdb, nsh, ncent, nqn, nch,
     &     dbl_mb(k_eg0), dbl_mb(k_eg1), ionerg, dbl_mb(k_ed0), 
     &     dbl_mb(k_ed1), chcent, chmmt, iondip, g_il)

      state0 = 0
      state1 = 0
      call mewfat_print_channel_property(nch, state0, state1, 
     &     sx, dbl_mb(k_eg0), dbl_mb(k_eg1), ionerg, dbl_mb(k_ed0), 
     &     dbl_mb(k_ed1), chcent, chmmt, iondip)


      call mewfat_erase_channelprop(l_eg0, l_eg1, l_ed0, l_ed1)

      rtdb_stat = rtdb_parallel(.true.)
      
      end
c=================================================
      

c=================================================
c=================================================
      subroutine mewfat_load_intgl0(rtdb, nsh, ncent, nqn, nch,
     &     gserg0, gserg1, ionerg, edipm0, edipm1, chcent, chmmt,
     &     iondip, g_il)

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
      double precision gserg0(nch)    ! [output]
      double precision gserg1(nch)    ! [output]
      double precision ionerg(nch)    ! [output]
      double precision edipm0(3,nch)  ! [output]
      double precision edipm1(3,nch)  ! [output]
      double precision chcent(3,nch)  ! [output]
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
     &     call errquit('MEWFAT: The subroutine ''mewfat_load_intgl'''//
     &     ' cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)


      if (.not. rtdb_get(rtdb, 'mewfat:gserg0', mt_dbl, nch, gserg0))
     &     call errquit('MEWFAT: Failed reading gserg0 from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'mewfat:gserg1', mt_dbl, nch, gserg1))
     &     call errquit('MEWFAT: Failed reading gserg1 from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'mewfat:ionerg', mt_dbl, nch, ionerg))
     &     call errquit('MEWFAT: Failed reading ionerg from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'mewfat:edipm0', mt_dbl, 3*nch, edipm0))
     &     call errquit('MEWFAT: Failed reading edipm0 from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'mewfat:edipm1', mt_dbl, 3*nch, edipm1))
     &     call errquit('MEWFAT: Failed reading edipm1 from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'mewfat:chcent', mt_dbl, 3*nch, chcent))
     &     call errquit('MEWFAT: Failed reading chcent from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'mewfat:chmmt', mt_dbl, 3*nch, chmmt))
     &     call errquit('MEWFAT: Failed reading chmmt from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'mewfat:iondip', mt_dbl, 3*nch, iondip))
     &     call errquit('MEWFAT: Failed reading iondip from RTDB.',
     &     1, RTDB_ERR)
      

      nil_ = nsh * ncent * nqn
      if (.not. ma_alloc_get(mt_dcpl, nil_, 'intgl_', l_il_, k_il_))
     &     call errquit('MEWFAT: Cannot allocate intgl_.', 1, MA_ERR)
      

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('MEWFAT: Cannot read wfat:intgl_name from '//
     &     'RTDB.', 1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (.not. isthere) then
         call errquit('MEWFAT: Cannot find the MEWFAT integral '//
     &        'file '''//trim(adjustl(fname))//'''. This file is '//
     &        'needed when using ''restart'' as the startup mode. '//
     &        'If it is not specified in the input file, make sure '//
     &        'it is located in the current directory and is named '//
     &        'as ''<file_prefix>.mwigl'' where file_prefix is the '//
     &        'file prefix of NWChem files for the present '//
     &        'calculation. Otherwise, it can be given through the '//
     &        'input file using the ''intgl'' input parameter '//
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
      write(luout, '(7x, "No. of ionization channels = ", i4)') nch
      

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
     &     call errquit('MEWFAT: Cannot deallocate intgl_.', 1, MA_ERR)

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
      
