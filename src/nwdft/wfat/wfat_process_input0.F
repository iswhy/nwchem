c=======================================================
c=======================================================
      subroutine wfat_process_input0(rtdb)

      implicit none
      integer rtdb

      call wfat_process_quadgrid(rtdb)

      call wfat_process_mypos(rtdb)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_quadgrid(rtdb)

      implicit none
#include "wfat_const.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer rtdb

      logical r
      integer ozn, gs
      double precision thr, b2a

      parameter (b2a=0.52917721067d0)


c     ==== gridset ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:gqual', mt_int, 1, gs))
     &     gs = def_gqual
      call wfat_ival_inrange(gs, 1, 4, 'WFAT', '''gridset''', 'i', r)
      if (.not. rtdb_put(rtdb, 'wfat:gqual', mt_int, 1, gs))
     &     call errquit('WFAT: Failed writing gqual to RTDB',
     &     1, RTDB_ERR)

      
c     ==== omegaz ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:ozn', mt_int, 1, ozn))
     &     ozn = DEF_OZN
      call wfat_ival_atleast(ozn, 1, 'WFAT', '''omegaz''', 'i', r)
      if (.not. rtdb_put(rtdb, 'wfat:ozn', mt_int, 1, ozn))
     &     call errquit('WFAT: Failed writing ozn to RTDB',
     &     1, RTDB_ERR)


c     ==== cgridthr ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:cgrid_thr', mt_dbl, 1, thr))
     &     thr = DEF_CGRID_THR
c     DEF_CGRID_THR and wfat:inp:cgrid_thr are in angstrom.
      call wfat_dval_atleast(thr, 0.0d0, 'g', 'WFAT', '''cgridthr''',
     &     'e', 'i', r)
      if (.not. rtdb_put(rtdb, 'wfat:cgrid_thr', mt_dbl, 1, thr/b2a))
     &     call errquit('WFAT: Failed writing cgrid_thr to RTDB',
     &     1, RTDB_ERR)      

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_process_mypos(rtdb)

      implicit none
c#include "wfat_const.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "wfat_input.fh"
      integer rtdb

      integer nch, n
      integer l_s, l_p0, l_p, l_oi, l_si, l_u
      integer k_s, k_p0, k_p, k_oi, k_si, k_u


c     ==== Get lengths ====
      if (wfat_level == 'oe') then
         if (.not. rtdb_get(rtdb, 'wfat:n_orbid', mt_int, 1, nch))
     &        call errquit('WFAT: Failed reading wfat:n_orbid from RTDB'
     &        , 1, INPUT_ERR)
      else if (wfat_level == 'me') then
         nch = 1         ! When using mewfat, nch has to be set manually because orbid_need is not .false. so wfat:n_orbid is not recorded into rtdb.
      endif
      if (.not. rtdb_get(rtdb, 'wfat:inp:chpos_n', mt_int, 1, n))
     &     return
      if (n <= 0) call errquit('WFAT: Bad value for wfat:inp:chpos_n')

      if (n > nch)
     &     call errquit('WFAT: The number of manually specified '//
     &     'channel centers exceed the number of channels.',
     &     78, INPUT_ERR)


      if (.not. MA_alloc_get(MT_BYTE, maxidlen*nch, 'chid', l_s, k_s))
     &     call errquit('WFAT: chid allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, 3*n, 'pos0', l_p0, k_p0))
     &     call errquit('WFAT: pos0 allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, 3*nch, 'pos', l_p, k_p))
     &     call errquit('WFAT: pos allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_INT, nch, 'orbid', l_oi, k_oi))
     &     call errquit('WFAT: orbid allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_INT, nch, 'spinid', l_si, k_si))
     &     call errquit('WFAT: spinid allocation failed', 59, MA_ERR)
      
      
c     ==== Get inputs from RTDB ====
      if (.not. rtdb_cget(rtdb, 'wfat:inp:chpos_id', maxidlen*n,
     &     byte_mb(k_s)))
     &     call errquit('WFAT: Failed reading wfat:inp:chpos_id from '//
     &     'rtdb.', 1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:inp:chpos', mt_dbl, 3*n,
     &     dbl_mb(k_p0)))
     &     call errquit('WFAT: Failed reading wfat:inp:chpos from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:inp:orbid', mt_int, nch,
     &     int_mb(k_oi)))
     &     call errquit('WFAT: Failed reading inp:orbid from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:inp:spinid', mt_int, nch, 
     &     int_mb(k_si))) call errquit('WFAT: Failed writing '//
     &     'inp:spinid to rtdb', 1, INPUT_ERR)


c     ==== Process the inputs ====
      if (.not. MA_alloc_get(MT_LOG, nch, 'umanpos', l_u, k_u))
     &     call errquit('WFAT: umanpos allocation failed', 59, MA_ERR)
      log_mb(k_u:k_u+nch-1) = .false.

      call wfat_process_mypos0(n, byte_mb(k_s), dbl_mb(k_p0), nch,
     &     int_mb(k_oi), int_mb(k_si), log_mb(k_u), dbl_mb(k_p))
      
      
c     ==== Write to RTDB ====
      if (.not. rtdb_put(rtdb, 'wfat:use_manpos', mt_log, nch,
     &     log_mb(k_u)))
     &     call errquit('WFAT: Failed writing wfat:umanpos to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:manpos', mt_dbl, 3*nch,
     &     dbl_mb(k_p)))
     &     call errquit('WFAT: Failed writing wfat:manpos to rtdb',
     &     1, INPUT_ERR)


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_s))
     &     call errquit('WFAT: chid deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_p0))
     &     call errquit('WFAT: pos0 deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_p))
     &     call errquit('WFAT: pos deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_oi))
     &     call errquit('WFAT: orbid deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_si))
     &     call errquit('WFAT: spinid deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_u))
     &     call errquit('WFAT: umanpos deallocation failed', 1, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_mypos0(n, chid, pos0, nch, oid, sid,
     &     umanpos, pos)

      implicit none
#include "wfat_const.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer n
      character(len=maxidlen) chid(n)
      double precision pos0(3,n)
      integer nch
      integer oid(nch)
      integer sid(nch)
      logical umanpos(nch)
      double precision pos(3,nch)

      integer i, ic, info
      

      pos = DEF_MPOS
      
      do i = 1, n
         if (len_trim(chid(i)) > 0) then
            write(luout, '("chid = ", a)') trim(adjustl(chid(i)))

            if (trim(adjustl(chid(i))) == 'all') then
               umanpos(:) = .true.
               pos(1,:) = pos0(1,i)
               pos(2,:) = pos0(2,i)
               pos(3,:) = pos0(3,i)
               exit
            else
               call wfat_read_orb_id(nch, oid, sid, chid(i), ic, info)
               if (info /= 0) then
                  write(luout, '(4x, "info = ", i3)') info
                  write(luout, '(4x, "Bad entry = ", a)')
     &                 trim(adjustl(chid(i)))
                  if (info == 5) then
                     call errquit('WFAT: The above entry of input '//
     &                    'directive mychpos contains an '//
     &                    'unspecified ionized channel.', 8, INPUT_ERR)
                  else
                     call errquit('WFAT: The above entry of input '//
     &                    'mychpos does not conform to the '//
     &                    'correct format.', 8, INPUT_ERR)
                  endif
               endif
               umanpos(ic) = .true.
               pos(:,ic) = pos0(:,i)
            endif
         endif
      enddo
      
      end
c=======================================================
