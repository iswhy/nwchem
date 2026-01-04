c=======================================================
c=======================================================
      subroutine wfat_process_input0(rtdb)

      implicit none
#include "wfat_const.fh"
#include "errquit.fh"
#include "inp.fh"
      integer rtdb

      character*1 w
      logical r
      integer ival, k
      double precision dval
      character*1000 aval
      character*5 si_opt(3)
      parameter (w='i')


c     ==== Print ====
      call wfat_process_print(rtdb)
      
c     ==== Field ====
      call wfat_process_field(rtdb)

c     ==== quadgrid ====
      call wfat_process_quadgrid(rtdb)

c     ==== mychpos ====
      call wfat_process_mypos(rtdb)

c     ==== angles ====
      call wfat_process_angle(rtdb, 'beta')
      call wfat_process_angle(rtdb, 'gamma')
      
c     ==== kx and ky ====
      call wfat_process_momentum(rtdb)

c     ==== intgl ====
      call wfat_process_intgl(rtdb)

c     ==== lmax ====
      call wfat_process_idata(rtdb, 'lmax', .false., def_wfat_lmax,
     &     ival)
            
c     ==== nlgl ====
      call wfat_process_idata(rtdb, 'nlgl', .false., def_nlgl, ival)
      call wfat_ival_atleast(ival, 0, 'WFAT', '''nlgl''', w, r)
      
c     ==== degenthr ====
      call wfat_process_ddata(rtdb, 'dgnthr', .false., def_dgnthr,
     &     dval)
      call wfat_dval_atleast(dval, 0.0d0, 'g', 'WFAT', '''degenthr''',
     &     'e', w, r)

c     ==== diagthr ====
      call wfat_process_ddata(rtdb, 'dthr', .false., def_wfat_dthr,
     &     dval)
      call wfat_dval_atleast(dval, 0.0d0, 'g', 'WFAT', '''diagthr''', 
     &     'e', w, r)
      
c     ==== si_type ====
      si_opt = (/ def_si_type//' ', 'mo   ', 'dyson' /)
      call wfat_process_adata(rtdb, 'si_type', .false., def_si_type,
     &     aval)
      if (.not. inp_match(3, .false., aval, si_opt, k))
     &     call errquit('WFAT: Choose among none, mo, or dyson for '//
     &     'the value of the si_type directive.', 25, INPUT_ERR)
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_adata(rtdb, key, md, defval, aval)
c     DESCRIPTION:
c     Checks for an integer-valued entry of the form wfat:inp:<key> in the RTDB. If
c     found, it will create a new RTDB entry named wfat:<key> with the same value
c     as wfat:inp:<key>.
c       rtdb: RTDB handle
c       key: name of input
c       md: if .true., the input is mandatory, otherwise, optional
c       defval: default value that will be used if input is missing and md=.false.
c               Regardless of md, a value for md must be supplied.
c       aval: the value of the input
      implicit none
#include "errquit.fh"
#include "rtdb.fh"      
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb             ![input]
      character*(*) key        ![input]
      logical md               ![input]
      character*(*) defval     ![input]
      character*(*) aval       ![output]

      if (.not. rtdb_cget(rtdb, 'wfat:inp:'//key, 1, aval)) then
         if (md) then
            call errquit('Mandatory WFAT input directive '//key//' '//
     &           'cannot be found. Either supply it in the input '//
     &           'file or, if it is a restart job, make sure that it '//
     &           'was defined in the previous job the present one is '//
     &           'restarting from.', 30, INPUT_ERR)
         else
            write(luout, 100) key, defval
            aval = defval
         endif
      endif
      if (.not. rtdb_cput(rtdb, 'wfat:'//key, 1, trim(adjustl(aval))) )
     &     call errquit('WFAT: Failed writing '//key//' to rtdb',
     &     1, RTDB_ERR)

 100  format(
     &     /, 3x, '*** INFO ***',
     &     /, 5x, 'Input directive ', a, ' is not found. The default ',
     &     'value of ', a, ' will be assumed.')

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_process_ldata(rtdb, key, md, defval, lval)
c     DESCRIPTION:
c     Checks for a logical-valued entry of the form wfat:inp:<key> in the RTDB. If
c     found, it will create a new RTDB entry named wfat:<key> with the same value
c     as wfat:inp:<key>.
c       rtdb: RTDB handle
c       key: name of input
c       md: if .true., the input is mandatory, otherwise, optional
c       defval: default value that will be used if input is missing and md=.false.
c               Regardless of md, a value for md must be supplied.
c       lval: the value of the input
      implicit none
#include "errquit.fh"
#include "rtdb.fh"      
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb       ![input]
      character*(*) key  ![input]
      logical md         ![input]
      logical defval     ![input]
      logical lval       ![output]

      character*5 defs
      
      if (.not. rtdb_get(rtdb, 'wfat:inp:'//key, mt_log, 1, lval)) then
         if (md) then
            call errquit('Mandatory WFAT input directive '//key//' '//
     &           'cannot be found. Either supply it in the input '//
     &           'file or, if it is a restart job, make sure that it '//
     &           'was defined in the previous job the present one is '//
     &           'restarting from.', 30, INPUT_ERR)
         else
            if (defval) then
               defs = 'True '
            else
               defs = 'False'
            endif
            write(luout, 100) key, trim(defs)
            lval = defval
         endif
      endif
      if (.not. rtdb_put(rtdb, 'wfat:'//key, mt_log, 1, lval))
     &     call errquit('WFAT: Failed writing '//key//' to rtdb',
     &     1, RTDB_ERR)

 100  format(
     &     /, 3x, '*** INFO ***',
     &     /, 5x, 'Input directive ', a, ' is not found. The default ',
     &     'value of ', a, ' will be assumed.')

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_process_idata(rtdb, key, md, defval, ival)
c     DESCRIPTION:
c     Checks for an integer-valued entry of the form wfat:inp:<key> in the RTDB. If
c     found, it will create a new RTDB entry named wfat:<key> with the same value
c     as wfat:inp:<key>.
c       rtdb: RTDB handle
c       key: name of input
c       md: if .true., the input is mandatory, otherwise, optional
c       defval: default value that will be used if input is missing and md=.false.
c               Regardless of md, a value for md must be supplied.
c       ival: the value of the input
      implicit none
#include "errquit.fh"
#include "rtdb.fh"      
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb       ![input]
      character*(*) key  ![input]
      logical md         ![input]
      integer defval     ![input]
      integer ival       ![output]

      if (.not. rtdb_get(rtdb, 'wfat:inp:'//key, mt_int, 1, ival)) then
         if (md) then
            call errquit('Mandatory WFAT input directive '//key//' '//
     &           'cannot be found. Either supply it in the input '//
     &           'file or, if it is a restart job, make sure that it '//
     &           'was defined in the previous job the present one is '//
     &           'restarting from.', 30, INPUT_ERR)
         else
            write(luout, 100) key, defval
            ival = defval
         endif
      endif
      if (.not. rtdb_put(rtdb, 'wfat:'//key, mt_int, 1, ival))
     &     call errquit('WFAT: Failed writing '//key//' to rtdb',
     &     1, RTDB_ERR)

 100  format(
     &     /, 3x, '*** INFO ***',
     &     /, 5x, 'Input directive ', a, ' is not found. The default ',
     &     'value of ', i4, ' will be assumed.')

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_ddata(rtdb, key, md, defval, dval)
c     DESCRIPTION:
c     Checks for a decimal-valued entry of the form wfat:inp:<key> in the RTDB. If
c     found, it will create a new RTDB entry named wfat:<key> with the same value
c     as wfat:inp:<key>.
c       rtdb: RTDB handle
c       key: name of input
c       md: if .true., the input is mandatory, otherwise, optional
c       defval: default value that will be used if input is missing and md=.false.
c               Regardless of md, a value for defval must be supplied.
c       dval: the value of the input
      implicit none
#include "errquit.fh"
#include "rtdb.fh"      
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb       ![input]
      character*(*) key  ![input]
      logical md         ![input]
      double precision defval     ![input]
      double precision dval       ![output]

      if (.not. rtdb_get(rtdb, 'wfat:inp:'//key, mt_dbl, 1, dval)) then
         if (md) then
            call errquit('Mandatory WFAT input directive '//key//' '//
     &           'cannot be found. Either supply it in the input '//
     &           'file or, if it is a restart job, make sure that it '//
     &           'was defined in the previous job the present one is '//
     &           'restarting from.', 30, INPUT_ERR)
         else
            write(luout, 100) key, defval
            dval = defval
         endif
      endif
      if (.not. rtdb_put(rtdb, 'wfat:'//key, mt_dbl, 1, dval))
     &     call errquit('WFAT: Failed writing '//key//' to rtdb',
     &     1, RTDB_ERR)

 100  format(
     &     /, 3x, '*** INFO ***',
     &     /, 5x, 'Input directive ', a, ' is not found. The default ',
     &     'value of ', es13.6, ' will be assumed.')

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_print(rtdb)

      implicit none
#include "errquit.fh"
      integer rtdb

      logical sr, rt, yi, avb, avg, tmd, tmda, md, def, test1, test2
      parameter (md=.false., def=.false.)

      call wfat_process_ldata(rtdb, 'p_srfc', md, def, sr)
      call wfat_process_ldata(rtdb, 'p_rate', md, def, rt)
      call wfat_process_ldata(rtdb, 'p_yild', md, def, yi)
      call wfat_process_ldata(rtdb, 'p_avbeta', md, def, avb)
      call wfat_process_ldata(rtdb, 'p_avgamma', md, def, avg)
      call wfat_process_ldata(rtdb, 'p_tmd', md, def, tmd)
      call wfat_process_ldata(rtdb, 'p_tmd_atom', md, def, tmda)

      test1 = (rt .or. tmd .or. tmda)
      test2 = (yi .or. avb .or. avg)
      if (test1 .and. test2)
     &     call errquit('WFAT: If any of the keywords rate, tmd, or '//
     &     'tmd_atom are specified for the print directive, then none'//
     &     ' of the keywords yield, avbeta, or avgamma can be '//
     &     'specified simultaneously, and vice versa.', 63, INPUT_ERR)

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
      if (n <= 0) call errquit('WFAT: Bad value for wfat:inp:chpos_n',
     &     22, INPUT_ERR)

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


c=======================================================
c=======================================================
      subroutine wfat_process_angle(rtdb, c)

      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh" 
      integer rtdb
      character*(*) c

      character*1 w
      logical r
      integer i, na
      double precision a0, an, ulim
      parameter( w='i' )


      if (c == 'beta') then
         ulim = 180.0d0
      else if (c== 'gamma') then
         ulim = 360.0d0
      else
         call errquit('WFAT: Bad value for the argument c of '//
     &        'wfat_process_angle', 50, UERR)
      endif      
         

c     ==== Start angle ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:'//c//'_start', mt_dbl, 1, a0))
     &     call errquit('WFAT: Cannot read inp:'//c//'_start from rtdb',
     &     35, INPUT_ERR)
      call wfat_dval_inrange(a0, 0.0d0, 'q', ulim, 'q', 'WFAT', 
     &     'first '//c, 'f', w, r)
      if (.not. rtdb_put(rtdb, 'wfat:'//c//'_start', mt_dbl, 1, a0))
     &     call errquit('WFAT: Cannot write '//c//'_start to rtdb', 
     &     35, INPUT_ERR)

c     ==== No. of angles ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:n'//c, mt_int, 1, na))
     &     call errquit('WFAT: Cannot read n'//c//' from rtdb', 
     &     35, INPUT_ERR)
      call wfat_ival_atleast(na, 1, 'WFAT', '''n'//c//'''', 'i', r)
      if (.not. rtdb_put(rtdb, 'wfat:n'//c, mt_int, 1, na))
     &     call errquit('WFAT: Cannot write n'//c//' to rtdb', 
     &     35, INPUT_ERR)

c     ==== Last angle ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:'//c//'_last', mt_dbl, 1, an))
     &     call errquit('WFAT: Cannot read inp:'//c//'_last from rtdb',
     &     35, INPUT_ERR)
      if (na == 1) then
         call wfat_dval_inrange(an, a0, 'q', ulim, 'q', 'WFAT', 
     &        'last '//c, 'f', w, r)
      else  ! na > 1
         call wfat_dval_inrange(an, a0+1.0d-12, 'q', ulim, 'q', 
     &        'WFAT', 'last '//c, 'f', 'o', r)
         if (.not. r) call errquit('WFAT: If there are more than one '//
     &        'angles, the start and last angles should not be too '//
     &        'close.', 92, INPUT_ERR)
      endif
      if (.not. rtdb_put(rtdb, 'wfat:'//c//'_last', mt_dbl, 1, an))
     &     call errquit('WFAT: Cannot write '//c//'_last to rtdb', 
     &     35, INPUT_ERR)
      
      

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_process_momentum(rtdb)

      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh" 
      integer rtdb

      character*1 cc(2)
      character*1 c
      logical r, p_tmd
      integer i, nk


      cc = (/ 'x', 'y' /)

      if (.not. rtdb_get(rtdb, 'wfat:inp:p_tmd', mt_log, 1, p_tmd))
     &     p_tmd = .false.

      if (p_tmd) then
         do i = 1, 2
            c = cc(i)
            
            call wfat_rtdb_fduplicate(rtdb, 'wfat:inp:k'//c//'_start',
     &           'wfat:k'//c//'_start')
            call wfat_rtdb_fduplicate(rtdb, 'wfat:inp:k'//c//'_last',
     &           'wfat:k'//c//'_last')
            if (.not. rtdb_get(rtdb, 'wfat:inp:nk'//c, mt_int, 1, nk))
     &           call errquit('WFAT: Cannot read nk'//c//' from rtdb', 
     &           35, INPUT_ERR)
            call wfat_ival_atleast(nk, 1, 'WFAT', '''nk'//c//'''', 'i',
     &           r)
            if (.not. rtdb_put(rtdb, 'wfat:nk'//c, mt_int, 1, nk))
     &           call errquit('WFAT: Cannot write nk'//c//' to rtdb', 
     &           35, INPUT_ERR)
         enddo
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_intgl(rtdb)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      
      character*4096 fname
      logical loadint

c     ==== Save file ====
      if (rtdb_cget(rtdb, 'wfat:inp:intgl_sname', 1, fname)) then
         call wfat_rtdb_aduplicate(rtdb, 'wfat:inp:intgl_sname',
     &        'wfat:intgl_sname')
      else
         call wfat_def_intglsname(rtdb, fname)
         if (.not. rtdb_cput(rtdb, 'wfat:intgl_sname', 1,
     &        trim(adjustl(fname))))
     &        call errquit('WFAT: Failed writing intgl_sname to RTDB.',
     &        31, RTDB_ERR)
      endif
c     Note that in the above conditional, the default filename for the saved
c     integral file will be used if wfat:inp:intgl_sname cannot be found in
c     the RTDB. However, in a restart job, even though the saved integral
c     file is not specified in the input, wfat:inp:intgl_sname may still be
c     found in the RTDB if the previous job specifies a saved file name.


c     ==== Load file ====
      if (rtdb_cget(rtdb, 'wfat:inp:intgl_lname', 1, fname)) then
         call wfat_rtdb_aduplicate(rtdb, 'wfat:inp:intgl_lname',
     &        'wfat:intgl_lname')
         loadint = .true.
      else
         loadint = .false.
      endif

      if (.not. rtdb_put(rtdb, 'wfat:load_intgl', mt_log, 1, loadint)) 
     &     call errquit('WFAT: Cannot write load_intgl from RTDB',
     &     1, RTDB_ERR)

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_process_field(rtdb)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "wfat_field_param.fh"
      integer rtdb

      character*1 w
      character*2 tp(max_nfield), is
      logical p_srfc, p_rate, p_yild, p_avbeta, p_avgamma, p_tmd
      logical p_tmd_atom, need_dc, need_td, r
      integer nf, i, tg, nt
      double precision mx(max_nfield), fr(max_nfield), ce(max_nfield)
      double precision wi(max_nfield), cp(max_nfield), dt, tmax

      parameter(w='i')


c     ==== Get print switches ====
      call wfat_prep_printsw(rtdb, 'inp:', p_srfc, p_rate, p_yild,
     &     p_avbeta, p_avgamma, p_tmd, p_tmd_atom)
      need_dc = (p_rate .or. p_tmd .or. p_tmd_atom)
      need_td = (p_yild .or. p_avbeta .or. p_avgamma)

      
c     ==== Get the number of field blocks ====
      if (.not. rtdb_get(rtdb, 'wfat:field:n', mt_int, 1, nf)) nf = 0


c     ==== Check for consistencies with the print directive ====
      if (nf > 0) then
         if (.not. rtdb_get(rtdb, 'wfat:field:typeg', mt_int, 1, tg))
     &        call errquit('WFAT: Failed reading field:typeg from rtdb',
     &        1, RTDB_ERR)
         if (need_dc .and. tg /= 1)
     &        call errquit('WFAT: When rate or tmd is requested, the '//
     &        'field must be DC.', 30, UERR)
         if (need_td .and. tg == 1)
     &        call errquit('WFAT: When yield, averaged or not, is '//
     &        'requested, the field must not be DC.', 30, UERR)
         if (.not. need_dc .and. .not. need_td)
     &        write(luout, 201)
      else if (nf == 0) then
         if (need_dc) call errquit('WFAT: One field block of type dc '//
     &        'is needed when the rate and/or tmd are specified for '//
     &        'the print directive.', 90, INPUT_ERR)
         if (need_td) call errquit('WFAT: At least one time-dependent'//
     &        ' field block (gaussian, hann, or one from an external '//
     &        'file) is needed when the yields, averaged or not, are '//
     &        'specified for the print directive.', 90, INPUT_ERR)
      endif
      

 201  format(/,
     &      3x, '*** ATTENTION ***', /, 
     &      3x, 'One or more field blocks are found but ', 
     &      'will be ignored because the input directive', /, 
     &      3x, '''print'' does not specify rate, yield, or', 
     &      'TMD.')      


c     ==== Fill local vectors ====
      do i = 1, nf
         write(is, '(i2)') i
         is = adjustl(is)

c        >>>> Field type <<<<
         if (.not. rtdb_cget(rtdb, 'wfat:inp:field:type:'//is, 1, 
     &        tp(i))) call errquit('WFAT: Failed reading '//
     &        'inp:field:type:'//is//' from rtdb', 1, RTDB_ERR)

         if (tg==1 .or. tg==3) then
c           >>>> Amplitude <<<<
            if (.not. rtdb_get(rtdb, 'wfat:inp:field:max:'//is, mt_dbl, 
     &           1, mx(i))) call errquit('WFAT: Failed reading '//
     &           'inp:field:max:'//is//' from rtdb', 1, RTDB_ERR)
            call wfat_dval_atleast(mx(i), 0.0d0, 'q', 'WFAT',
     &           'field '//is//' max', 'e', w, r)
         endif
         
         if (tg==3) then
c           >>>> Frequency <<<<
            if (.not. rtdb_get(rtdb, 'wfat:inp:field:freq:'//is, mt_dbl,
     &           1, fr(i))) call errquit('WFAT: Failed reading '//
     &           'inp:field:freq:'//is//' from rtdb', 1, RTDB_ERR)
            call wfat_dval_atleast(fr(i), 0.0d0, 'q', 'WFAT', 
     &           'field '//is//' frequency', 'e', w, r)
c           >>>> Center <<<<
            if (.not. rtdb_get(rtdb, 'wfat:inp:field:center:'//is, 
     &           mt_dbl, 1, ce(i))) call errquit('WFAT: Failed '//
     &           'reading inp:field:center:'//is//' from rtdb',
     &           1, RTDB_ERR)
            call wfat_dval_atleast(ce(i), 0.0d0, 'q', 'WFAT',
     &           'field '//is//' center', 'e', w, r)
c           >>>> Width <<<<
            if (.not. rtdb_get(rtdb, 'wfat:inp:field:width:'//is, 
     &           mt_dbl, 1, wi(i))) call errquit('WFAT: Failed '//
     &           'reading inp:field:width:'//is//' from rtdb',
     &           1, RTDB_ERR)
            call wfat_dval_atleast(wi(i), 0.0d0, 'g', 'WFAT',
     &           'field '//is//' width', 'e', w, r)
c           >>>> CEP <<<<
            if (.not. rtdb_get(rtdb, 'wfat:inp:field:cep:'//is, mt_dbl, 
     &           1, cp(i))) call errquit('WFAT: Failed reading '//
     &           'inp:field:cep:'//is//' from rtdb', 1, RTDB_ERR)
         endif
      enddo


c     ==== Write the above local vectors to RTDB ====
c     >>>> Field type <<<<
      if (.not. rtdb_cput(rtdb, 'wfat:field:type', nf, tp(1:nf)))
     &     call errquit('WFAT: Failed writing field:type to rtdb',
     &     1, RTDB_ERR)
      
      if (tg==1 .or. tg==3) then
c        >>>> Amplitude <<<<
         if (.not. rtdb_put(rtdb, 'wfat:field:max', mt_dbl, nf,
     &        mx(1:nf))) call errquit('WFAT: Failed writing '//
     &        'field:max to rtdb', 1, RTDB_ERR)
      endif

      if (tg==3) then
c        >>>> Frequency <<<<
         if (.not. rtdb_put(rtdb, 'wfat:field:freq', mt_dbl, nf,
     &        fr(1:nf))) call errquit('WFAT: Failed writing '//
     &        'field:freq to rtdb', 1, RTDB_ERR)
c        >>>> Center <<<<
         if (.not. rtdb_put(rtdb, 'wfat:field:center', mt_dbl, nf,
     &        ce(1:nf))) call errquit('WFAT: Failed writing '//
     &        'field:center to rtdb', 1, RTDB_ERR)
c        >>>> Width <<<<
         if (.not. rtdb_put(rtdb, 'wfat:field:width', mt_dbl, nf,
     &        wi(1:nf))) call errquit('WFAT: Failed writing '//
     &        'field:width to rtdb', 1, RTDB_ERR)
c        >>>> CEP <<<<
         if (.not. rtdb_put(rtdb, 'wfat:field:cep', mt_dbl, nf,
     &        cp(1:nf))) call errquit('WFAT: Failed writing '//
     &        'field:cep to rtdb', 1, RTDB_ERR)
      endif
      

      if (tg==2 .or. tg==3) then
c        >>>> dt <<<<
         call wfat_process_ddata(rtdb, 'field:dt', .true., 0.0d0, dt)
         call wfat_dval_atleast(dt, 0.0d0, 'g', 'WFAT', '''dt''', 
     &        'e', w, r)         
c        >>>> tmax <<<<
         call wfat_process_ddata(rtdb, 'field:tmax', .true., 0.0d0,
     &        tmax)
         call wfat_dval_atleast(tmax, 0.0d0, 'g', 'WFAT', '''tmax''',
     &        'e', w, r)
      endif
         

      if (tg==3) then
         nt = int(tmax/dt) + 1
         if (.not. rtdb_put(rtdb, 'wfat:field:ntime', mt_int, 1, nt))
     &        call errquit('WFAT: Failed writing field:ntime to rtdb',
     &        1, RTDB_ERR)
      endif
      
      end
c=======================================================
