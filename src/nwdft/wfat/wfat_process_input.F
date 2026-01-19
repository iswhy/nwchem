c=======================================================
c=======================================================
      subroutine wfat_process_input(rtdb)

      implicit none
#include "errquit.fh"
#include "wfat_input.fh"
#include "rtdb.fh"
#include "stdio.fh"      
      integer rtdb
      

c     ==== Obtain the level of WFAT (OE or ME) ====
      if (.not. rtdb_cget(rtdb,'task:theory', 1, wfat_level)) 
     &     call errquit('WFAT: Failed reading theory', 0, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'wfat:level', 1, trim(wfat_level)))
     &     call errquit('WFAT: Error writing waft_level to RTDB',
     &     1, RTDB_ERR)


c     ==== Check the existence of inputs ====
      call wfat_exist_input(rtdb)

c     ==== Check the sanity of inputs ====
      call wfat_sanity_input(rtdb)

c     ==== Put the checked inputs to the RTDB ====
      call wfat_putrtdb_input(rtdb)

c     ==== Process inputs using the new algorithm ====
      call wfat_process_input0(rtdb)      
            
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_exist_input(rtdb)

      implicit none
#include "errquit.fh"
#include "wfat_input.fh"
#include "stdio.fh"      
      integer rtdb

      logical print_srfc, print_rate, print_yild
      logical print_avbeta, print_avgamma, print_tmd, print_tmd_atom

      
      call wfat_prep_printsw(rtdb, 'inp:', print_srfc, print_rate, 
     &     print_yild, print_avbeta, print_avgamma, print_tmd,
     &     print_tmd_atom)
cOLD      call util_get_rtdb_state(rtdb, os, oc, or)


c     ==========================================================
c     ==================== EXISTENCE CHECKS ====================
c     ==========================================================
      
c     ==== Checking existence of CONDITIONALLY OPTIONAL INPUTS ====
c     == nktmd ==
      if (.not. nktmd_exist .and. (print_tmd .or. print_tmd_atom)) then
         nktmd = def_nktmd
         nktmd_stat = .true.
      else if (nktmd_exist .and. .not. (print_tmd .or. print_tmd_atom))
     &        then
         write(luout, 201) 'nktmd',
     &        '* print_tmd       ',
     &        '* print_tmd_atom  ', ' ', ' ', ' ', ' ', ' '
      endif
      
c     ==== Checking existence of OPTIONAL INPUTS ====
      if (.not. omgpos_exist) then
         autopos = .true.
         omgpos = 0.0d0
         omgpos_stat = .true.
      endif
      if (.not. wfat1_exist) then
         wfat1 = DEF_WFAT1
         wfat1_stat = .true.
      endif
      if (.not. wfat_mask_exist) then
         wfat_mask = DEF_WFAT_MASK
         wfat_mask_param = (/ -1.0d6, -1.0d6 /)
         wfat_mask_stat = .true.
      endif
      if (wfat_level == 'oe') then
         if (.not. movec_exist) then
            call wfat_def_movecname('oe', rtdb, wfat_movec_name)
            movec_stat = .true.
         endif
      else if (wfat_level == 'me') then
         if (.not. movec_exist) then
            call wfat_def_movecname('me0', rtdb, wfat_movec_name)
            movec_stat = .true.
         endif
         if (.not. movec1_exist) then
            call wfat_def_movecname('me1', rtdb, wfat_movec1_name)
            movec1_stat = .true.
         endif
         if (wfat1) call errquit('MEWFAT: At the moment, MEWFAT '//
     &        'cannot be used with the first order correction. '//
     &        'Please remove the ''wfat1'' keyword from your input '//
     &        'file.', 62, UERR)
      endif
      if (.not. prtb_type_exist) then
         prtb_type = DEF_PRTB_type
         prtb_type_stat = .true.
      endif
      if (.not. prtb_ord_exist) then
         if (prtb_type=='bw') prtb_ord = DEF_PRTB_ORD_BW
         if (prtb_type=='rs') prtb_ord = DEF_PRTB_ORD_RS
         prtb_ord_stat = .true.
      endif
      if (.not. prtb_bt_exist) then
         prtb_bt = DEF_PRTB_BT
         prtb_bt_stat = .true.
      endif
      if (.not. prtb_gm_exist) then
         prtb_gm = DEF_PRTB_GM
         prtb_gm_stat = .true.
      endif
      if (.not. prtb_dgnthr_exist) then
         prtb_dgnthr = DEF_PRTB_DGNTHR
         prtb_dgnthr_stat = .true.
      endif
      if (.not. prtb_offthr_exist) then
         prtb_offthr = DEF_PRTB_OFFTHR
         prtb_offthr_stat = .true.
      endif
      if (.not. prtb_e1thr_exist) then
         prtb_e1thr = DEF_PRTB_E1THR
         prtb_e1thr_stat = .true.
      endif
      if (.not. prtb_itrthr_exist) then
         prtb_itrthr = DEF_PRTB_ITRTHR
         prtb_itrthr_stat = .true.
      endif
      if (.not. prtb_maxitr_exist) then
         prtb_maxitr = DEF_PRTB_MAXITR
         prtb_maxitr_stat = .true.
      endif
      if (.not. prtb_ordthr_exist) then
         prtb_ordthr = DEF_PRTB_ORDTHR
         prtb_ordthr_stat = .true.
      endif
      if (.not. prtb_maxord_exist) then
         prtb_maxord = DEF_PRTB_MAXORD
         prtb_maxord_stat = .true.
      endif
      
      
c     =======================================================
c     ==================== FORMAT CHECKS ====================
c     =======================================================
            
c     ==== Checking format of CONDITIONALLY OPTIONAL INPUTS ====
      if (.not. nktmd_stat .and. (print_tmd .or. print_tmd_atom)) then
         call errquit('WFAT: The format of input directive '//
     &        '''nktmd'' is incorrect. The correct format '//
     &        'should be ''nktmd <integer value>''.',
     &        0, INPUT_ERR)
      endif

c     ==== Checking format of OPTIONAL INPUTS ====
      if (.not. omgpos_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''omgpos''. Please check if it is specified '//
     &        'following the format ''omgpos <string ''auto''> || '//
     &        '<float valuex  float valuey  float valuez>''',
     &        0, INPUT_ERR)
      endif
      if (wfat_level == 'oe') then
         if (.not. movec_stat) then
            call errquit('WFAT: Error reading optional input '//
     &           'parameter ''vectors''. Please check if it is '//
     &           'specified following the format '//
     &           '''vectors "<string value>"'' '//
     &           '(do not skip the double quotes).',
     &           0, INPUT_ERR)
         endif
      else if (wfat_level == 'me') then
         if (.not. (movec_stat .and. movec1_stat)) then
            call errquit('WFAT: Error reading optional input '//
     &           'parameter ''vectors''. Please check if it is '//
     &           'specified following the format '//
     &           '''vectors "<string value1>" ["<string value2>"]'' '//
     &           '(do not skip the double quotes).',
     &           0, INPUT_ERR)
         endif
      endif
      if (.not. wfat1_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''wfat1''. Please check if it is specified following '//
     &        'the format ''wfat1''.',
     &        0, INPUT_ERR)
      endif

      if (.not. wfat_mask_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''tailmask''. Please check if it is specified '//
     &        'following the format ''tailmask <float value> '//
     &        '<float value>''.',
     &        0, INPUT_ERR)
      endif
      if (.not. prtb_type_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        'that controls the type of perturbation analysis. '//
     &        'Please check if it is correctly specified inside '//
     &        'the ''perturb'' directive. The valid option is '//
     &        '''bw'' or ''rs''.', 0, INPUT_ERR)
      endif
      if (.not. prtb_ord_stat) then
         write(luout, 301) 'order'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_bt_stat) then
         write(luout, 301) 'beta'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_gm_stat) then
         write(luout, 301) 'gamma'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_dgnthr_stat) then
         write(luout, 301) 'dgnthr'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_offthr_stat) then
         write(luout, 301) 'offthr'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_e1thr_stat) then
         write(luout, 301) 'e1thr'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_itrthr_stat) then
         write(luout, 301) 'itrthr'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_maxitr_stat) then
         write(luout, 301) 'maxitr'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_ordthr_stat) then
         write(luout, 301) 'ordthr'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif
      if (.not. prtb_maxord_stat) then
         write(luout, 301) 'maxord'
         call errquit('WFAT: Input error', 0, INPUT_ERR)
      endif

      
 301  format(/
     &       'ERROR: Cannot read optional input parameter ''', a,
     &       ''' in ''perturb'' directive. ', /,
     &       'Please check if it is correctly specified inside ',
     &       'the ''perturb'' directive.')

 302  format(/
     &     'ERROR: In the ''print'' directive, suppose A={',
     &     'print_rate, print_tmd, print_tmd_atom} ', /,
     &     7x, 'and B={print_yild, print_avbeta, print_avgamma}, ',
     &     'then only the members of one set (A or B) ', /,
     &     7x, 'may be specified at a time. For example, if both ',
     &     '''print_avbeta'' and ''print_rate'' are ', /,
     &     7x, 'specified, this error will happen. The reason of ',
     &     'this behavior is that the tasks in A run ', /,
     &     7x, 'on a DC field, whereas those in B run on a time-',
     &     'dependent field, and a given job can only ',/
     &     7x, 'perform either DC or time-dependent calculation.')
      
      
c     ====================================================
c     ================= NECESSITY CHECKS =================
c     ====================================================
      if (wfat_level == 'oe') then
         call wfat_need_input(rtdb)
      else if (wfat_level == 'me') then 
cdebug         call mewfat_dummy()
         call mewfat_need_input(rtdb)
      endif


 201  format(/,
     &      3x, '*** ATTENTION ***', /, 
     &      3x, 'Input directive ''', a,''' is found but ', 
     &      'will be ignored because the input directive', /, 
     &      3x, '''print'' does not contain any of the ', 
     &      'following: ', /,
     &      5x, a18, /,
     &      5x, a18, /,
     &      5x, a18, /,
     &      5x, a18, /,
     &      5x, a18, /,
     &      5x, a18, /,
     &      5x, a18)

 203  format(/,
     &     3x, '*** ATTENTION ***', /
     &     3x, 'Input parameter that controls the type of ',
     &     'perturbation is found but will be ignored because the ', /
     &     3x, 'present calculation is a restarted job in which case ',
     &     'the value of this parameter is obtained from RTDB.')      
      
 204  format(/,
     &     3x, '*** ATTENTION ***', /
     &     3x, 'Input parameter ''', a, ''' belonging to ''perturb'' ',
     &     'directive is found but will be ignored because the ', /
     &     3x, 'present calculation is a restarted job in which case ',
     &     'the value of this parameter is obtained from RTDB.')      

 205  format(/,
     &      'WFAT: Cannot find input directive ''', a,'''. ', 
     &      'This directive is required when the input directive ', /, 
     &      6x, '''print'' contains any of the following keywords: ', /,
     &      5x, a, /,
     &      5x, a, /,
     &      5x, a, /,
     &      5x, a, /, 
     &      5x, a, /,
     &      5x, a, /,
     &      5x, a)

 206  format(3x, '*** ATTENTION ***', /
     &     3x, 'Input directive ''', a, ''' is found ',
     &     'but will be ignored because either DC fields are being ',
     &     'used or no fields are specified.')
      
      end
c=======================================================


cc c=======================================================
cc c=======================================================
cc       subroutine mewfat_need_input()
cc       implicit none
cc #include "errquit.fh"
cc 
cc       call errquit('stopping in mewfat_need_input', 26, UERR)
cc       end
cc c=======================================================

      
c=======================================================
c=======================================================
      subroutine mewfat_need_input(rtdb)

      implicit none
#include "wfat_input.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer rtdb
      
      logical sr, rt, yi, avb, avg, tmd, tmda
      character*1 m


      nktmd_need = .false.          ! needed when OE or ME
      omgpos_need = .false.         ! needed when OE or ME
      lmax_need = .false.           ! needed when OE or ME
      movec_need = .false.          ! needed when OE or ME
      movec1_need = .false.         ! needed when ME
      wfat1_need = .false.          ! needed when OE or ME
      wfat_mask_need = .false.      ! needed when oe or me
      prtb_type_need = .false.      ! needed when OE
      prtb_ord_need = .false.       ! needed when OE
      prtb_bt_need = .false.        ! needed when OE
      prtb_gm_need = .false.        ! needed when OE
      prtb_dgnthr_need = .false.    ! needed when OE
      prtb_offthr_need = .false.    ! needed when OE
      prtb_e1thr_need = .false.     ! needed when OE
      prtb_itrthr_need = .false.    ! needed when OE
      prtb_maxitr_need = .false.    ! needed when OE
      prtb_ordthr_need = .false.    ! needed when OE
      prtb_maxord_need = .false.    ! needed when OE


      call wfat_prep_printsw(rtdb, 'inp:', sr, rt, yi, avb, avg, tmd, 
     &     tmda)
      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
     &     .false., .false., nktmd_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     omgpos_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     lmax_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat_mask_need)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_need_input(rtdb)

      implicit none
#include "wfat_input.fh"
#include "stdio.fh"
      integer rtdb
      
      logical sr, rt, yi, avb, avg, tmd, tmda
      character*1 m

      
      nktmd_need = .false.              ! needed when oe or me
      omgpos_need = .false.             ! needed when oe or me
      lmax_need = .false.               ! needed when oe or me
      movec_need = .false.              ! needed when oe or me
      movec1_need = .false.             ! needed when me
      wfat1_need = .false.              ! needed when oe or me
      wfat_mask_need = .false.          ! needed when oe or me
      prtb_type_need = .false.      ! needed when oe
      prtb_ord_need = .false.       ! needed when oe
      prtb_bt_need = .false.        ! needed when oe
      prtb_gm_need = .false.        ! needed when oe
      prtb_dgnthr_need = .false.    ! needed when oe
      prtb_offthr_need = .false.    ! needed when oe
      prtb_e1thr_need = .false.     ! needed when oe
      prtb_itrthr_need = .false.    ! needed when oe
      prtb_maxitr_need = .false.    ! needed when oe
      prtb_ordthr_need = .false.    ! needed when oe
      prtb_maxord_need = .false.    ! needed when oe


      call wfat_prep_printsw(rtdb, 'inp:', sr, rt, yi, avb, avg, tmd, 
     &     tmda)      
      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
     &     .false., .false., nktmd_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     omgpos_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     lmax_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat_mask_need)
      call wfat_need_input0(rt, tmd, tmda, .false., .false., .false., 
     &     .false., prtb_type_need)
      call wfat_need_input0(rt, tmd, tmda, .false., .false., .false., 
     &     .false., prtb_ord_need)
      call wfat_need_input0(rt, tmd, tmda, .false., .false., .false., 
     &     .false., prtb_bt_need)
      call wfat_need_input0(rt, tmd, tmda, .false., .false., .false., 
     &     .false., prtb_gm_need)
      call wfat_need_input0(rt, tmd, tmda, .false., .false., .false., 
     &     .false., prtb_dgnthr_need)
      call wfat_need_input0(rt, tmd, tmda, .false., .false., .false., 
     &     .false., prtb_offthr_need)
      if (prtb_type == 'bw') then
         call wfat_need_input0(rt, tmd, tmda, .false., .false., .false.,
     &        .false., prtb_e1thr_need)
         call wfat_need_input0(rt, tmd, tmda, .false., .false., .false.,
     &        .false., prtb_itrthr_need)
         call wfat_need_input0(rt, tmd, tmda, .false., .false., .false.,
     &        .false., prtb_maxitr_need)
         call wfat_need_input0(rt, tmd, tmda, .false., .false., .false.,
     &        .false., prtb_ordthr_need)
         call wfat_need_input0(rt, tmd, tmda, .false., .false., .false.,
     &        .false., prtb_maxord_need)
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_need_input0(a, b, c, d, e, f, g, out)

      implicit none
      logical a, b, c, d, e, f, g
      logical out

      if (a .or. b .or. c .or. d .or. e .or. f .or. g) then
         out = .true.
      else
         out = .false.
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_sanity_input(rtdb)

      implicit none
#include "wfat_input.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"      
      integer rtdb

      integer i, is
      double precision dt_00, tmax_00
      
      character*1 w
      logical r
      parameter ( w='i' )
      integer k

      logical print_srfc, print_rate, print_yild, print_avbeta
      logical print_avgamma, print_tmd, print_tmd_atom

      
      call wfat_prep_printsw(rtdb, 'inp:', print_srfc, print_rate, 
     &     print_yild, print_avbeta, print_avgamma, print_tmd,
     &     print_tmd_atom)
            
      if (nktmd_need) then
         call wfat_ival_atleast(nktmd, 1, 'WFAT', '''nktmd''', w, r)
      endif
      
      if (wfat_mask_need .and. wfat_mask) then
         call wfat_dval_atleast(wfat_mask_param(1), 0.0d0, 'g', 'WFAT',
     &        '''tailmask:param_1''', 'e', w, r)
         call wfat_dval_atleast(wfat_mask_param(2), 0.0d0, 'g', 'WFAT',
     &        '''tailmask:param_2''', 'e', w, r)
      endif
      
      if (prtb_ord_need) then
         if (prtb_type == 'rs') then
            call wfat_ival_inrange(prtb_ord, 1, 1, 'WFAT', 
     &           '''perturb:order''', 'o', r)
            if (.not. r)
     &           call errquit('WFAT: When ''rs'' perturbation type '//
     &           'is chosen, the value of input parameter '//
     &           '''perturb:order'' can only be 1.', 1, INPUT_ERR)
         else if (prtb_type == 'bw') then
            call wfat_ival_inrange(prtb_ord, 1, prtb_maxord, 'WFAT', 
     &           '''perturb:order''', w, r)
         endif
      endif

      if (prtb_bt_need) then
         call wfat_dval_inrange(prtb_bt, 0.0d0, 'q', 180.0d0, 'q',
     &        'WFAT', '''perturb:beta''', 'f', w, r)
      endif

      if (prtb_gm_need) then
         call wfat_dval_inrange(prtb_gm, 0.0d0, 'q', 360.0d0, 'q',
     &        'WFAT', '''perturb:beta''', 'f', w, r)
      endif

      if (prtb_dgnthr_need) then
         call wfat_dval_atleast(prtb_dgnthr, 0.0d0, 'q', 'WFAT',
     &        '''perturb:dgnthr''', 'e', w, r)
      endif

      if (prtb_offthr_need) then
         call wfat_dval_atleast(prtb_offthr, 0.0d0, 'q', 'WFAT',
     &        '''perturb:offthr''', 'e', w, r)
      endif

      if (prtb_e1thr_need) then
         call wfat_dval_atleast(prtb_e1thr, 0.0d0, 'q', 'WFAT',
     &        '''perturb:e1thr''', 'e', w, r)
      endif

      if (prtb_itrthr_need) then
         call wfat_dval_atleast(prtb_itrthr, 0.0d0, 'q', 'WFAT',
     &        '''perturb:itrthr''', 'e', w, r)
      endif

      if (prtb_maxitr_need) then
         call wfat_ival_atleast(prtb_maxitr, 10, 'WFAT',
     &        '''perturb:maxitr''', w, r)
      endif

      if (prtb_ordthr_need) then
         call wfat_dval_atleast(prtb_ordthr, 0.0d0, 'q', 'WFAT',
     &        '''perturb:ordthr''', 'e', w, r)
      endif

      if (prtb_maxord_need) then
         call wfat_ival_atleast(prtb_maxord, 1, 'WFAT',
     &        '''perturb:maxord''', w, r)
      endif
      
cOLD 500  format(/,
cOLD     &     3x, '*** ATTENTION ***', /,
cOLD     &     3x, 'NWChem detected that the following situations ',
cOLD     &     'are true:', /,
cOLD     &     3x, '  * ionization rate is requested, ', /,
cOLD     &     3x, '  * there are more than one channels, and ', /,
cOLD     &     3x, '  * unit field factor is chosen.', /,
cOLD     &     3x, 'While the above situations do not lead to ',
cOLD     &     'a runtime error, the calculated total ionization ', /,
cOLD     &     3x, 'rate is erroneous. This is because the rate of ',
cOLD     &     'several channels are added together to compute total ', /,
cOLD     &     3x, 'rate. And each of these rates has a specific ',
cOLD     &     'field factor that depends on their orbital energies. ', /,
cOLD     &     3x, 'Summing these rates while some or all of them have ',
cOLD     &     'their corresponding field factor artificially set to ', /,
cOLD     &     3x, 'unity produces an unphysical situation.')
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_putrtdb_input(rtdb)

      implicit none
#include "stdio.fh"      
#include "rtdb.fh"
#include "errquit.fh"
#include "wfat_input.fh"
#include "mafdecls.fh"      
      integer rtdb

      character*2 td
      logical nosave
      integer nt_00, er(2)
      double precision dt_00, tmax_00


cOLD      call util_get_rtdb_state(rtdb, os, oc, or)
cOLDc     .and. os means write that input into RTDB only when it is a start-type job.
cOLDc     An example of input that will only be written into the RTDB when doing a
cOLDc     start-type job (os=.true.) is lmax. If it is written into the RTDB when
cOLDc     the job is not of start-type, meaning restart-type, it will overwrite lmax
cOLDc     that is already in the RTDB.

      

c forcheck      if (.not. rtdb_print(rtdb, .false.))
c forcheck     &     call errquit('WFAT: Printing rtdb failed', 1, RTDB_ERR)

      
      if (nktmd_need) then
         if (.not. rtdb_put(rtdb, 'wfat:nktmd', mt_int, 1, nktmd))
     &        call errquit('WFAT: Failed writing nktmd to rtdb',
     &        1, INPUT_ERR)
      endif
      
      if (omgpos_need) then
         if (.not. rtdb_put(rtdb, 'wfat:autopos', mt_log, 1, autopos))
     &        call errquit('WFAT: Failed writing autopos to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:omgpos', mt_dbl, 3, omgpos))
     &        call errquit('WFAT: Failed writing omgpos to rtdb',
     &        1, INPUT_ERR)
      endif

      if (wfat_level == 'oe') then
         if (movec_need) then
            if (.not. rtdb_cput(rtdb, 'wfat:movec_name', 1,
     &           trim(adjustl(wfat_movec_name))))
     &           call errquit('WFAT: Failed writing movec_name to RTDB',
     &           1, RTDB_ERR)
         endif
      else if (wfat_level == 'me') then
         if (movec_need) then
            if (.not. rtdb_cput(rtdb, 'mewfat:movec0_name', 1,
     &           trim(adjustl(wfat_movec_name))))
     &           call errquit('WFAT: Failed writing movec0_name to '//
     &           'RTDB', 1, RTDB_ERR)
         endif
         if (movec1_need) then
            if (.not. rtdb_cput(rtdb, 'mewfat:movec1_name', 1,
     &           trim(adjustl(wfat_movec1_name))))
     &           call errquit('WFAT: Failed writing movec1_name to '//
     &           'RTDB', 1, RTDB_ERR)
         endif
      endif

cOLD      if (prtb_type_need .and. os) then
      if (prtb_type_need) then
         if (.not. rtdb_cput(rtdb, 'wfat:perturb_type', 1, prtb_type))
     &        call errquit('WFAT: Failed writing perturb_type to rtdb',
     &        1, INPUT_ERR)
      endif

cOLD      if (prtb_ord_need .and. os) then
      if (prtb_ord_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_order', mt_int, 1,
     &        prtb_ord))
     &        call errquit('WFAT: Failed writing perturb_order to rtdb',
     &        1, INPUT_ERR)
      endif

      if (prtb_bt_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_beta', mt_dbl, 1,
     &        prtb_bt))
     &        call errquit('WFAT: Failed writing perturb_beta to rtdb',
     &        1, INPUT_ERR)
      endif
      
      if (prtb_gm_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_gamma', mt_dbl, 1,
     &        prtb_gm))
     &        call errquit('WFAT: Failed writing perturb_gamma to rtdb',
     &        1, INPUT_ERR)
      endif

      if (prtb_dgnthr_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_dgnthr', mt_dbl, 1,
     &        prtb_dgnthr))
     &        call errquit('WFAT: Failed writing perturb_dgnthr to rtdb'
     &        , 1, INPUT_ERR)
      endif

      if (prtb_offthr_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_offthr', mt_dbl, 1,
     &        prtb_offthr))
     &        call errquit('WFAT: Failed writing perturb_offthr to rtdb'
     &        , 1, INPUT_ERR)
      endif

      if (prtb_e1thr_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_e1thr', mt_dbl, 1,
     &        prtb_e1thr))
     &        call errquit('WFAT: Failed writing perturb_e1thr to rtdb'
     &        , 1, INPUT_ERR)
      endif

      if (prtb_itrthr_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_itrthr', mt_dbl, 1,
     &        prtb_itrthr))
     &        call errquit('WFAT: Failed writing perturb_itrthr to rtdb'
     &        , 1, INPUT_ERR)
      endif

      if (prtb_maxitr_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_maxitr', mt_int, 1,
     &        prtb_maxitr))
     &        call errquit('WFAT: Failed writing perturb_maxitr to rtdb'
     &        , 1, INPUT_ERR)
      endif

      if (prtb_ordthr_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_ordthr', mt_dbl, 1,
     &        prtb_ordthr))
     &        call errquit('WFAT: Failed writing perturb_ordthr to rtdb'
     &        , 1, INPUT_ERR)
      endif

      if (prtb_maxord_need) then
         if (.not. rtdb_put(rtdb, 'wfat:perturb_maxord', mt_int, 1,
     &        prtb_maxord))
     &        call errquit('WFAT: Failed writing perturb_maxord to rtdb'
     &        , 1, INPUT_ERR)
      endif

      if (wfat1_need) then
         if (.not. rtdb_put(rtdb, 'wfat:wfat1', mt_log, 1, wfat1))
     &        call errquit('WFAT: Failed writing wfat1 to RTDB', 1,
     &        RTDB_ERR)
      endif

      if (wfat_mask_need) then
         if (.not. rtdb_put(rtdb, 'wfat:mask', mt_log, 1, wfat_mask))
     &        call errquit('WFAT: Failed writing mask to RTDB', 1,
     &        RTDB_ERR)
         if (wfat_mask) then
            if (.not. rtdb_put(rtdb, 'wfat:mask_g', mt_dbl, 2,
     &           wfat_mask_param))
     &           call errquit('WFAT: Failed writing mask_g to RTDB', 1,
     &           RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:mask_j', mt_dbl, 2,
     &           wfat_mask_param))
     &           call errquit('WFAT: Failed writing mask_j to RTDB', 1,
     &           RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:mask_k', mt_dbl, 2,
     &           wfat_mask_param))
     &           call errquit('WFAT: Failed writing mask_k to RTDB', 1,
     &           RTDB_ERR)
         endif
      endif

      if (s2a1e) then
         if (.not. rtdb_put(rtdb, 'wfat:2a1e_nc', mt_int, 1, wfat2a_nc))
     &        call errquit('WFAT: Failed writing 2a1e_nc to RTDB', 1,
     &        RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:2a1e_q', mt_dbl, 2, wfat2a_q))
     &        call errquit('WFAT: Failed writing 2a1e_q to RTDB', 1,
     &        RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:internuc', mt_dbl, 1,
     &        wfat2a_inter))
     &        call errquit('WFAT: Failed writing internuc to RTDB', 1,
     &        RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:2a1e_m', mt_int, 1, wfat2a_mag))
     &        call errquit('WFAT: Failed writing 2a1e_m to RTDB', 1,
     &        RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:2a1e_xz', mt_int, 1,
     &        wfat2a_xzr))
     &        call errquit('WFAT: Failed writing 2a1e_xz to RTDB', 1,
     &        RTDB_ERR)
         if (.not. rtdb_cput(rtdb, 'wfat:2a1e_pot', 1, wfat2a_pot))
     &        call errquit('WFAT: Failed writing 2a1e_pot to RTDB', 1,
     &        RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:2a1e_potparam', mt_dbl, 2,
     &        wfat2a_potparam))
     &        call errquit('WFAT: Failed writing 2a1e_potparam to RTDB',
     &        1, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:2a1e_nel', mt_int, 1,
     &        wfat2a_nel))
     &        call errquit('WFAT: Failed writing 2a1e_nel to RTDB', 1,
     &        RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:2a1e_nt', mt_int, 3, wfat2a_nt))
     &        call errquit('WFAT: Failed writing 2a1e_nt to RTDB', 1,
     &        RTDB_ERR)
         if (wfat2a_fedvr) then
            if (.not. rtdb_put(rtdb, 'wfat:fedvr', mt_log, 1,
     &           wfat2a_fedvr))
     &           call errquit('WFAT: Failed writing fedvr to RTDB', 1,
     &           RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:2a1e_nfe', mt_int, 1,
     &           wfat2a_nfe))
     &           call errquit('WFAT: Failed writing 2a1e_nfe to RTDB', 
     &           1, RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:2a1e_fenodes', mt_dbl,
     &           wfat2a_nfe, wfat2a_fenodes))
     &           call errquit('WFAT: Failed writing 2a1e_fenodes to '//
     &           'RTDB', 1, RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:2a1e_ndvr', mt_int,
     &           wfat2a_nfe, wfat2a_ndvr))
     &           call errquit('WFAT: Failed writing 2a1e_ndvr to '//
     &           'RTDB', 1, RTDB_ERR)
         else
            if (.not. rtdb_put(rtdb, 'wfat:2a1e_nr', mt_int, 3,
     &           wfat2a_nr))
     &           call errquit('WFAT: Failed writing 2a1e_nr to RTDB',
     &           1, RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:2a1e_inke', mt_dbl, 1, 
     &           wfat2a_kein))
     &           call errquit('WFAT: Failed writing 2a1e_inke to RTDB',
     &           1, RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:2a1e_multke', mt_dbl, 1, 
     &           wfat2a_kemult))
     &           call errquit('WFAT: Failed writing 2a1e_multke to '//
     &           'RTDB', 1, RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:2a1e_maxit', mt_int, 1,
     &           wfat2a_maxit))
     &           call errquit('WFAT: Failed writing 2a1e_maxit to RTDB',
     &           1, RTDB_ERR)
         endif
      endif      
      
      end
c=======================================================


