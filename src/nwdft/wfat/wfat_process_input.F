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
      if (.not. merg_exist) then
         wfat_merg(:) = DEF_MERG
         merg_stat = .true.
      else
         call wfat_process_myerg(ORBID_CHARLEN, wfat_merg_s, n_orbid,
     &        wfat_orbid0, wfat_spinid, wfat_merg, wfat_merg_l)
         merg_stat = .true.
      endif
      if (.not. wfat_mask_exist) then
         wfat_mask = DEF_WFAT_MASK
         wfat_mask_param = (/ -1.0d6, -1.0d6 /)
         wfat_mask_stat = .true.
      endif
      if (.not. ffc_exist) then
         ffc_type0(1:n_orbid) = DEF_FFC_TYPE
         ergffc(1:n_orbid) = DEF_ERG_FFC
         ffc_stat = .true.
      else if (ffc_exist .and. .not. print_rate) then
         write(luout, 201) 'field_factor',
     &        '* print_rate      ',
     &        ' ', ' ', ' ', ' ', ' ', ' '
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
      if (.not. orbid_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''movec_id''. Please check if it is '//
     &        'specified following the format '//
     &        '''movec_id a <integer so#1> [<integer so#2> ... ] '//
     &        '[b <integer so#1> [<integer so#2> ... ]]''. For '//
     &        'example, the following input ''movec_id a 0 -1 -2 '//
     &        'b -1'' corresponds to the following list of molecular '//
     &        'orbitals to be ionized: HOMO, HOMO-1, and HOMO-2 '//
     &        'having alpha-spin and HOMO-1 having beta-spin. Note: '//
     &        'At the moment only zero and negative integers are '//
     &        'supported for the spatial orbitals (so#1, so#2, etc.)',
     &        0, INPUT_ERR)
      endif
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

      if (.not. merg_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''myenergy''. Please check if it is specified '//
     &        'following the format ''myenergy '//
     &        '"spin, oid" erg ...'' where spin is ''a'' '//
     &        'or ''b'', oid is movec_id, and erg is the '//
     &        'manually given ionization energy.', 0, INPUT_ERR)
      endif
      
      if (.not. wfat_mask_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''tailmask''. Please check if it is specified '//
     &        'following the format ''tailmask <float value> '//
     &        '<float value>''.',
     &        0, INPUT_ERR)
      endif
      if (.not. ffc_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''field_factor''. Please check if it is specified '//
     &        'following the format ''field_factor <value> '//
     &        '<value> ... <value>'' where ''value'' is either a '//
     &        'floating point number, a string ''orb'', or a string '//
     &        '''unit''. There must be as many <value> fields as '//
     &        'there are the total number of ionized channels.',
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


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_acdc_field(nf, tp, tdf)
cOLD
cOLD      implicit none
cOLD#include "errquit.fh"
cOLD#include "wfat_input.fh"
cOLD      integer nf                         ! [input],
cOLD      character(len=FLD_TP_LEN) tp(nf)   ! [input],
cOLD      character*2 tdf                    ! [output],
cOLD
cOLD      integer i
cOLD
cOLD      tdf = 'no'
cOLD      do i = 1, nf
cOLD         if (tp(i) == 'dc') then
cOLD            tdf = 'ti'
cOLD         endif
cOLD
cOLD         if (tdf=='ti' .and. (tp(i)=='gaussian' .or. tp(i)=='hann' .or.
cOLD     &        tp(i)=='external'))
cOLD     &        call errquit('WFAT: When there is a ''dc'' type '//
cOLD     &        'field in the ''field'' directive, no other field '//
cOLD     &        'types are allowed to coexist', 1, INPUT_ERR)
cOLD         
cOLD         if (tp(i)=='gaussian' .or. tp(i)=='hann' .or.
cOLD     &        tp(i)=='external') then
cOLD            tdf = 'td'
cOLD         endif
cOLD      enddo
cOLD
cOLD      end
cOLDc=======================================================


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


cOLD      field_need = .false.          ! needed when OE or ME
cOLD      beta_tmd_need = .false.       ! needed when OE or ME
cOLD      gamma_tmd_need = .false.      ! needed when OE or ME
      nktmd_need = .false.          ! needed when OE or ME
      orbid_need = .false.          ! needed when OE
cOLD      gqual_need = .false.          ! needed when OE or ME
      omgpos_need = .false.         ! needed when OE or ME
cOLD      dgnthr_need = .false.         ! needed when oe or me
      lmax_need = .false.           ! needed when OE or ME
      movec_need = .false.          ! needed when OE or ME
      movec1_need = .false.         ! needed when ME
      wfat1_need = .false.          ! needed when OE or ME
      merg_need = .false.           ! needed when oe or me
cOLD      nlgl_need = .false.           ! needed when oe or me
      ffc_need = .false.            ! needed when oe or me
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
            
cOLD      call wfat_need_input0(rt, tmd, tmda, yi, avb, avg, .false.,
cOLD     &     field_need)
cOLD      call wfat_need_input0(yi, avb, avg, .false., .false., .false.,
cOLD     &     .false., dt_need)
cOLD      call wfat_need_input0(yi, avb, avg, .false., .false., .false.,
cOLD     &     .false., tmax_need)
cOLD      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
cOLD     &     .false., .false., beta_tmd_need)
cOLD      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
cOLD     &     .false., .false., gamma_tmd_need)
      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
     &     .false., .false., nktmd_need)

cOLD      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
cOLD     &     gqual_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     omgpos_need)
cOLD      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
cOLD     &     dgnthr_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     lmax_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     merg_need)
cOLD      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
cOLD     &     nlgl_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     ffc_need)
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

      
cOLD      field_need = .false.              ! needed when oe or me
cOLD      beta_tmd_need = .false.           ! needed when oe or me
cOLD      gamma_tmd_need = .false.          ! needed when oe or me
      nktmd_need = .false.              ! needed when oe or me
      orbid_need = .false.              ! needed when oe
cOLD      gqual_need = .false.              ! needed when oe or me
      omgpos_need = .false.             ! needed when oe or me
cOLD      dgnthr_need = .false.             ! needed when oe or me
      lmax_need = .false.               ! needed when oe or me
      movec_need = .false.              ! needed when oe or me
      movec1_need = .false.             ! needed when me
      wfat1_need = .false.              ! needed when oe or me
      merg_need = .false.               ! needed when oe or me
cOLD      nlgl_need = .false.               ! needed when oe or me
      ffc_need = .false.                ! needed when oe or me
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
      
cOLD      call wfat_need_input0(rt, tmd, tmda, yi, avb, avg, .false.,
cOLD     &     field_need)
cOLD      call wfat_need_input0(yi, avb, avg, .false., .false., .false.,
cOLD     &     .false., dt_need)
cOLD      call wfat_need_input0(yi, avb, avg, .false., .false., .false.,
cOLD     &     .false., tmax_need)
cOLD      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
cOLD     &     .false., .false., beta_tmd_need)
cOLD      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
cOLD     &     .false., .false., gamma_tmd_need)
      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
     &     .false., .false., nktmd_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     orbid_need)
cOLD      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
cOLD     &     gqual_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     omgpos_need)
cOLD      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
cOLD     &     dgnthr_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     lmax_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     merg_need)
cOLD      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
cOLD     &     nlgl_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     ffc_need)
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

      character*4 nch_s
      integer i, is, ffcount
      double precision dt_00, tmax_00
      
      character*1 w
      logical r, unit_ffc
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

      unit_ffc = .false.
      ffcount = 0
      do i = 1, n_orbid
         if (.not. unit_ffc) unit_ffc = (ffc_type0(i) == 1)
         if (ffc_type0(i)==1 .or. ffc_type0(i)==2 .or. ffc_type0(i)==3)
     &        ffcount = ffcount + 1
      enddo
      write(nch_s, '(i4)') n_orbid
      if (ffcount /= n_orbid)
     &     call errquit('WFAT: The number of entries given for the '//
     &     'input keyword ''field_factor'' must be the same as the '//
     &     'number of ionized channels, which is '//nch_s//'.',
     &     36, INPUT_ERR)
      if (ffc_need .and. unit_ffc .and. print_rate .and. n_orbid > 1)
     &     then
         write(luout, 500)
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
      
 500  format(/,
     &     3x, '*** ATTENTION ***', /,
     &     3x, 'NWChem detected that the following situations ',
     &     'are true:', /,
     &     3x, '  * ionization rate is requested, ', /,
     &     3x, '  * there are more than one channels, and ', /,
     &     3x, '  * unit field factor is chosen.', /,
     &     3x, 'While the above situations do not lead to ',
     &     'a runtime error, the calculated total ionization ', /,
     &     3x, 'rate is erroneous. This is because the rate of ',
     &     'several channels are added together to compute total ', /,
     &     3x, 'rate. And each of these rates has a specific ',
     &     'field factor that depends on their orbital energies. ', /,
     &     3x, 'Summing these rates while some or all of them have ',
     &     'their corresponding field factor artificially set to ', /,
     &     3x, 'unity produces an unphysical situation.')
      end
c=======================================================


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_sanity_field(nf, tp, fr, ce, wi, mx)
cOLD
cOLD      implicit none
cOLD#include "wfat_input.fh"
cOLD      integer nf
cOLD      character(len=FLD_TP_LEN) tp(nf)
cOLD      double precision fr(max_nfield)
cOLD      double precision ce(max_nfield)
cOLD      double precision wi(max_nfield)
cOLD      double precision mx(max_nfield)
cOLD
cOLD      character*2 s
cOLD      character*3 is
cOLD      character*1 w
cOLD      logical r
cOLD      integer i
cOLD      parameter ( w='i' )
cOLD      
cOLD      
cOLD      do i = 1, nf
cOLD         write(is, '(i3)') i
cOLD         call wfat_ordinal_suffix(i, s)
cOLD         if (tp(i) /= 'dc' .and. tp(i) /= 'external') then
cOLD            call wfat_dval_atleast(fr(i), 0.0d0, 'q', 'WFAT', ''''//
cOLD     &           trim(adjustl(is))//s//' field frequency''', 'e', w, r)
cOLD            call wfat_dval_atleast(ce(i), 0.0d0, 'q', 'WFAT', ''''//
cOLD     &           trim(adjustl(is))//s//' field center''', 'e', w, r)
cOLD            call wfat_dval_atleast(wi(i), 0.0d0, 'g', 'WFAT', ''''//
cOLD     &           trim(adjustl(is))//s//' field width''', 'e', w, r)
cOLD         endif
cOLD         if (tp(i) /= 'external')
cOLD     &        call wfat_dval_atleast(mx(i), 0.0d0, 'q', 'WFAT', ''''//
cOLD     &        trim(adjustl(is))//s//' field max''', 'e', w, r)
cOLD      enddo
cOLD
cOLD      end
cOLDc=======================================================


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


      
cOLD      if (field_need) then
cOLD         call wfat_acdc_field(nfield, fld_tp, td)
cOLD         if (.not. rtdb_cput(rtdb, 'wfat:tdfield', 1, td))
cOLD     &        call errquit('WFAT: Failed writing field:tdfield '
cOLD     &        //'to rtdb', 1, INPUT_ERR)
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:nfield', mt_int, 1,
cOLD     &        nfield)) call errquit('WFAT: Failed writing field:nfield '
cOLD     &        //'to rtdb', 1, INPUT_ERR)
cOLD         if (.not. rtdb_cput(rtdb, 'wfat:field:type', nfield, fld_tp))
cOLD     &        call errquit('WFAT: Failed writing field:type '//
cOLD     &        'to rtdb', 1, INPUT_ERR)
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:freq', mt_dbl, nfield,
cOLD     &        fld_fr)) call errquit('WFAT: Failed writing field:freq '//
cOLD     &        'to rtdb', 1, INPUT_ERR)
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:center', mt_dbl, nfield,
cOLD     &        fld_ce)) call errquit('WFAT: Failed writing field:center '
cOLD     &        //'to rtdb', 1, INPUT_ERR)
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:width', mt_dbl, nfield,
cOLD     &        fld_wi)) call errquit('WFAT: Failed writing field:width '
cOLD     &        //'to rtdb', 1, INPUT_ERR)
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:max', mt_dbl, nfield,
cOLD     &        fld_mx)) call errquit('WFAT: Failed writing field:max '//
cOLD     &        'to rtdb', 1, INPUT_ERR)
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:cep', mt_dbl, nfield,
cOLD     &        fld_cp)) call errquit('WFAT: Failed writing field:cep '//
cOLD     &        'to rtdb', 1, INPUT_ERR)
cOLD
cOLD         if (.not. rtdb_put(rtdb, 'wfat:is_extfield', mt_log, 1,
cOLD     &        is_extfield)) call errquit('WFAT: Failed writing '//
cOLD     &        'is_extfield to rtdb', 1, INPUT_ERR)
cOLD         if (is_extfield) then
cOLD            if (.not. rtdb_cput(rtdb, 'wfat:ext_fname', 1, 
cOLD     &           trim(adjustl(extfield_name)) )) 
cOLD     &           call errquit('WFAT: Failed writing ext_fname to rtdb',
cOLD     &           1, INPUT_ERR)
cOLD            er = (/ ext_start, ext_end /)
cOLD            if (.not. rtdb_put(rtdb, 'wfat:ext_range', mt_int, 2,
cOLD     &           er)) call errquit('WFAT: Failed writing ext_range '//
cOLD     &           'to rtdb', 1, INPUT_ERR)
cOLD         endif
cOLD      endif

cOLD      if (dt_need) then
cOLD         if (is_extfield) then
cOLD            dt_00 = dt_ext
cOLD            nt_00 = ntime_ext
cOLD         else
cOLD            dt_00 = wfat_dt
cOLD            nt_00 = int(wfat_tmax/wfat_dt) + 1
cOLD         endif
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:dt', mt_dbl, 1, dt_00))
cOLD     &        call errquit('WFAT: Failed writing dt to rtdb',
cOLD     &        1, INPUT_ERR)
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:ntime', mt_int, 1, nt_00))
cOLD     &        call errquit('WFAT: Failed writing ntime to rtdb',
cOLD     &        1, INPUT_ERR)
cOLD      endif
cOLD
cOLD      if (tmax_need) then
cOLD         if (is_extfield) then
cOLD            tmax_00 = tmax_ext
cOLD         else
cOLD            tmax_00 = wfat_tmax
cOLD         endif
cOLD         if (.not. rtdb_put(rtdb, 'wfat:field:tmax', mt_dbl, 1,
cOLD     &        tmax_00))
cOLD     &        call errquit('WFAT: Failed writing tmax to rtdb',
cOLD     &        1, INPUT_ERR)
cOLD      endif
            

      if (nktmd_need) then
         if (.not. rtdb_put(rtdb, 'wfat:nktmd', mt_int, 1, nktmd))
     &        call errquit('WFAT: Failed writing nktmd to rtdb',
     &        1, INPUT_ERR)
      endif

cOLD      if (beta_tmd_need) then
cOLD         if (.not. rtdb_put(rtdb, 'wfat:beta_tmd', mt_dbl, 1, beta_tmd))
cOLD     &        call errquit('WFAT: Failed writing beta_tmd to rtdb',
cOLD     &        1, INPUT_ERR)
cOLD      endif
cOLD
cOLD      if (gamma_tmd_need) then
cOLD         if (.not. rtdb_put(rtdb, 'wfat:gamma_tmd', mt_dbl, 1,
cOLD     &        gamma_tmd))
cOLD     &        call errquit('WFAT: Failed writing gamma_tmd to rtdb',
cOLD     &        1, INPUT_ERR)
cOLD      endif

cOLD      if (lmax_need) then
cOLD         if (.not. rtdb_put(rtdb, 'wfat:lmax', mt_int, 1, wfat_lmax))
cOLD     &        call errquit('WFAT: Failed writing lmax to rtdb',
cOLD     &        1, INPUT_ERR)
cOLD      endif

cOLD      if (gqual_need) then
cOLD         if (.not. rtdb_put(rtdb, 'wfat:gqual', mt_int, 1, wfat_gqual))
cOLD     &        call errquit('WFAT: Failed writing gqual to rtdb',
cOLD     &        1, INPUT_ERR)
cOLD      endif
      
      if (omgpos_need) then
         if (.not. rtdb_put(rtdb, 'wfat:autopos', mt_log, 1, autopos))
     &        call errquit('WFAT: Failed writing autopos to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:omgpos', mt_dbl, 3, omgpos))
     &        call errquit('WFAT: Failed writing omgpos to rtdb',
     &        1, INPUT_ERR)
      endif

cOLD      if (intgl_need) then
cOLD         if (.not. rtdb_cput(rtdb, 'wfat:intgl_name', 1,
cOLD     &        trim(adjustl(wfat_intgl_name))))
cOLD     &        call errquit('WFAT: Failed writing intgl_name to RTDB.',
cOLD     &        1, INPUT_ERR)
cOLD      endif

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

      if (merg_need) then
         if (.not. rtdb_put(rtdb, 'wfat:use_manerg', mt_log, n_orbid,
     &        wfat_merg_l(1:n_orbid)))
     &        call errquit('WFAT: Failed writing use_manerg to RTDB',
     &        1, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:manerg', mt_dbl, n_orbid,
     &        wfat_merg(1:n_orbid)))
     &        call errquit('WFAT: Failed writing manerg to RTDB',
     &        1, RTDB_ERR)
      endif

cOLD      if (nlgl_need) then
cOLD         if (.not. rtdb_put(rtdb, 'wfat:nlgl', mt_int, 1, wfat_nlgl))
cOLD     &        call errquit('WFAT: Failed writing nlgl to RTDB', 1,
cOLD     &        RTDB_ERR)
cOLD      endif

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

      if (ffc_need) then
         if (.not. rtdb_put(rtdb, 'wfat:ffc_type', mt_int, n_orbid,
     &        ffc_type0))
     &        call errquit('WFAT: Failed writing ffc_type to RTDB', 1,
     &        RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:erg_ffc', mt_dbl, n_orbid,
     &        ergffc))
     &        call errquit('WFAT: Failed writing erg_ffc to RTDB', 1,
     &        RTDB_ERR)
      endif

cOLD      if (.true.) then
cOLD         if (.not. rtdb_put(rtdb, 'wfat:ozn', mt_int, 1, DEF_OZN))
cOLD     &        call errquit('WFAT: Failed writing ozn to RTDB', 1,
cOLD     &        RTDB_ERR)
cOLD      endif
cOLD
cOLD      if (.true.) then
cOLD         if (.not. rtdb_put(rtdb, 'wfat:cgrid_thr', mt_dbl, 1,
cOLD     &        DEF_CGRID_THR))
cOLD     &        call errquit('WFAT: Failed writing cgrid_thr to RTDB', 1,
cOLD     &        RTDB_ERR)
cOLD      endif
      
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


