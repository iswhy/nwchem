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

      character*2 tdfield
      logical calc_td, calc_dc, pp, gg, os, oc, or

      
      call util_get_rtdb_state(rtdb, os, oc, or)


c     ==========================================================
c     ==================== EXISTENCE CHECKS ====================
c     ==========================================================

c     ==== Checking existence of MANDATORY INPUTS ====
      calc_dc = (print_rate .or. print_tmd .or. print_tmd_atom)
      calc_td = (print_yild .or. print_avbeta .or. print_avgamma)
      if (print_exist) then
         if (calc_dc .and. calc_td) then
            write(luout, 302)
            call errquit('WFAT: Conflict in the ''print'' keywords',
     &           1, INPUT_ERR)
         endif
      else
         call errquit('WFAT: Cannot find the mandatory directive '//
     &        '''print''.', 0, INPUT_ERR)
      endif

c     ==== Checking existence of CONDITIONALLY MANDATORY INPUTS ====
      call wfat_acdc_field(nfield, fld_tp, tdfield)
      pp = (print_rate .or. print_yild .or. print_avbeta .or. 
     &     print_avgamma .or. print_tmd .or. print_tmd_atom)
      if (.not. field_exist .and. pp) then
         write(luout, 205) 'field', '* print_rate', '* print_yild', 
     &        '* print_avbeta', '* print_avgamma', '* print_tmd',
     &        '* print_tmd_atom', ' '
         call errquit('INPUT ERROR: Cannot find input directive '//
     &        '''field''.', 0, INPUT_ERR)
      else if (field_exist) then
         if (calc_dc .and. tdfield=='td') then
            call errquit('WFAT: When either of ''print_rate'', '//
     &           '''print_tmd'', or ''print_tmd_atom'' is invoked, '//
     &           'the field type must ONLY be ''dc''', 1, INPUT_ERR)
         else if (calc_td .and. tdfield=='ti') then
            call errquit('WFAT: When either of ''print_yild'', '//
     &           '''print_avbeta'', or ''print_avgamma'' is invoked, '//
     &           'no field with type ''dc'' is allowed', 1, INPUT_ERR)
         else if (.not. calc_dc .and. .not. calc_td) then
            write(luout, 201) 'field',
     &           '* print_rate      ',
     &           '* print_yild      ',
     &           '* print_avbeta    ',
     &           '* print_avgamma   ',
     &           '* print_tmd       ',
     &           '* print_tmd_atom  ', ' '
         endif
      endif

c     == dt ==
      if (.not. dt_exist .and. tdfield=='td') then
         call errquit('WFAT: Input directive ''dt'' is required when '//
     &        'non-DC fields are used', 1, INPUT_ERR)
      else if (dt_exist .and. (tdfield=='ti' .or. tdfield=='no')) then
         write(luout, 206) 'dt'
      endif

c     == tmax ==
      if (.not. tmax_exist .and. tdfield=='td') then
         call errquit('WFAT: Input directive ''tmax'' is required '//
     &        'when non-DC fields are used', 1, INPUT_ERR)
      else if (tmax_exist .and. (tdfield=='ti' .or. tdfield=='no')) then
         write(luout, 206) 'tmax' 
      endif
      
c     == beta ==
      gg = (print_srfc .or. print_rate .or. print_yild .or.
     &      print_avbeta .or. print_avgamma)
      if (.not. beta_exist .and. gg) then
         write(luout, 205) 'beta', '* print_srfc', '* print_rate', 
     &        '* print_yild', '* print_avbeta', '* print_avgamma', ' ',
     &        ' '
         call errquit('INPUT ERROR: Cannot find input directive '//
     &        '''beta''. ', 0, INPUT_ERR)
      else if (beta_exist .and. .not. gg) then
         write(luout, 201) 'beta',
     &        '* print_rate      ',
     &        '* print_yild      ',
     &        '* print_avbeta    ',
     &        '* print_avgamma   ',         
     &        '* print_srfc      ', ' ', ' '
      endif
c     == gamma ==
      if (.not. gamma_exist .and. gg) then
         write(luout, 205) 'gamma', '* print_srfc', '* print_rate',
     &        '* print_yild', '* print_avbeta', '* print_avgamma', ' ',
     &        ' '
         call errquit('INPUT ERROR: Cannot find input directive '//
     &        '''gamma''. ', 0, INPUT_ERR)
      else if (gamma_exist .and. .not. gg) then
         write(luout, 201) 'gamma',
     &        '* print_rate      ',
     &        '* print_yild      ',
     &        '* print_avbeta    ',
     &        '* print_avgamma   ',
     &        '* print_srfc      ', ' ', ' '
      endif
c     == beta_tmd ==
      if (.not. beta_tmd_exist .and.
     &     (print_tmd .or. print_tmd_atom)) then
         write(luout, 205) 'beta_tmd', '* print_tmd',
     &        '* print_tmd_atom', ' ', ' ', ' ', ' ', ' '
         call errquit('INPUT ERROR: Cannot find input directive '//
     &        '''beta_tmd''. ', 0, INPUT_ERR)
      else if (beta_tmd_exist .and. .not.
     &        (print_tmd .or. print_tmd_atom)) then
         write(luout, 201) 'beta_tmd',
     &        '* print_tmd       ',
     &        '* print_tmd_atom  ', ' ', ' ', ' ', ' ', ' '
      endif
c     == gamma_tmd ==
      if (.not. gamma_tmd_exist .and.
     &     (print_tmd .or. print_tmd_atom)) then
         write(luout, 205) 'gamma_tmd', '* print_tmd',
     &        '* print_tmd_atom', ' ', ' ', ' ', ' ', ' '
         call errquit('INPUT ERROR: Cannot find input directive '//
     &        '''gamma_tmd''. ', 0, INPUT_ERR)
      else if (gamma_tmd_exist .and. .not.
     &        (print_tmd .or. print_tmd_atom)) then
         write(luout, 201) 'gamma_tmd',
     &        '* print_tmd       ',
     &        '* print_tmd_atom  ', ' ', ' ', ' ', ' ', ' '
      endif
      
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
      if (.not. orbid_exist) then
         if (wfat_level == 'oe') then
            n_orbid = 2
            wfat_spinid(1:n_orbid) = (/ 1,  1 /)
            wfat_orbid0(1:n_orbid) = (/ 0, -1 /)
         else if (wfat_level == 'me') then
            n_orbid = 1
         endif
         orbid_stat = .true.
      else if (orbid_exist .and. or) then
         write(luout, 202) 'movec_id', 'movec_id'
      endif
      if (.not. paraqn_exist) then
         n_npara = 3
         npara(1:n_npara) = 0
         mpara(1:n_npara) = (/ 0, -1, +1 /)
         paraqn_stat = .true.
      else if (paraqn_exist .and. or) then
         write(luout, 202) 'paraqn', 'paraqn'
      endif
      if (.not. lmax_exist) then
         wfat_lmax = def_wfat_lmax
         lmax_stat = .true.
      else if (lmax_exist .and. or) then
         write(luout, 202) 'lmax', 'lmax'
      endif
cOLD      if (.not. gqual_exist) then
cOLD         wfat_gqual = def_gqual
cOLD         gqual_stat = .true.
cOLD      else if (gqual_exist .and. or) then
cOLD         write(luout, 202) 'gqual', 'gqual'
cOLD      endif
      if (.not. omgpos_exist) then
         autopos = .true.
         omgpos = 0.0d0
         omgpos_stat = .true.
      else if (omgpos_exist .and. or) then
         write(luout, 202) 'omgpos', 'omgpos'
      endif
      if (.not. dgnthr_exist) then
         dgnthr = DEF_DGNTHR
         dgnthr_stat = .true.
      else if (dgnthr_exist .and. or) then
         write(luout, 202) 'degenthr', 'degenthr'
      endif      
      if (.not. intgl_exist) then
         call wfat_def_intglname(rtdb, wfat_level, wfat_intgl_name)
         intgl_stat = .true.
      endif
      if (.not. wfat1_exist) then
         wfat1 = DEF_WFAT1
         wfat1_stat = .true.
      else if (wfat1_exist .and. or) then
         write(luout, 202) 'wfat1', 'wfat1'
      endif
      if (.not. w1x_exist) then
         wfat1_exclude = DEF_W1X
         w1x_stat = .true.
      else if (w1x_exist) then
         if (or) then
            write(luout, 202) 'wfat1_exclude', 'wfat1_exclude'
         else
            call wfat_extract_w1x(CHID_CHARLEN,
     &           MAX_N_NPARA*MAX_N_ORBID, wfat1_exclude, n_orbid,
     &           wfat_orbid0, wfat_spinid, n_npara, npara, mpara,
     &           wfat1_exclude_l)
            w1x_stat = .true.
         endif
      endif
      if (.not. mdip_exist) then
         wfat_mdip(:,:) = DEF_MDIP
         mdip_stat = .true.
      else
         if (or) then
            write(luout, 202) 'mydipole', 'mydipole'
         else
            call wfat_process_mydip(ORBID_CHARLEN, wfat_mdip_s, n_orbid,
     &           wfat_orbid0, wfat_spinid, wfat_mdip, wfat_mdip_l)
            mdip_stat = .true.
         endif
      endif
      if (.not. merg_exist) then
         wfat_merg(:) = DEF_MERG
         merg_stat = .true.
      else
         if (or) then
            write(luout, 202) 'myenergy', 'myenergy'
         else
            call wfat_process_myerg(ORBID_CHARLEN, wfat_merg_s, n_orbid,
     &           wfat_orbid0, wfat_spinid, wfat_merg, wfat_merg_l)
            merg_stat = .true.
         endif
      endif
      if (.not. nlgl_exist) then
         wfat_nlgl = DEF_NLGL
         nlgl_stat = .true.
      else if (nlgl_exist .and. or) then
         write(luout, 202) 'nlgl', 'nlgl'
      endif
      if (.not. dthr_exist) then
         wfat_dthr = DEF_WFAT_DTHR
         dthr_stat = .true.
      else if (dthr_exist .and. or) then
         write(luout, 202) 'diagthr', 'diagthr'
      endif
      if (.not. wfat_mask_exist) then
         wfat_mask = DEF_WFAT_MASK
         wfat_mask_param = (/ -1.0d6, -1.0d6 /)
         wfat_mask_stat = .true.
      else if (wfat_mask_exist .and. or) then
         write(luout, 202) 'tailmask', 'tailmask'
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
         if (.not. sitype_exist) then
            si_type = def_si_type
            sitype_stat = .true.
         else if (sitype_exist .and. or) then
            write(luout, 202) 'si_type', 'si_type'
         endif
         if (wfat1) call errquit('MEWFAT: At the moment, MEWFAT '//
     &        'cannot be used with the first order correction. '//
     &        'Please remove the ''wfat1'' keyword from your input '//
     &        'file.', 62, UERR)
      endif
      if (.not. prtb_type_exist) then
         prtb_type = DEF_PRTB_type
         prtb_type_stat = .true.
cOLD      else if (prtb_type_exist .and. or) then
cOLD         write(luout, 203)
      endif
      if (.not. prtb_ord_exist) then
         if (prtb_type=='bw') prtb_ord = DEF_PRTB_ORD_BW
         if (prtb_type=='rs') prtb_ord = DEF_PRTB_ORD_RS
         prtb_ord_stat = .true.
cOLD      else if (prtb_ord_exist .and. or) then
cOLD         write(luout, 204) 'order', 'order'
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
      
c     ==== Checking format of MANDATORY INPUTS ====
      if (.not. print_stat) then
         call errquit('WFAT: Error reading mandatory input parameter '//
     &        '''print''. Please check if it is specified following '//
     &        'the format ''print [srfc] [rate] [yield] [avbeta] '//
     &        '[avgamma] [tmd] [tmd_atom]''. At least one of the '//
     &        'seven options must be present.', 0, INPUT_ERR)
      endif
      
c     ==== Checking format of CONDITIONALLY MANDATORY INPUTS ====
      pp = (print_rate .or. print_yild .or. print_avbeta .or. 
     &     print_avgamma .or. print_tmd .or. print_tmd_atom)
      if (.not. field_stat .and. pp) then
         call errquit('WFAT: The format of input directive '//
     &        '''field'' is incorrect', 0, INPUT_ERR)
      endif

c     == dt ==
      if (.not. dt_stat .and. tdfield=='td')
     &     call errquit('WFAT: The format of input directive '//
     &     '''dt'' is incorrect', 0, INPUT_ERR)
c     == tmax ==
      if (.not. tmax_stat .and. tdfield=='td')
     &     call errquit('WFAT: The format of input directive '//
     &     '''tmax'' is incorrect', 0, INPUT_ERR)      
c     == beta ==
      gg = (print_srfc .or. print_rate .or. print_yild .or.
     &     print_avbeta .or. print_avgamma)
      if (.not. beta_stat .and. gg) then
         call errquit('WFAT: The format of input directive '//
     &        '''beta'' is incorrect. The correct format '//
     &        'should be '//
     &        '''beta <double start_angle> [<integer npoints> '//
     &        '<double last_angle>]''. If only ''start_angle'' is '//
     &        'given or when ''npoints=1'', ''beta'' is assumed to '//
     &        'be equal to ''start_angle''.', 0, INPUT_ERR)
      endif
c     == gamma ==
      if (.not. gamma_stat .and. gg) then
         call errquit('WFAT: The format of input directive '//
     &        '''gamma'' is incorrect. The correct format '//
     &        'should be '//
     &        '''gamma <double start_angle> [<integer npoints> '//
     &        '<double last_angle>]''. If only ''start_angle'' is '//
     &        'given or when ''npoints=1'', ''gamma'' is assumed to '//
     &        'be equal to ''start_angle''.', 0, INPUT_ERR)
      endif
c     == beta_tmd ==
      if (.not. beta_tmd_stat .and.
     &     (print_tmd .or. print_tmd_atom)) then
         call errquit('WFAT: The format of input directive '//
     &        '''beta_tmd'' is incorrect. The correct format '//
     &        'should be ''beta_tmd <float value>''.',
     &        0, INPUT_ERR)
      endif
c     == gamma_tmd ==
      if (.not. gamma_tmd_stat .and.
     &     (print_tmd .or. print_tmd_atom)) then
         call errquit('WFAT: The format of input directive '//
     &        '''gamma_tmd'' is incorrect. The correct format '//
     &        'should be ''gamma_tmd <float value>''.',
     &        0, INPUT_ERR)
      endif

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
      if (.not. paraqn_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''paraqn''. Please check if it is '//
     &        'specified following the format '//
     &        '''paraqn "n<integer i>" <integer m#1> [<integer m#2> '//
     &        '...] ["n<integer j>" <integer m#1> [<integer m#2> ... '//
     &        ']] .... For example, the following input '//
     &        '''paraqn "n0" 0 -1 1 "n1" 0'' corresponds to the '//
     &        'following list of (n,m) parabolic quantum numbers: '//
     &        '(0,0); (0,-1); (0,+1); (1,0).',
     &        0, INPUT_ERR)
      endif      
cOLD      if (.not. gqual_stat) then
cOLD         call errquit('WFAT: Error reading optional input parameter '//
cOLD     &        '''gridqual''. Please check if it is specified '//
cOLD     &        'following the format ''gridqual <integer value>''.',
cOLD     &        0, INPUT_ERR)
cOLD      endif
      if (.not. omgpos_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''omgpos''. Please check if it is specified '//
     &        'following the format ''omgpos <string ''auto''> || '//
     &        '<float valuex  float valuey  float valuez>''',
     &        0, INPUT_ERR)
      endif
      if (.not. dgnthr_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''degenthr''. Please check if it is specified '//
     &        'following the format ''degenthr <float value>''.',
     &        0, INPUT_ERR)
      endif
      if (.not. sitype_stat .and. wfat_level == 'me') then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''si_type''. Please check if it is specified '//
     &        'following the format ''si_type <value>'', where '//
     &        '<value> is any of ''none'', ''mo'', and ''dyson''.',
     &        0, INPUT_ERR)
      endif
      if (.not. lmax_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''lmax''. Please check if it is specified following '//
     &        'the format ''lmax <integer value>''.',
     &        0, INPUT_ERR)
      endif
      if (.not. intgl_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''intgl''. Please check if it is specified '//
     &        'following the format ''intgl "<string value>"'' '//
     &        '(do not skip the double quotes).',
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
      if (.not. w1x_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''wfat1_exclude''. Please check if it is specified '//
     &        'following the format ''wfat1_exclude '//
     &        '"spin, oid" ...'' where spin is ''a'' or ''b'', '//
     &        'oid is movec_id, and n and m are the parabolic '//
     &        'quantum number.', 0, INPUT_ERR)
      endif

      if (.not. mdip_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''mydipole''. Please check if it is specified '//
     &        'following the format ''mydipole '//
     &        '"spin, oid" dipx dipy dipz ...'' where spin is ''a'' '//
     &        'or ''b'', oid is movec_id, and dipx dipy dipz are the '//
     &        'three Cartesian components of dipole.', 0, INPUT_ERR)
      endif

      if (.not. merg_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''myenergy''. Please check if it is specified '//
     &        'following the format ''myenergy '//
     &        '"spin, oid" erg ...'' where spin is ''a'' '//
     &        'or ''b'', oid is movec_id, and erg is the '//
     &        'manually given ionization energy.', 0, INPUT_ERR)
      endif
      
      if (.not. nlgl_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''nlgl''. Please check if it is specified following '//
     &        'the format ''nlgl <integer value>''.',
     &        0, INPUT_ERR)
      endif
      if (.not. dthr_stat) then
         call errquit('WFAT: Error reading optional input parameter '//
     &        '''diagthr''. Please check if it is specified '//
     &        'following the format ''diagthr <float value>''.',
     &        0, INPUT_ERR)
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
         call wfat_need_input()
      else if (wfat_level == 'me') then 
cdebug         call mewfat_dummy()
         call mewfat_need_input()
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

 202  format(/,
     &     3x, '*** ATTENTION ***', /
     &     3x, 'Input directive ''', a, ''' is found ',
     &     'but will be ignored because the present calculation ', /
     &     3x, 'is a restarted job in which case the information ',
     &     'about ''', a, ''' is obtained from RTDB.')

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


c=======================================================
c=======================================================
      subroutine wfat_acdc_field(nf, tp, tdf)

      implicit none
#include "errquit.fh"
#include "wfat_input.fh"
      integer nf                         ! [input],
      character(len=FLD_TP_LEN) tp(nf)   ! [input],
      character*2 tdf                    ! [output],

      integer i

      tdf = 'no'
      do i = 1, nf
         if (tp(i) == 'dc') then
            tdf = 'ti'
         endif

         if (tdf=='ti' .and. (tp(i)=='gaussian' .or. tp(i)=='hann' .or.
     &        tp(i)=='external'))
     &        call errquit('WFAT: When there is a ''dc'' type '//
     &        'field in the ''field'' directive, no other field '//
     &        'types are allowed to coexist', 1, INPUT_ERR)
         
         if (tp(i)=='gaussian' .or. tp(i)=='hann' .or.
     &        tp(i)=='external') then
            tdf = 'td'
         endif
      enddo

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
      subroutine mewfat_need_input()

      implicit none
#include "wfat_input.fh"
#include "stdio.fh"
#include "errquit.fh"      
      logical sr, rt, yi, avb, avg, tmd, tmda
      character*1 m


      field_need = .false.          ! needed when OE or ME
      beta_need = .false.           ! needed when OE or ME
      gamma_need = .false.          ! needed when OE or ME
      beta_tmd_need = .false.       ! needed when OE or ME
      gamma_tmd_need = .false.      ! needed when OE or ME
      nktmd_need = .false.          ! needed when OE or ME
      orbid_need = .false.          ! needed when OE
      paraqn_need = .false.         ! needed when OE or ME
cOLD      gqual_need = .false.          ! needed when OE or ME
      omgpos_need = .false.         ! needed when OE or ME
      dgnthr_need = .false.         ! needed when oe or me
      lmax_need = .false.           ! needed when OE or ME
      intgl_need = .false.          ! needed when OE or ME
      movec_need = .false.          ! needed when OE or ME
      movec1_need = .false.         ! needed when ME
      sitype_need = .false.         ! needed when ME
      wfat1_need = .false.          ! needed when OE or ME
      w1x_need = .false.            ! needed when oe or me
      mdip_need = .false.           ! needed when oe or me
      merg_need = .false.           ! needed when oe or me
      nlgl_need = .false.           ! needed when oe or me
      dthr_need = .false.           ! needed when oe or me
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


      sr = print_srfc
      rt = print_rate
      yi = print_yild
      avb = print_avbeta
      avg = print_avgamma
      tmd = print_tmd
      tmda = print_tmd_atom

            
      call wfat_need_input0(rt, tmd, tmda, yi, avb, avg, .false.,
     &     field_need)
      call wfat_need_input0(yi, avb, avg, .false., .false., .false.,
     &     .false., dt_need)
      call wfat_need_input0(yi, avb, avg, .false., .false., .false.,
     &     .false., tmax_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     beta_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     gamma_need)
c TMD notyet      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
c TMD notyet     &     .false., .false., beta_tmd_need)
c TMD notyet      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
c TMD notyet     &     .false., .false., gamma_tmd_need)
c TMD notyet      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
c TMD notyet     &     .false., .false., nktmd_need)

      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     paraqn_need)
cOLD      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
cOLD     &     gqual_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     omgpos_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     dgnthr_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     lmax_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     intgl_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     sitype_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     w1x_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     mdip_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     merg_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     nlgl_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     dthr_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     ffc_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat_mask_need)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_need_input()

      implicit none
#include "wfat_input.fh"
#include "stdio.fh"
      logical sr, rt, yi, avb, avg, tmd, tmda
      character*1 m

      
      field_need = .false.              ! needed when oe or me
      beta_need = .false.               ! needed when oe or me
      gamma_need = .false.              ! needed when oe or me
      beta_tmd_need = .false.           ! needed when oe or me
      gamma_tmd_need = .false.          ! needed when oe or me
      nktmd_need = .false.              ! needed when oe or me
      orbid_need = .false.              ! needed when oe
      paraqn_need = .false.             ! needed when oe or me
cOLD      gqual_need = .false.              ! needed when oe or me
      omgpos_need = .false.             ! needed when oe or me
      dgnthr_need = .false.             ! needed when oe or me
      lmax_need = .false.               ! needed when oe or me
      intgl_need = .false.              ! needed when oe or me
      movec_need = .false.              ! needed when oe or me
      movec1_need = .false.             ! needed when me
      wfat1_need = .false.              ! needed when oe or me
      w1x_need = .false.                ! needed when oe or me
      mdip_need = .false.               ! needed when oe or me
      merg_need = .false.               ! needed when oe or me
      nlgl_need = .false.               ! needed when oe or me
      dthr_need = .false.               ! needed when oe or me
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

      
      sr = print_srfc
      rt = print_rate
      yi = print_yild
      avb = print_avbeta
      avg = print_avgamma
      tmd = print_tmd
      tmda = print_tmd_atom

      
      call wfat_need_input0(rt, tmd, tmda, yi, avb, avg, .false.,
     &     field_need)
      call wfat_need_input0(yi, avb, avg, .false., .false., .false.,
     &     .false., dt_need)
      call wfat_need_input0(yi, avb, avg, .false., .false., .false.,
     &     .false., tmax_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     beta_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     gamma_need)
      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
     &     .false., .false., beta_tmd_need)
      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
     &     .false., .false., gamma_tmd_need)
      call wfat_need_input0(tmd, tmda, .false., .false., .false., 
     &     .false., .false., nktmd_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     orbid_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     paraqn_need)
cOLD      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
cOLD     &     gqual_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     omgpos_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     dgnthr_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     lmax_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     intgl_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     movec_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     wfat1_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     w1x_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     mdip_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     merg_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     nlgl_need)
      call wfat_need_input0(sr, rt, tmd, tmda, yi, avb, avg,
     &     dthr_need)
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
      character*1 spin_str
      character*10 ss1, ss2
      character*3 n_npara_c, MAX_N_NPARA_c
      character*5 n_orbid_c, MAX_N_ORBID_c

      integer nsp, numelc(2), i, is, nrr, ntt, ffcount
      double precision dt_00, tmax_00
      
      character*1 w
      logical r, unit_ffc
      parameter ( w='i' )
      integer k
      

      if (field_need) then
         call wfat_sanity_field(nfield, fld_tp, fld_fr, fld_ce, fld_wi,
     &        fld_mx)
      endif

      if (dt_need) then
         if (is_extfield) then
            dt_00 = dt_ext
         else
            dt_00 = wfat_dt
         endif
         call wfat_dval_atleast(dt_00, 0.0d0, 'g', 'WFAT',
     &        '''dt''', 'e', w, r)
      endif

      if (tmax_need) then
         if (is_extfield) then
            tmax_00 = tmax_ext
         else
            tmax_00 = wfat_tmax
         endif
         call wfat_dval_atleast(tmax_00, 0.0d0, 'g', 'WFAT',
     &        '''tmax''', 'e', w, r)
      endif
      
      if (beta_need) then
         call wfat_dval_inrange(wfat_beta_start, 0.0d0, 'q', 180.0d0, 
     &        'q', 'WFAT', '''start_angle'' of ''beta''', 'f', w, r)
         call wfat_dval_inrange(wfat_beta_last, 0.0d0, 'q', 180.0d0,
     &        'q', 'WFAT', '''last_angle'' of ''beta''', 'f', w, r)
         call wfat_ival_atleast(wfat_nbeta, 1, 'WFAT',
     &        '''npoints'' of ''beta''', w, r)
      endif

      if (gamma_need) then
         call wfat_dval_inrange(wfat_gamma_start, 0.0d0, 'q', 360.0d0, 
     &       'q', 'WFAT', '''start_angle'' of ''gamma''', 'f', w, r)
         call wfat_dval_inrange(wfat_gamma_last, 0.0d0, 'q', 360.0d0, 
     &        'q', 'WFAT', '''last_angle'' of ''gamma''', 'f', w, r)
         call wfat_ival_atleast(wfat_ngamma, 1, 'WFAT',
     &        '''npoints'' of ''gamma''', w, r)
      endif
      
      if (orbid_need) then
         write(n_orbid_c, '(i5)') n_orbid
         write(MAX_N_ORBID_c, '(i5)') MAX_N_ORBID
         if (n_orbid >= MAX_N_ORBID) then
            call errquit('WFAT: The number of molecular orbitals '//
     &           'specified in ''movec_id'' (='//
     &           trim(adjustl(n_orbid_c))//
     &           ') must be less than or equal to the maximum value, '//
     &           'which is '//trim(adjustl(MAX_N_ORBID_c))//'.',
     &           1, INPUT_ERR)
         endif

c        ==== Check spin ID ====
         if (s2a1e) then
            nsp = 2
         else
            if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, nsp))
     &           call errquit('WFAT: Error in getting ipol from RTDB',
     &           1, RTDB_ERR)
         endif
         call wfat_ivec_inrange(n_orbid, wfat_spinid, 1, nsp, 'WFAT', 
     &        '''spin_id''', 'o', r, k)
         if (.not. r)
     &        call errquit('WFAT: One spin channel specified has '//
     &        'an illegal value. The allowed value is ''a'' when the '//
     &        'system is spin-unpolarized. If it is spin-polarized, ' //
     &        'the possible spin is either ''a'' or ''b''.',
     &        1, INPUT_ERR)

c        ==== Check orbitals ID ====
         if (s2a1e) then
            numelc = (/ 1, 0 /)
         else
            if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, numelc))
     &           call errquit('WFAT: Error in getting noc from RTDB',
     &           1, RTDB_ERR)
         endif
         do i = 1, n_orbid
            is = wfat_spinid(i)
            if (is == 1) spin_str = 'a'
            if (is == 2) spin_str = 'b'
            if (s2a1e) then
c               if (.not. rtdb_get(rtdb, 'wfat:2a1e_nr', mt_int, 1, nrr))
c     &              call errquit('WFAT: Error in getting 2a1e_nr '//
c     &              'from RTDB', 1, RTDB_ERR)
c               if (.not. rtdb_get(rtdb, 'wfat:2a1e_nt', mt_int, 1, ntt))
c     &              call errquit('WFAT: Error in getting 2a1e_nt '//
c     &              'from RTDB', 1, RTDB_ERR)

               ntt = wfat2a_nt(1)
               if (wfat2a_fedvr) then
                  call wfat_get_fedvrnr(wfat2a_nfe, wfat2a_ndvr, nrr)
               else
                  nrr = wfat2a_nr(1)
               endif
               call wfat_ival_inrange(wfat_orbid0(i), 1, nrr*ntt, 
     &              'WFAT', '''orb_id''', 'o', r)
               write(ss1, '(i10)') 1
               write(ss2, '(i10)') nrr*ntt
            else
               call wfat_ival_inrange(wfat_orbid0(i), -(numelc(is)-1), 
     &              0, 'WFAT', '''orb_id''', 'o', r)
               write(ss1, '(i10)') -(numelc(is)-1)
               write(ss2, '(i10)') 0
            endif
            if (.not. r) then
               call errquit('WFAT: One spatial orbital in '//
     &              spin_str//'-spin channel lies outside of the '//  
     &              'allowed range, which is in between and '//
     &              'including '//trim(adjustl(ss1))//' and '//
     &              trim(adjustl(ss2))//'.', 1, INPUT_ERR)
            endif
         enddo
      endif

c     ==== Check parabolic quantum numbers ====
      if (paraqn_need) then
         write(n_npara_c, '(i3)') n_npara
         write(MAX_N_NPARA_c, '(i3)') MAX_N_NPARA
         if (n_npara >= MAX_N_NPARA) then
            call errquit('WFAT: The number of integers specified for'//
     &           '''npara'' (='//trim(adjustl(n_npara_c))//') '//
     &           'must be less than or equal to the maximum value, '//
     &           'which is '//trim(adjustl(MAX_N_NPARA_c))//'.',
     &           1, INPUT_ERR)
         endif
         call wfat_ivec_atleast(n_npara, npara, 0, 'WFAT', '''npara''',
     &        w, r, k)
      endif

      if (nktmd_need) then
         call wfat_ival_atleast(nktmd, 1, 'WFAT', '''nktmd''', w, r)
      endif

      if (beta_tmd_need) then
         call wfat_dval_inrange(beta_tmd, 0.0d0, 'q', 180.0d0, 'q',
     &        'WFAT', '''beta_tmd''', 'f', w, r)
      endif

      if (gamma_tmd_need) then
         call wfat_dval_inrange(gamma_tmd, 0.0d0, 'q', 360.0d0, 'q',
     &        'WFAT', '''gamma_tmd''', 'f', w, r)
      endif

cOLD      if (gqual_need) then
cOLD         call wfat_ival_inrange(wfat_gqual, 1, 4, 'WFAT', 
cOLD     &        '''gridqual''', w, r)
cOLD      endif

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
      
      if (dgnthr_need) then
         call wfat_dval_atleast(dgnthr, 0.0d0, 'g', 'WFAT', 
     &        '''degenthr''', 'e', w, r)
      endif

      if (sitype_need) then
         if (si_type/='none' .and. si_type/='mo' .and. si_type/='dyson')
     &        call errquit('WFAT: The value of si_type ('//si_type// 
     &        ') does not match the available options, which are '//
     &        '''none'', ''mo'', or ''dyson''.', 31, INPUT_ERR)
      endif

      if (nlgl_need) then
         call wfat_ival_atleast(wfat_nlgl, 0, 'WFAT', '''nlgl''', w, r)
      endif

      if (dthr_need) then
         call wfat_dval_atleast(wfat_dthr, 0.0d0, 'g', 'WFAT', 
     &        '''diagthr''', 'e', w, r)
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


c=======================================================
c=======================================================
      subroutine wfat_sanity_field(nf, tp, fr, ce, wi, mx)

      implicit none
#include "wfat_input.fh"
      integer nf
      character(len=FLD_TP_LEN) tp(nf)
      double precision fr(max_nfield)
      double precision ce(max_nfield)
      double precision wi(max_nfield)
      double precision mx(max_nfield)

      character*2 s
      character*3 is
      character*1 w
      logical r
      integer i
      parameter ( w='i' )
      
      
      do i = 1, nf
         write(is, '(i3)') i
         call wfat_ordinal_suffix(i, s)
         if (tp(i) /= 'dc' .and. tp(i) /= 'external') then
            call wfat_dval_atleast(fr(i), 0.0d0, 'q', 'WFAT', ''''//
     &           trim(adjustl(is))//s//' field frequency''', 'e', w, r)
            call wfat_dval_atleast(ce(i), 0.0d0, 'q', 'WFAT', ''''//
     &           trim(adjustl(is))//s//' field center''', 'e', w, r)
            call wfat_dval_atleast(wi(i), 0.0d0, 'g', 'WFAT', ''''//
     &           trim(adjustl(is))//s//' field width''', 'e', w, r)
         endif
         if (tp(i) /= 'external')
     &        call wfat_dval_atleast(mx(i), 0.0d0, 'q', 'WFAT', ''''//
     &        trim(adjustl(is))//s//' field max''', 'e', w, r)
      enddo

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
      logical os, oc, or, nosave
      integer nt_00, er(2), nw1x_
      double precision dt_00, tmax_00


      call util_get_rtdb_state(rtdb, os, oc, or)
c     .and. os means write that input into RTDB only when it is a start-type job.
c     An example of input that will only be written into the RTDB when doing a
c     start-type job (os=.true.) is lmax. If it is written into the RTDB when
c     the job is not of start-type, meaning restart-type, it will overwrite lmax
c     that is already in the RTDB.

      

c forcheck      if (.not. rtdb_print(rtdb, .false.))
c forcheck     &     call errquit('WFAT: Printing rtdb failed', 1, RTDB_ERR)


      if (.not. rtdb_put(rtdb, 'wfat:print_srfc', mt_log, 1, 
     &     print_srfc))
     &     call errquit('WFAT: Failed writing print_srfc to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:print_rate', mt_log, 1, 
     &     print_rate))
     &     call errquit('WFAT: Failed writing print_rate to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:print_yild', mt_log, 1, 
     &     print_yild))
     &     call errquit('WFAT: Failed writing print_yild to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:print_avbeta', mt_log, 1, 
     &     print_avbeta))
     &     call errquit('WFAT: Failed writing print_avbeta to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:print_avgamma', mt_log, 1, 
     &     print_avgamma))
     &     call errquit('WFAT: Failed writing print_avgamma to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:print_tmd', mt_log, 1, 
     &     print_tmd))
     &     call errquit('WFAT: Failed writing print_tmd to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:print_tmd_atom', mt_log, 1, 
     &     print_tmd_atom))
     &     call errquit('WFAT: Failed writing print_tmd_atom to rtdb',
     &     1, INPUT_ERR)

      
      if (field_need) then
         call wfat_acdc_field(nfield, fld_tp, td)
         if (.not. rtdb_cput(rtdb, 'wfat:tdfield', 1, td))
     &        call errquit('WFAT: Failed writing field:tdfield '
     &        //'to rtdb', 1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:field:nfield', mt_int, 1,
     &        nfield)) call errquit('WFAT: Failed writing field:nfield '
     &        //'to rtdb', 1, INPUT_ERR)
         if (.not. rtdb_cput(rtdb, 'wfat:field:type', nfield, fld_tp))
     &        call errquit('WFAT: Failed writing field:type '//
     &        'to rtdb', 1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:field:freq', mt_dbl, nfield,
     &        fld_fr)) call errquit('WFAT: Failed writing field:freq '//
     &        'to rtdb', 1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:field:center', mt_dbl, nfield,
     &        fld_ce)) call errquit('WFAT: Failed writing field:center '
     &        //'to rtdb', 1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:field:width', mt_dbl, nfield,
     &        fld_wi)) call errquit('WFAT: Failed writing field:width '
     &        //'to rtdb', 1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:field:max', mt_dbl, nfield,
     &        fld_mx)) call errquit('WFAT: Failed writing field:max '//
     &        'to rtdb', 1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:field:cep', mt_dbl, nfield,
     &        fld_cp)) call errquit('WFAT: Failed writing field:cep '//
     &        'to rtdb', 1, INPUT_ERR)

         if (.not. rtdb_put(rtdb, 'wfat:is_extfield', mt_log, 1,
     &        is_extfield)) call errquit('WFAT: Failed writing '//
     &        'is_extfield to rtdb', 1, INPUT_ERR)
         if (is_extfield) then
            if (.not. rtdb_cput(rtdb, 'wfat:ext_fname', 1, 
     &           trim(adjustl(extfield_name)) )) 
     &           call errquit('WFAT: Failed writing ext_fname to rtdb',
     &           1, INPUT_ERR)
            er = (/ ext_start, ext_end /)
            if (.not. rtdb_put(rtdb, 'wfat:ext_range', mt_int, 2,
     &           er)) call errquit('WFAT: Failed writing ext_range '//
     &           'to rtdb', 1, INPUT_ERR)
         endif
      endif

      if (dt_need) then
         if (is_extfield) then
            dt_00 = dt_ext
            nt_00 = ntime_ext
         else
            dt_00 = wfat_dt
            nt_00 = int(wfat_tmax/wfat_dt) + 1
         endif
         if (.not. rtdb_put(rtdb, 'wfat:field:dt', mt_dbl, 1, dt_00))
     &        call errquit('WFAT: Failed writing dt to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:field:ntime', mt_int, 1, nt_00))
     &        call errquit('WFAT: Failed writing ntime to rtdb',
     &        1, INPUT_ERR)
      endif

      if (tmax_need) then
         if (is_extfield) then
            tmax_00 = tmax_ext
         else
            tmax_00 = wfat_tmax
         endif
         if (.not. rtdb_put(rtdb, 'wfat:field:tmax', mt_dbl, 1,
     &        tmax_00))
     &        call errquit('WFAT: Failed writing tmax to rtdb',
     &        1, INPUT_ERR)
      endif
      
      if (beta_need) then
         if (.not. rtdb_put(rtdb, 'wfat:beta_start', mt_dbl, 1,
     &        wfat_beta_start))
     &        call errquit('WFAT: Failed writing beta_start to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:beta_last', mt_dbl, 1,
     &        wfat_beta_last))
     &        call errquit('WFAT: Failed writing beta_last to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:nbeta', mt_int, 1, wfat_nbeta))
     &        call errquit('WFAT: Failed writing nbeta to rtdb',
     &        1, INPUT_ERR)
      endif
      
      if (gamma_need) then
         if (.not. rtdb_put(rtdb, 'wfat:gamma_start', mt_dbl, 1,
     &        wfat_gamma_start))
     &        call errquit('WFAT: Failed writing gamma_start to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:gamma_last', mt_dbl, 1,
     &        wfat_gamma_last))
     &        call errquit('WFAT: Failed writing gamma_last to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:ngamma', mt_int, 1,
     &        wfat_ngamma))
     &        call errquit('WFAT: Failed writing ngamma to rtdb',
     &        1, INPUT_ERR)
      endif

      if (orbid_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:n_orbid', mt_int, 1, n_orbid))
     &        call errquit('WFAT: Failed writing n_orbid to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:orbid', mt_int, n_orbid,
     &        wfat_orbid0))
     &        call errquit('WFAT: Failed writing orbid to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:spinid', mt_int, n_orbid,
     &        wfat_spinid))
     &        call errquit('WFAT: Failed writing spinid to rtdb',
     &        1, INPUT_ERR)
      endif

      if (paraqn_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:n_npara', mt_int, 1, n_npara))
     &        call errquit('WFAT: Failed writing n_npara to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:npara', mt_int, n_npara, npara))
     &        call errquit('WFAT: Failed writing npara to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:mpara', mt_int, n_npara, mpara))
     &        call errquit('WFAT: Failed writing mpara to rtdb',
     &        1, INPUT_ERR)
      endif

      if (nktmd_need) then
         if (.not. rtdb_put(rtdb, 'wfat:nktmd', mt_int, 1, nktmd))
     &        call errquit('WFAT: Failed writing nktmd to rtdb',
     &        1, INPUT_ERR)
      endif

      if (beta_tmd_need) then
         if (.not. rtdb_put(rtdb, 'wfat:beta_tmd', mt_dbl, 1, beta_tmd))
     &        call errquit('WFAT: Failed writing beta_tmd to rtdb',
     &        1, INPUT_ERR)
      endif

      if (gamma_tmd_need) then
         if (.not. rtdb_put(rtdb, 'wfat:gamma_tmd', mt_dbl, 1,
     &        gamma_tmd))
     &        call errquit('WFAT: Failed writing gamma_tmd to rtdb',
     &        1, INPUT_ERR)
      endif

      if (lmax_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:lmax', mt_int, 1, wfat_lmax))
     &        call errquit('WFAT: Failed writing lmax to rtdb',
     &        1, INPUT_ERR)
      endif

cOLD      if (gqual_need .and. os) then
cOLD         if (.not. rtdb_put(rtdb, 'wfat:gqual', mt_int, 1, wfat_gqual))
cOLD     &        call errquit('WFAT: Failed writing gqual to rtdb',
cOLD     &        1, INPUT_ERR)
cOLD      endif
      
      if (omgpos_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:autopos', mt_log, 1, autopos))
     &        call errquit('WFAT: Failed writing autopos to rtdb',
     &        1, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:omgpos', mt_dbl, 3, omgpos))
     &        call errquit('WFAT: Failed writing omgpos to rtdb',
     &        1, INPUT_ERR)
      endif

      if (dgnthr_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:dgnthr', mt_dbl, 1, dgnthr))
     &        call errquit('WFAT: Failed writing dgnthr to rtdb',
     &        1, INPUT_ERR)
      endif
      
      if (intgl_need) then
         if (os .and. .not. or) then            
            if (wfat_intgl_name == 'nosave') then
               nosave = .true.
            else
               nosave = .false.
            endif
            
            if (.not. rtdb_put(rtdb, 'wfat:intgl_nosave', mt_log, 1,
     &           nosave))
     &           call errquit('WFAT: Failed writing intgl_nosave '//
     &           'to RTDB.', 1, INPUT_ERR)
            if (wfat_intgl_name /= 'nosave') then
               if (.not. rtdb_cput(rtdb, 'wfat:intgl_name', 1,
     &              trim(adjustl(wfat_intgl_name))))
     &              call errquit('WFAT: Failed writing intgl_name to '//
     &              'RTDB.', 1, INPUT_ERR)
            endif
         else if (.not. os .and. or) then
            if (wfat_intgl_name == 'nosave') then
               call errquit('WFAT: The startup mode is ''restart''. '//
     &              'For this startup type, the ''intgl'' directive '//
     &              'is only used to provide the file name '//
     &              'containing WFAT integrals to load from (not to '//
     &              'save to). Thus, the keyword ''nosave'' in this '//
     &              'regard is irrelevant. To prevent this error, '//
     &              'either provide a valid WFAT integral file name '//
     &              'or simply remove ''intgl'' directive from the '//
     &              'input file, for which the file name will '//
     &              'default to ''<file_prefix>.<m|o>wigl'' where '//
     &              'file_prefix is the file prefix of NWChem files '//
     &              'for the present calculation. The user is '//
     &              'responsible in assuring the existence of this '//
     &              'integral file.', 1, INPUT_ERR)
            else
               if (.not. rtdb_cput(rtdb, 'wfat:intgl_name', 1,
     &              trim(adjustl(wfat_intgl_name))))
     &              call errquit('WFAT: Failed writing intgl_name to '//
     &              'RTDB.', 1, INPUT_ERR)
            endif
         endif
      endif

      if (wfat_level == 'oe') then
         if (movec_need .and. os) then
            if (.not. rtdb_cput(rtdb, 'wfat:movec_name', 1,
     &           trim(adjustl(wfat_movec_name))))
     &           call errquit('WFAT: Failed writing movec_name to RTDB',
     &           1, RTDB_ERR)
         endif
      else if (wfat_level == 'me') then
         if (movec_need .and. os) then
            if (.not. rtdb_cput(rtdb, 'mewfat:movec0_name', 1,
     &           trim(adjustl(wfat_movec_name))))
     &           call errquit('WFAT: Failed writing movec0_name to '//
     &           'RTDB', 1, RTDB_ERR)
         endif
         if (movec1_need .and. os) then
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

      if (sitype_need .and. os) then
         if (.not. rtdb_cput(rtdb, 'mewfat:si', 1,
     &        trim(adjustl(si_type))))
     &        call errquit('WFAT: Failed writing si to RTDB', 1,
     &        RTDB_ERR)
      endif

      if (wfat1_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:wfat1', mt_log, 1, wfat1))
     &        call errquit('WFAT: Failed writing wfat1 to RTDB', 1,
     &        RTDB_ERR)
      endif

      if (w1x_need .and. os) then
         nw1x_ = n_npara * n_orbid
         if (.not. rtdb_put(rtdb, 'wfat:wfat1_exclude', mt_log, nw1x_,
     &        wfat1_exclude_l(1:nw1x_)))
     &        call errquit('WFAT: Failed writing wfat1_exclude to RTDB',
     &        1, RTDB_ERR)
      endif

      if (mdip_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:use_mandip', mt_log, n_orbid,
     &        wfat_mdip_l(1:n_orbid)))
     &        call errquit('WFAT: Failed writing use_mandip to RTDB',
     &        1, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:mandip', mt_dbl, 3*n_orbid,
     &        wfat_mdip(1:3, 1:n_orbid)))
     &        call errquit('WFAT: Failed writing mandip to RTDB',
     &        1, RTDB_ERR)
      endif

      if (merg_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:use_manerg', mt_log, n_orbid,
     &        wfat_merg_l(1:n_orbid)))
     &        call errquit('WFAT: Failed writing use_manerg to RTDB',
     &        1, RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'wfat:manerg', mt_dbl, n_orbid,
     &        wfat_merg(1:n_orbid)))
     &        call errquit('WFAT: Failed writing manerg to RTDB',
     &        1, RTDB_ERR)
      endif

      if (nlgl_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:nlgl', mt_int, 1, wfat_nlgl))
     &        call errquit('WFAT: Failed writing nlgl to RTDB', 1,
     &        RTDB_ERR)
      endif

      if (dthr_need .and. os) then
         if (.not. rtdb_put(rtdb, 'wfat:dthr', mt_dbl, 1, wfat_dthr))
     &        call errquit('WFAT: Failed writing wfat_dthr to RTDB', 1,
     &        RTDB_ERR)
      endif

      if (wfat_mask_need .and. os) then
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


