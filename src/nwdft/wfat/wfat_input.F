c=======================================================
c=======================================================
      subroutine wfat_input(rtdb)

      implicit none
#include "stdio.fh"      
#include "global.fh"      
      integer rtdb

      if (ga_nodeid() .ne. 0) return      
      
c     ==== Read WFAT inputs from input file ====
      call wfat_info_banner(luout, .true.)
      write(luout, '(5x, "Reading the WFAT inputs ...")')
      call wfat_read_input(rtdb)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_read_input(rtdb)

      implicit none
#include "errquit.fh"
#include "wfat_input.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"      
      integer rtdb
      
      logical tt, ll1, ll2, ll3
      character*50 field, tail_unit
      character*255 ss
      integer i1, i2, k
      double precision ang2bohr

      parameter(ang2bohr=1.8897261255d0)
      
      logical wfat_inp_paraqn, wfat_inp_movec, wfat_inp_omgpos
      logical wfat_inp_field, wfat_inp_perturb, wfat_inp_dt
      logical wfat_inp_tmax, wfat_inp_ffc, wfat_inp_w1x, wfat_inp_s2a1e
      logical wfat_inp_mydip, wfat_inp_myerg
      external wfat_inp_paraqn, wfat_inp_movec, wfat_inp_omgpos
      external wfat_inp_field, wfat_inp_perturb, wfat_inp_dt
      external wfat_inp_tmax, wfat_inp_ffc, wfat_inp_w1x, wfat_inp_s2a1e
      external wfat_inp_mydip, wfat_inp_myerg

      
      WFAT_PRINT_OPTIONS = (/ 'srfc', 'rate', 'yield', 'avbeta', 
     &                        'avgamma', 'tmd', 'tmd_atom' /)


c     ==== Initial status of WFAT inputs ====
      print_exist = .false.
      print_stat = .false.

      nfield = 0
      fld_tp = 'n/a'        ! Initiation.
      field_exist = .false.
      field_stat = .false.

      dt_exist = .false.
      dt_stat = .false.
      
      tmax_exist = .false.
      tmax_stat = .false.
      
      beta_exist = .false.
      beta_stat = .false.

      gamma_exist = .false.
      gamma_stat = .false.
      
      orbid_exist = .false.
      orbid_stat = .false.

      paraqn_exist = .false.
      paraqn_stat = .false.
      
      nktmd_exist = .false.
      nktmd_stat = .false.

      beta_tmd_exist = .false.
      beta_tmd_stat = .false.

      gamma_tmd_exist = .false.
      gamma_tmd_stat = .false.

      lmax_exist = .false.
      lmax_stat = .false.

cOLD      gqual_exist = .false.
cOLD      gqual_stat = .false.

      omgpos_exist = .false.
      omgpos_stat = .false.

      sitype_exist = .false.
      sitype_stat = .false.

      dgnthr_exist = .false.
      dgnthr_stat =  .false.
      
      intgl_exist = .false.
      intgl_stat = .false.

      movec_exist = .false.
      movec_stat = .false.
      movec1_exist = .false.
      movec1_stat = .false.

      wfat1_exist = .false.
      wfat1_stat = .false.

      w1x_exist = .false.
      w1x_stat = .false.

      mdip_exist = .false.
      mdip_stat = .false.

      merg_exist = .false.
      merg_stat = .false.
      
      nlgl_exist = .false.
      nlgl_stat = .false.

      dthr_exist = .false.
      dthr_stat = .false.

      ffc_exist = .false.
      ffc_stat = .false.

      wfat_mask_exist = .false.
      wfat_mask_stat = .false.

      perturb_exist = .false.
      perturb_stat = .false.
      
      prtb_type_exist = .false.
      prtb_type_stat = .false.
      
      prtb_ord_exist = .false.
      prtb_ord_stat = .false.
      
      prtb_bt_exist = .false.
      prtb_bt_stat = .false.
      
      prtb_gm_exist = .false.
      prtb_gm_stat = .false.      

      prtb_dgnthr_exist = .false.
      prtb_dgnthr_stat = .false. 
      
      prtb_offthr_exist = .false.
      prtb_offthr_stat = .false. 
      
      prtb_e1thr_exist = .false. 
      prtb_e1thr_stat = .false.
      
      prtb_itrthr_exist = .false.
      prtb_itrthr_stat = .false. 
      
      prtb_maxitr_exist = .false.
      prtb_maxitr_stat = .false. 
      
      prtb_ordthr_exist = .false.
      prtb_ordthr_stat = .false. 
      
      prtb_maxord_exist = .false.
      prtb_maxord_stat = .false.

      print_srfc = .false.
      print_rate = .false.
      print_yild = .false.
      print_avbeta = .false.
      print_avgamma = .false.
      print_tmd = .false.
      print_tmd_atom = .false.

      s2a1e_exist = .false.
      s2a1e_stat = .true.

      
c     ========================================================
c     ==================== INPUTS PARSING ====================
c     ========================================================

      do while (inp_read()) 
         if (inp_a(field)) then

cfordebug            write(luout, '("input field : ", a)') field

c           ==== ABSOLUTELY MANDATORY INPUTS ====
            if (inp_compare(.false., 'print', field)) then
               print_exist = .true.
               
               if (.not. inp_line(ss))
     &              call errquit('WFAT: Error in reading the '//
     &              '''print'' line.', 1, INPUT_ERR)
               i1 = 0       ! i1 must be initiated to zero before first call of inp_strtok
               tt = inp_strtok(ss, ' ', i1, i2) !! The 'print' keyword
               
               do while (inp_strtok(ss, ' ', i1, i2))        !! The rest of 'print' line
                  tt = inp_match(WFAT_NPRINT, .false.,
     &                 ss(i1:i2), WFAT_PRINT_OPTIONS, k)
                  
                  if (WFAT_PRINT_OPTIONS(k) == 'srfc') then
                     print_srfc = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'rate') then
                     print_rate = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'yield') then
                     print_yild = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'avbeta') then
                     print_avbeta = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'avgamma') then
                     print_avgamma = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'tmd') then
                     print_tmd = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'tmd_atom') then
                     print_tmd_atom = .true.
                  else
                     print_stat = .false.
                     exit
                  endif
                  print_stat = .true.
               enddo

c           ==== CONDITIONALLY MANDATORY INPUTS ====    
cOLD            else if (inp_compare(.false., 'fieldmag', field)) then
cOLD               fieldmag_exist = .true.
cOLD               fieldmag_stat = inp_f(fieldmag)

            else if (inp_compare(.false., 'field', field)) then
               nfield = nfield + 1
               field_exist = .true.
               field_stat = wfat_inp_field(nfield, fld_tp, 
     &              fld_fr, fld_fr0, fld_fr_un, 
     &              fld_ce, fld_ce0, fld_ce_un, 
     &              fld_wi, fld_wi0, fld_wi_un,
     &              fld_mx, fld_mx0, fld_mx_un,
     &              fld_cp, fld_cp0, fld_cp_un)
               if (is_extfield) then
                  dt_exist = .true.
                  dt_stat = .true.
                  tmax_exist = .true.
                  tmax_stat = .true.
               endif

            else if (inp_compare(.false., 'dt', field)) then
               dt_exist = .true.
               dt_stat = wfat_inp_dt(wfat_dt, wfat_dt0, dt_un)

            else if (inp_compare(.false., 'tmax', field)) then
               tmax_exist = .true.
               tmax_stat = wfat_inp_tmax(wfat_tmax, wfat_tmax0, tmax_un)
               
            else if (inp_compare(.false., 'beta', field)) then
               beta_exist = .true.
               beta_start_stat = inp_f(wfat_beta_start)
               beta_n_stat = inp_i(wfat_nbeta)
               if (.not. beta_n_stat .or. wfat_nbeta==1) then
                  wfat_nbeta = 1
                  wfat_beta_last = wfat_beta_start
                  beta_last_stat = .true.
               else
                  beta_last_stat = inp_f(wfat_beta_last)
               endif
               beta_stat = beta_start_stat 
            else if (inp_compare(.false., 'gamma', field)) then
               gamma_exist = .true.
               gamma_start_stat = inp_f(wfat_gamma_start)
               gamma_n_stat = inp_i(wfat_ngamma)
               if (.not. gamma_n_stat .or. wfat_ngamma==1) then
                  wfat_ngamma = 1
                  wfat_gamma_last = wfat_gamma_start
                  gamma_last_stat = .true.
               else
                  gamma_last_stat = inp_f(wfat_gamma_last)
               endif
               gamma_stat = gamma_start_stat
            else if (inp_compare(.false., 'beta_tmd', field)) then
               beta_tmd_exist = .true.
               beta_tmd_stat = inp_f(beta_tmd)
            else if (inp_compare(.false., 'gamma_tmd', field)) then
               gamma_tmd_exist = .true.
               gamma_tmd_stat = inp_f(gamma_tmd)
               
c           ==== CONDITIONALLY OPTIONAL INPUTS ====
            else if (inp_compare(.false., 'nktmd', field)) then
               nktmd_exist = .true.
               nktmd_stat = inp_i(nktmd)
               
c           ==== OPTIONAL INPUTS ====
            else if (inp_compare(.false., 'movec_id', field)) then
               orbid_exist = .true.
               orbid_stat = wfat_inp_movec(rtdb, MAX_N_ORBID, 
     &                      wfat_orbid0, wfat_spinid, n_orbid)
            else if(inp_compare(.false., 'paraqn', field)) then
               paraqn_exist = .true.
               paraqn_stat = wfat_inp_paraqn(MAX_N_NPARA, npara, mpara,
     &                       n_npara)
            else if (inp_compare(.false., 'lmax', field)) then
               lmax_exist = .true.
               lmax_stat = inp_i(wfat_lmax)
            else if (inp_compare(.false., 'quadgrid', field)) then
               call wfat_parse_quadgrid(rtdb)
cOLD            else if (inp_compare(.false., 'gridqual', field)) then
cOLD               gqual_exist = .true.
cOLD               gqual_stat = inp_i(wfat_gqual)
            else if (inp_compare(.false., 'omgpos', field)) then
               omgpos_exist = .true.
               omgpos_stat = wfat_inp_omgpos(autopos, omgpos)
            else if (inp_compare(.false., 'si_type', field)) then
               sitype_exist = .true.
               sitype_stat = inp_a(si_type)
            else if (inp_compare(.false., 'degenthr', field)) then
               dgnthr_exist = .true.
               dgnthr_stat = inp_f(dgnthr)
            else if (inp_compare(.false., 'intgl', field)) then
               intgl_exist = .true.
               intgl_stat = inp_a(wfat_intgl_name)
            else if (inp_compare(.false., 'vectors', field)) then   ! Needed for restarting MEWFAT task.
               movec_exist = .true.
               movec_stat = inp_a(wfat_movec_name)
               if (inp_a(wfat_movec1_name)) then
                  movec1_exist = .true.
                  movec1_stat = .true.
               endif
               call inp_prev_field()
            else if (inp_compare(.false., 'wfat1', field)) then
               wfat1_exist = .true.
               wfat1_stat = .true.
               wfat1 = .true.
            else if (inp_compare(.false., 'wfat1_exclude', field)) then
               w1x_exist = .true.
               w1x_stat = wfat_inp_w1x(CHID_CHARLEN,
     &              MAX_N_NPARA*MAX_N_ORBID, wfat1_exclude, w1x_len)
            else if (inp_compare(.false., 'mydipole', field)) then
               mdip_exist = .true.
               mdip_stat = wfat_inp_mydip(ORBID_CHARLEN,
     &              MAX_N_ORBID, wfat_mdip_s, wfat_mdip, mdip_len)
            else if (inp_compare(.false., 'myenergy', field)) then
               merg_exist = .true.
               merg_stat = wfat_inp_myerg(ORBID_CHARLEN,
     &              MAX_N_ORBID, wfat_merg_s, wfat_merg, merg_len)
            else if (inp_compare(.false., 'mychpos', field)) then
               call wfat_inp_mypos(rtdb, MAX_N_ORBID)
            else if (inp_compare(.false., 'nlgl', field)) then
               nlgl_exist = .true.
               nlgl_stat = inp_i(wfat_nlgl)
            else if (inp_compare(.false., 'diagthr', field)) then
               dthr_exist = .true.
               dthr_stat = inp_f(wfat_dthr)
            else if (inp_compare(.false., 'tailmask', field)) then
               wfat_mask_exist = .true.
               ll1 = inp_f(wfat_mask_param(1))
               ll2 = inp_f(wfat_mask_param(2))
               ll3 = inp_a(tail_unit)
               if (.not. ll3 .or. inp_compare(.false., 'angstrom',
     &              tail_unit)) then
                  wfat_mask_param = wfat_mask_param * ang2bohr
               else if (inp_compare(.false., 'bohr', tail_unit)) then
                  !do nothing
               else
                  call errquit('WFAT: The option for the unit of '//
     &                 'tailmask is angstrom or bohr', 90, INPUT_ERR)
               endif
               wfat_mask_stat = (ll1 .and. ll2)
               if (wfat_mask_stat) wfat_mask = .true.
            else if (inp_compare(.false., 'field_factor', field)) then
               ffc_exist = .true.
               ffc_type0 = -9999
               ffc_stat = wfat_inp_ffc(ffc_type0, ergffc)
            else if (inp_compare(.false., 'perturb', field)) then
               perturb_exist = .true.
               perturb_stat = wfat_inp_perturb(
     &              prtb_type,   prtb_type_exist,   prtb_type_stat,
     &              prtb_ord,    prtb_ord_exist,    prtb_ord_stat,
     &              prtb_bt,     prtb_bt_exist,     prtb_bt_stat,
     &              prtb_gm,     prtb_gm_exist,     prtb_gm_stat,
     &              prtb_dgnthr, prtb_dgnthr_exist, prtb_dgnthr_stat,
     &              prtb_offthr, prtb_offthr_exist, prtb_offthr_stat,
     &              prtb_e1thr,  prtb_e1thr_exist,  prtb_e1thr_stat,
     &              prtb_itrthr, prtb_itrthr_exist, prtb_itrthr_stat,
     &              prtb_maxitr, prtb_maxitr_exist, prtb_maxitr_stat,
     &              prtb_ordthr, prtb_ordthr_exist, prtb_ordthr_stat,
     &              prtb_maxord, prtb_maxord_exist, prtb_maxord_stat)
            else if (inp_compare(.false., 'system_2a1e', field)) then
               s2a1e_exist = .true.
               s2a1e_stat = wfat_inp_s2a1e(wfat2a_nc, wfat2a_q, 
     &              wfat2a_inter, wfat2a_mag, wfat2a_xzr, wfat2a_pot, 
     &              wfat2a_potparam, wfat2a_nel, wfat2a_fedvr, 
     &              wfat2a_nfe, wfat2a_fenodes, wfat2a_ndvr, wfat2a_nr, 
     &              wfat2a_kein, wfat2a_kemult, wfat2a_maxit, wfat2a_nt)
               s2a1e = s2a1e_stat

c           ==== End ====
            else if (inp_compare(.false., 'end', field)) then
               exit
            else
               call errquit('WFAT: The specified keyword '''
     &              //trim(adjustl(field))//''' is unavailable.', 1,
     &              INPUT_ERR)
            endif
         else
            call errquit('WFAT: The first keyword of an input line '//
     &           'must be a valid character type.',
     &           1, INPUT_ERR)
         endif
      enddo


c     ==== Obtain the theory of previous task ====
      if (s2a1e) then
         wfat_theory = '2a1e'
      else
         if (rtdb_cget(rtdb,'scf:scftype', 1, wfat_theory)) then
            call wfat_info_banner(luout, .true.)
            write(luout, '(5x, "Previous SCF type = ", a)')
     &           trim(wfat_theory)
            write(luout, '(5x, "If this is a ME-WFAT job, then this '//
     &           'information refers to the last charge state ", /,'//
     &           '5x, "SCF calculation (either the neutral or , '//
     &           'cation ME-WFAT needs two prior SCF calculations).")')
         else
            call errquit('task: no task input for theory?',0, RTDB_ERR)
         endif
      endif
      if (.not. rtdb_cput(rtdb, 'wfat:theory', 1, trim(wfat_theory)))
     &     call errquit('WFAT: Error writing theory to rtdb',
     &     1, RTDB_ERR)

      
c     ==== Check whether the current theory is DFT ====
      if (trim(adjustl(wfat_theory)) /= 'RHF' .and.
     &    trim(adjustl(wfat_theory)) /= 'UHF' .and.
     &    .not. s2a1e) then
         write(luout, '("WFAT: WARNING!")')
         write(luout, 100)
      endif
 100  format('      WFAT input directive is found while the SCF type ',
     &       'is not RHF. At the moment WFAT can only be used in ',
     &       'conjunction with single determinant methods to obtain ',
     &       'quantities such as molecular orbitals and effective ',
     &       'potential. All WFAT related inputs will thus be ignored.')

      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_field(i, tp, fr, fr0, fr_un,
     &     ce, ce0, ce_un, wi, wi0, wi_un, mx, mx0, mx_un,
     &     cp, cp0, cp_un)

c     field
c        type < dc | gaussian | hann | external <"laser_file".laser> >
c        frequency <double> [ au | ev | thz | nm | bohr ]
c        center <double> [ fs | au ]
c        width <double> [ fs | au ]
c        max <double> [ au | tw/cm2 ]
c        cep <double> [ rad | deg ]
c     end
      implicit none
#include "wfat_input.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      integer i                                        ! [input],
      character(len=FLD_TP_LEN) tp(max_nfield)         ! [output],
      double precision fr(max_nfield)                  ! [output],
      double precision fr0(max_nfield)                 ! [output],
      character(len=FLD_FR_UN_LEN) fr_un(max_nfield)   ! [output],
      double precision ce(max_nfield)                  ! [output],
      double precision ce0(max_nfield)                 ! [output],
      character(len=FLD_CE_UN_LEN) ce_un(max_nfield)   ! [output],
      double precision wi(max_nfield)                  ! [output],
      double precision wi0(max_nfield)                 ! [output],
      character(len=FLD_WI_UN_LEN) wi_un(max_nfield)   ! [output],
      double precision mx(max_nfield)                  ! [output],
      double precision mx0(max_nfield)                 ! [output],
      character(len=FLD_MX_UN_LEN) mx_un(max_nfield)   ! [output],
      double precision cp(max_nfield)                  ! [output],
      double precision cp0(max_nfield)                 ! [output],
      character(len=FLD_CP_UN_LEN) cp_un(max_nfield)   ! [output],

      character*50 field
      character*2 max_str
      logical tp_stat, fr_stat, ce_stat, wi_stat, mx_stat, cp_stat
      logical tp_exist, fr_exist, ce_exist, wi_exist, mx_exist, cp_exist
      integer k

      logical wfat_inp_externfield
      double precision wfat_convert_freq, wfat_convert_time,
     &     wfat_convert_efield, wfat_convert_angle
      external wfat_convert_freq, wfat_convert_time,
     &     wfat_convert_efield, wfat_convert_angle, wfat_inp_externfield


      wfat_inp_field = .false.
      is_extfield = .false.
      tp_exist = .false.
      fr_exist = .false.
      ce_exist = .false.
      wi_exist = .false.
      mx_exist = .false.
      cp_exist = .false.


      TP_OPTS = (/ 'dc', 'gaussian', 'hann', 'external' /)
      FR_UNIT = (/ 'au', 'ev', 'thz', 'nm', 'bohr' /)
      CE_UNIT = (/ 'fs', 'au' /)
      WI_UNIT = (/ 'fs', 'au' /)
      MX_UNIT = (/ 'au', 'tw/cm2' /)
      CP_UNIT = (/ 'rad', 'deg' /)
      

      write(max_str, '(i2)') max_nfield
      if (i > max_nfield)
     &     call errquit('WFAT: No. of fields exceeds the maximum, '//
     &     'which is '//trim(adjustl(max_str)), 1, INPUT_ERR)
     

      do while (inp_read())
         if (inp_a(field)) then

            if (inp_compare(.false., 'type', field)) then
               tp_exist = .true.
               tp_stat = inp_a(tp(i))
               if (.not. inp_match(TP_NOPTS, .false., tp(i), TP_OPTS,
     &              k)) call errquit('WFAT: The choice of field type '//
     &              'is unknown', 1, INPUT_ERR)
               if (tp(i) == 'external') then
                  if (i > 1) call errquit('WFAT: When there is a '//
     &                 'field with the type of ''external'', that '//
     &                 'field is the only field that can be specified'//
     &                 ' in the input file.', 1, INPUT_ERR)
                  is_extfield = wfat_inp_externfield()
               endif
               
            else if (inp_compare(.false., 'frequency', field)) then
               fr_exist = .true.
               fr_stat = inp_f(fr0(i))
               fr_un(i) = 'au'
               if (inp_a(field)) fr_un(i) = field
               fr(i) = wfat_convert_freq(fr_un(i), fr0(i), 'au')
               
            else if (inp_compare(.false., 'center', field)) then
               ce_exist = .true.
               ce_stat = inp_f(ce0(i))
               ce_un(i) = 'au'
               if (inp_a(field)) ce_un(i) = field
               ce(i) = wfat_convert_time(ce_un(i), ce0(i), 'au')
               
            else if (inp_compare(.false., 'width', field)) then
               wi_exist = .true.
               wi_stat = inp_f(wi0(i))
               wi_un(i) = 'au'
               if (inp_a(field)) wi_un(i) = field
               wi(i) = wfat_convert_time(wi_un(i), wi0(i), 'au')
               
            else if (inp_compare(.false., 'max', field)) then
               mx_exist = .true.
               mx_stat = inp_f(mx0(i))
               mx_un(i) = 'au'
               if (inp_a(field)) mx_un(i) = field
               mx(i) = wfat_convert_efield(mx_un(i), mx0(i), 'au')
               
            else if (inp_compare(.false., 'cep', field)) then
               cp_exist = .true.
               cp_stat = inp_f(cp0(i))
               cp_un(i) = 'rad'
               if (inp_a(field)) cp_un(i) = field
               cp(i) = wfat_convert_angle(cp_un(i), cp0(i), 'rad')

            else if (inp_compare(.false., 'end', field)) then
               exit
            else
               call errquit('WFAT: Unknown input directive '''//
     &              trim(adjustl(field))//''' is found under the '//
     &              '''field'' subdirective', 1, INPUT_ERR)
            endif

         endif
      enddo
      

c     ==== Checking unit choice availability ====
      if (.not. tp_exist) then
         call errquit('WFAT: The ''field'' subdirective '//
     &        '''type'' is missing', 1, INPUT_ERR)
      endif
      
      if (tp(i) /= 'external') then
         if (mx_exist) then
            if (.not. inp_match(MX_NUNIT, .false., mx_un(i),
     &           MX_UNIT, k))
     &           call errquit('WFAT: The unit choice for ''field'' '//
     &           'subdirective ''max'' is unavailable', 1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''max'' is missing', 1, INPUT_ERR)
         endif
      endif

      if (tp(i) /= 'dc' .and. tp(i) /= 'external') then
         if (fr_exist) then
            if (.not. inp_match(FR_NUNIT, .false., fr_un(i), FR_UNIT,
     &           k)) call errquit('WFAT: The unit choice for ''field'' '
     &           //'subdirective ''frequency'' is unavailable',
     &           1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''frequency'' is missing when it is necessary',
     &           1, INPUT_ERR)
         endif

         if (ce_exist) then
            if (.not. inp_match(CE_NUNIT, .false., ce_un(i), CE_UNIT, 
     &           k)) call errquit('WFAT: The unit choice for ''field'' '
     &           //'subdirective ''center'' is unavailable',
     &           1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''center'' is missing when it is necessary',
     &           1, INPUT_ERR)
         endif

         if (wi_exist) then
            if (.not. inp_match(WI_NUNIT, .false., wi_un(i), WI_UNIT, 
     &           k)) call errquit('WFAT: The unit choice for ''field'' '
     &           //'subdirective ''width'' is unavailable',
     &           1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''width'' is missing when it is necessary',
     &           1, INPUT_ERR)
         endif

         if (cp_exist) then
            if (.not. inp_match(CP_NUNIT, .false., cp_un(i), CP_UNIT, 
     &           k)) call errquit('WFAT: The unit choice for ''field'' '
     &           //'subdirective ''cep'' is unavailable', 1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''cep'' is missing when it is necessary',
     &           1, INPUT_ERR)
         endif
      endif
      

c     ==== When field type is 'dc' ====
      if (tp(i) == 'dc' .or. tp(i) == 'external') then
         fr0(i) = 0.0d0
         ce0(i) = 0.0d0
         wi0(i) = 0.0d0
         cp0(i) = 0.0d0
         fr(i) = 0.0d0
         ce(i) = 0.0d0
         wi(i) = 0.0d0
         cp(i) = 0.0d0
      endif
      if (tp(i) == 'external') then
         mx(i) = 0.0d0
         mx0(i) = 0.0d0
      endif
      

      if (tp(i) == 'dc') then
         wfat_inp_field = tp_stat .and. mx_stat
      else if (tp(i)=='gaussian' .or. tp(i)=='hann') then
         wfat_inp_field = tp_stat .and. fr_stat .and. ce_stat .and.
     &                    wi_stat .and. mx_stat .and. cp_stat
      else if (tp(i) == 'external') then
         wfat_inp_field = is_extfield
      else
         wfat_inp_field = .false.
      endif
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      logical function wfat_inp_externfield()
c     <string t_unit>
c     <string f_unit>
c     <integer data_start>  <integer data_end>
c     data:
c     <double t(1)>     <double f(1)>
c     <double t(2)>     <double f(2)>
c     ...               ...          
     
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "inp.fh"
#include "wfat_input.fh"
      character(len=EXT_TM_UN_LEN) t_unit
      character(len=EXT_FL_UN_LEN) f_unit
      character*5 ss
      integer k, h, i, j, ios
      double precision time, field

      double precision wfat_convert_time, wfat_convert_efield
      external wfat_convert_time, wfat_convert_efield


      wfat_inp_externfield = .false.
      TMEXT_UNIT = (/ 'fs', 'au' /)
      FLEXT_UNIT = (/ 'au', 'tw/cm2' /)
      fld_t_ext(:) = 0.0d0
      fld_f_ext(:) = 0.0d0

      
      if (inp_a(extfield_name)) then
         open(h, file=extfield_name, action='read')

c        ==== Time data unit ====
         read(h, '(a)') t_unit
         if (.not. inp_match(TMEXT_NUNIT, .false., t_unit, TMEXT_UNIT, 
     &        k)) call errquit('WFAT: The unit choice for time data of '
     &        //'the external field is unavailable', 1, INPUT_ERR)

c        ==== Field data unit ====
         read(h, '(a)') f_unit
         if (.not. inp_match(FLEXT_NUNIT, .false., f_unit, FLEXT_UNIT,
     &        k))
     &        call errquit('WFAT: The unit choice for field data of '//
     &        'the external field is unavailable', 1, INPUT_ERR)

c        ==== Data range ====
         read(h, *) ext_start, ext_end
         ntime_ext = ext_end - ext_start + 1
         if (ntime_ext > MAX_EXTERNFIELD) then
            write(luout, 100) ntime_ext, MAX_EXTERNFIELD
            call errquit('WFAT: Error inside the external field file.',
     &           1, INPUT_ERR)
         endif

c        ==== Field and time data ====
         read(h, '(a)') ss
         if (ss == 'data:') then
            j = 0
            do i = 1, MAX_EXTERNFIELD
               if (i >= ext_start .and. i <= ext_end) then
                  read(h, *, iostat=ios) time, field
                  if (ios /= 0) call errquit('WFAT: Bad data within '//
     &                 'the ''data:'' scope in the file containing'//
     &                 'the external field.', 1, INPUT_ERR)
                  j = j + 1
                  fld_t_ext(j) = wfat_convert_time(t_unit, time, 'au')
                  fld_f_ext(j) = wfat_convert_efield(f_unit, field,
     &                 'au')
               else if (i < ext_start) then
                  read(h, *)
               else if (i > ext_end) then
                  exit
               endif
            enddo
         else
            call errquit('WFAT: The keyword ''data:'' is not found '//
     &           'in the external field file.', 1, INPUT_ERR)
         endif
         
         close(h)

c        ==== dt and tmax ====
         dt_ext = fld_t_ext(2) - fld_t_ext(1)
         write(luout, '("dt_ext = ", es14.6)') dt_ext
         tmax_ext = fld_t_ext(ntime_ext) - fld_t_ext(1)
      else
         call errquit('WFAT: Filename of the file containing the '//
     &        'external field is incorrectly or not specified.',
     &        1, INPUT_ERR)
      endif


      write(luout, '(/, 3x, "*** ATTENTION ***")')
      write(luout, 101)

      wfat_inp_externfield = .true.

      
 100  format('ERROR: The no. of data points specified for the ',
     &     'external field, which is ', i6, /,
     &     7x, ', exceeds the maximum, which is ', i6, '.')

 101  format(
     &     5x, 'The field has been taken from an external file.',
     &     'In this case, the user is ', /,
     &     5x, 'responsible in ensuring that the time data values ',
     &     'are equally spaced.')
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_dt(dt, dt0, un)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "wfat_input.fh"
      double precision dt, dt0
      character(len=DT_UN_LEN) un

      integer k
      double precision wfat_convert_time
      external wfat_convert_time

      
      DT_UNIT = (/ 'au', 'fs' /)

      wfat_inp_dt = inp_f(dt0)
      if (inp_a(un)) then
         dt = wfat_convert_time(un, dt0, 'au')
         if (.not. inp_match(DT_NUNIT, .false., un, DT_UNIT, k)) 
     &        call errquit('WFAT: The unit choice for the directive '//
     &        '''dt'' is unavailable', 1, INPUT_ERR)
      else
         dt = dt0
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_tmax(tmax, tmax0, un)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "wfat_input.fh"
      double precision tmax, tmax0
      character(len=TMAX_UN_LEN) un

      integer k
      double precision wfat_convert_time
      external wfat_convert_time

      
      TMAX_UNIT = (/ 'au', 'fs' /)

      wfat_inp_tmax = inp_f(tmax0)
      if (inp_a(un)) then
         tmax = wfat_convert_time(un, tmax0, 'au')
         if (.not. inp_match(TMAX_NUNIT, .false., un, TMAX_UNIT, k)) 
     &        call errquit('WFAT: The unit choice for the directive '//
     &        '''tmax'' is unavailable', 1, INPUT_ERR)
      else
         tmax = tmax0
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_movec(rtdb, nelm, orbid0, spinid0, nor0)

      implicit none
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "wfat_input.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer rtdb
      integer nelm
      integer orbid0(nelm)
      integer spinid0(nelm)
      integer nor0

      character*4 nelm_str
      character*1 afield
      integer i, j, ior, ior1, ior2, is, io
      integer ifield
      integer orbid0_(nelm), spinid0_(nelm)

      logical wfat_check_invec
      external wfat_check_invec



      wfat_inp_movec = .false.

      ior = 1
      is = 0
c     ==== Loop for reading spins ====
      do i = 1, 2
         
         if (.not. inp_a(afield)) then
            if (is .ne. 0) then
               exit
            else
               call errquit('WFAT: No molecular orbitals to be '//
     &              'ionized have been specified. Specify at '//
     &              'least one orbital using the '//
     &              'format ''x 0 -1 -2 ...'' where ''x'' is either '//
     &              '''a'' or ''b''. The trailing list of integers '//
     &              'give the spatial orbitals having the x spin. At '//
     &              'the moment, only negative integers '//
     &              '(corresponding to occupied orbitals) are ',
     &              'supported.',1, INPUT_ERR)
            endif
         endif

         is = is + 1
         if (afield == 'a') then
            spinid0_(is) = 1
         else if (afield == 'b') then
            spinid0_(is) = 2
         else
            call errquit('WFAT: The spin specification for '//
     &           'input keyword ''movec_id'' can only take '//
     &           'either ''a'' or ''b''.', 1, INPUT_ERR)
         endif
         
         io = 0
c        ==== Loop for reading orbitals ====
         do j = 1, MAX_MOVEC_FIELDS
            if (inp_i(ifield)) then
               io = io + 1
               orbid0_(io) = ifield
            else
               if (io == 0)
     &              call errquit('WFAT: One spin channel has not '//
     &              'been given any orbitals.', 1, INPUT_ERR)
               exit
            endif
         enddo

         ior1 = ior
         ior2 = ior1 + io - 1
         nor0 = ior2

c     ==== Check duplicate in the spins ====
         if (is > 1) then
            if (wfat_check_invec(spinid0_(is), is-1, 
     &           spinid0_(1:is-1))) 
     &           call errquit('WFAT: Both of the specified spin '//
     &           'channels cannot be identical.', 1, INPUT_ERR)
         endif
         spinid0(ior1:ior2) = spinid0_(is)

c        ==== Check duplicate in orbital ====
         if (io > 1) then
            if (wfat_check_invec(orbid0_(io), io-1, 
     &           orbid0_(1:io-1)))
     &           call errquit('WFAT: Two of the specified spatial '//
     &           'orbitals that correspond to '//afield//
     &           '-spin are identical. '//
     &           'Each of the spatial orbitals for any given '//
     &           'spin channel must occur only once.',
     &           1, INPUT_ERR)
         endif
         orbid0(ior1:ior2) = orbid0_(1:io)
         

         ior = ior + io
         if (ior > nelm) then
            write(nelm_str, '(i4)')  nelm
            call errquit('WFAT: The number of the specified '//
     &           'molecular orbitals exceeds limit, which '//
     &           'is equal to '//nelm_str//'.', 1, INPUT_ERR)
         endif

      enddo

      if (.not. rtdb_put(rtdb, 'wfat:inp:nch', mt_int, 1, nor0))
     &     call errquit('WFAT: Failed writing wfat:inp:nch to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:orbid', mt_int, nor0, orbid0))
     &     call errquit('WFAT: Failed writing inp:orbid to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:spinid', mt_int, nor0, 
     &     spinid0)) call errquit('WFAT: Failed writing inp:spinid to'//
     &     ' rtdb', 1, INPUT_ERR)

      wfat_inp_movec = .true.

      end
c=======================================================

      
c=======================================================
c=======================================================
      logical function wfat_inp_paraqn(nelm, n0, m0, nqn0)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "wfat_input.fh"
#include "inp.fh"      
      integer nelm
      integer n0(nelm)
      integer m0(nelm)
      integer nqn0

      character*4 nelm_str, n_str
      character*5 afield
      logical cond1, cond2
      integer i, j, iqn, iqn1, iqn2, i_n, i_m
      integer ifield
      integer n0_(nelm), m0_(nelm)


      logical wfat_check_invec
      external wfat_check_invec



      wfat_inp_paraqn = .false.

      iqn = 1
      i_n = 0
c     ==== Loop for reading n quantum numbers ====
      do i = 1, MAX_PARAQN_NFIELDS
         
         if (.not. inp_a(afield)) then
            if (i_n .ne. 0) then
               exit
            else
               call errquit('WFAT: No n-parabolic quantum '//
     &              'number has been specified. Specify at '//
     &              'least one n-quantum number using the '//
     &              'format ''nx'' where ''x'' is zero or '//
     &              'positive integer.',1, INPUT_ERR)
            endif
         endif

         cond1 = (afield(1:1) == 'n')
         cond2 = (verify(trim(adjustl(afield(2:NMAX_NFIELD))),
     &           '0123456789') == 0)
         if (cond1 .and. cond2) then
            i_n = i_n + 1
            read(afield(2:NMAX_NFIELD), NFIELD_FORM) n0_(i_n)

            i_m = 0
c           ==== Loop for reading m quantum numbers ====
            do j = 1, MAX_PARAQN_MFIELDS
               if (inp_i(ifield)) then
                  i_m = i_m + 1
                  m0_(i_m) = ifield
               else
                  if (i_m == 0)
     &                 call errquit('WFAT: One n-parabolic quantum '//
     &                 'number has no m-parabolic quantum numbers.',
     &                 1, INPUT_ERR)
                  exit
               endif
            enddo

            iqn1 = iqn
            iqn2 = iqn1 + i_m - 1
            nqn0 = iqn2

c     ==== Check duplicate in n ====
            if (i_n > 1) then
               if (wfat_check_invec(n0_(i_n), i_n-1, 
     &              n0_(1:i_n-1)))
     &              call errquit('WFAT: Two of the specified n- '//
     &              'parabolic quantum numbers are identical. '//
     &              'Each of the n-quantum numbers must occur '//
     &              'only once.', 1, INPUT_ERR)
            endif
            n0(iqn1:iqn2) = n0_(i_n)

c     ==== Check duplicate in m ====
            if (i_m > 1) then
               if (wfat_check_invec(m0_(i_m), i_m-1, 
     &              m0_(1:i_m-1))) then
                  write(n_str, '(i4)') n0(i_n)
                  call errquit('WFAT: Two of the specified m-'//
     &                 'parabolic quantum numbers that correspond to '//
     &                 'n='//trim(adjustl(n_str))//' are identical. '//
     &                 'Each of the m-quantum numbers for any given '//
     &                 'n-quantum number must occur only once.',
     &                 1, INPUT_ERR)
               endif
            endif
            m0(iqn1:iqn2) = m0_(1:i_m)

            
            iqn = iqn + i_m
            if (iqn > nelm) then
               write(nelm_str, '(i4)')  nelm
               call errquit('WFAT: The number of the specified '//
     &              'parabolic quantum numbers exceeds limit, which '//
     &              'is equal to '//nelm_str//'.', 1, INPUT_ERR)
            endif
         else
            call errquit('WFAT: Wrong n-parabolic quantum number '//
     &           'input format. The correct format is nx where x is '//
     &           'zero or positive integer.', 1, INPUT_ERR)
         endif

      enddo

      wfat_inp_paraqn = .true.
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_ffc(ftp, ergf)
c     Input syntax:
c       field_factor <float energy> | unit | orb
      implicit none
#include "wfat_input.fh"
#include "inp.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer ftp(MAX_N_ORBID)
      double precision ergf(MAX_N_ORBID)

      character*255 field_a
      integer i
      double precision field_d

      wfat_inp_ffc = .false.
      
      ergf = 0.0d0
      do i = 1, MAX_N_ORBID
         if (inp_f(field_d)) then
            ergf(i) = field_d
            ftp(i) = 3
         else if (inp_a(field_a)) then
            ergf(i) = 0.0d0
            if (trim(adjustl(field_a)) == 'orb') then
               ftp(i) = 2
            else if (trim(adjustl(field_a)) == 'unit') then
               ftp(i) = 1
            else
               call errquit('WFAT: When a string is supplied for the '//
     &              'keyword ''field_factor'', the available options '//
     &              'are ''orb'' or ''unit''.', 82, UERR)
            endif
         else
            exit
         endif
      enddo

      wfat_inp_ffc = .true.
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_w1x(nc, maxn, w1x, n)

      implicit none
#include "inp.fh"
#include "errquit.fh"
      integer nc
      integer maxn
      character(len=nc) w1x(maxn)
      integer n

      integer i
      
c        wfat1_exclude \
c          "sid1 oid1 n1 m1" \
c          "sid2 oid2 n2 m2" \
c          ...

      wfat_inp_w1x = .false.
      w1x = ''
      n = 0
      do i = 1, maxn
         if (inp_a(w1x(i))) n = n + 1
      enddo

      if (n == 0) call errquit('WFAT: Input directive '//
     &     '''wfat1_exclude'' is given without entries.', 8, INPUT_ERR)

      wfat_inp_w1x = .true.

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_extract_w1x(nc, maxn, w1x, nch, orbid, 
     &     spinid, nqn, nid, mid, w1x_l)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer nc                    ! [input],
      integer maxn                  ! [input],
      character(len=nc) w1x(maxn)   ! [input],
      integer nch                   ! [input],
      integer orbid(nch)            ! [input],
      integer spinid(nch)           ! [input],
      integer nqn                   ! [input],
      integer nid(nqn)              ! [input],
      integer mid(nqn)              ! [input],
      logical w1x_l(nqn, nch)       ! [output],

      integer i, ir, ic, info
      

      w1x_l = .false.
      do i = 1, maxn
         if (len_trim(w1x(i)) > 0) then
            call wfat_read_channel_id(nch, orbid, spinid, nqn, nid, mid,
     &           w1x(i), ir, ic, info)
            if (info /= 0) then
               write(luout, '(4x, "info = ", i3)') info
               write(luout, '(4x, "Bad entry = ", a)')
     &              trim(adjustl(w1x(i)))
               if (info == 5) then
                  call errquit('WFAT: The above entry of input '//
     &                 'directive wfat1_exclude contains an '//
     &                 'unspecified ionized channel.', 8, INPUT_ERR)
               else if (info == 10) then
                  call errquit('WFAT: The above entry of input '//
     &                 'directive wfat1_exclude contains an '//
     &                 'unspecified parabolic quantum number.',
     &                 8, INPUT_ERR)
               else
                  call errquit('WFAT: The above entry of input '//
     &                 'wfat1_exclude does not conform to the '//
     &                 'correct format.', 8, INPUT_ERR)
               endif
            endif
            w1x_l(ir,ic) = .true.
         endif
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_mydip(nc, maxn, dip_s, dip, n)
c     DESCRIPTION:
c     Parse the user-specified unshifted moments of the ionized channels.
      implicit none
#include "inp.fh"
#include "errquit.fh"
#include "wfat_const.fh"
      integer nc
      integer maxn
      character(len=nc) dip_s(maxn)
      double precision dip(3,maxn)
      integer n

      integer i
      double precision tmp(3)
      

      wfat_inp_mydip = .false.

      if (inp_f(tmp(1))) then
c        mydipole  dipx dipy dipz
         n = maxn
         if (.not. inp_f(tmp(2))) call errquit('WFAT: Cannot read '//
     &        'the y component of manual dipole.', 8, INPUT_ERR)
         if (.not. inp_f(tmp(3))) call errquit('WFAT: Cannot read '//
     &        'the z component of manual dipole.', 8, INPUT_ERR)
         do i = 1, maxn
            dip_s(i) = 'all'
            dip(:,i) = tmp
         enddo
      else
c        mydipole \
c          "sid1 oid1"  dipx1 dipy1 dipz1 \
c          "sid2 oid2"  dipx2 dipy2 dipz2 \
c          ...
         dip_s = ''
         dip = DEF_MDIP
         n = 0
         do i = 1, maxn            
            if (inp_a(dip_s(i))) then
               if (.not. inp_f(tmp(1))) call errquit('WFAT: Cannot '//
     &        'read the x component of manual dipole.', 18, INPUT_ERR)
               if (.not. inp_f(tmp(2))) call errquit('WFAT: Cannot '//
     &        'read the y component of manual dipole.', 18, INPUT_ERR)
               if (.not. inp_f(tmp(3))) call errquit('WFAT: Cannot '//
     &        'read the z component of manual dipole.', 18, INPUT_ERR)
               n = n + 1
               dip(:,i) = tmp
            endif
         enddo
         if (n == 0) call errquit('WFAT: Input directive ''mydipole'''//
     &        'is given without entries.', 8, INPUT_ERR)
      endif

      wfat_inp_mydip = .true.

      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_myerg(nc, maxn, erg_s, erg, n)

      implicit none
#include "inp.fh"
#include "errquit.fh"
#include "wfat_const.fh"
      integer nc
      integer maxn
      character(len=nc) erg_s(maxn)
      double precision erg(maxn)
      integer n

      integer i
      double precision tmp
      

      wfat_inp_myerg = .false.

      if (inp_f(tmp)) then
c        myenergy  erg
         n = maxn
         do i = 1, maxn
            erg_s(i) = 'all'
            erg(i) = tmp
         enddo
      else
c        myenergy \
c          "sid1 oid1"  erg1 \
c          "sid2 oid2"  erg2 \
c          ...
         erg_s = ''
         erg = DEF_MERG
         n = 0
         do i = 1, maxn            
            if (inp_a(erg_s(i))) then
               if (.not. inp_f(tmp)) call errquit('WFAT: Cannot '//
     &        'read manual energy.', 18, INPUT_ERR)
               n = n + 1
               erg(i) = tmp
            endif
         enddo
         if (n == 0) call errquit('WFAT: Input directive ''myenergy'''//
     &        'is given without entries.', 8, INPUT_ERR)
      endif

      wfat_inp_myerg = .true.

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_mydip(nchar, dip_s, nch, orbid, spinid,
     &     dip, dip_l)
c     DESCRIPTION:
c     Rearranges contents of dip so that they are filled based on the
c     channel indices contained in dip_s.
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "wfat_const.fh"
      integer nchar                  ! [input],
      character(len=nchar) dip_s(*)  ! [input],
      integer nch                    ! [input],
      integer orbid(nch)             ! [input],
      integer spinid(nch)            ! [input],
      double precision dip(3,*)      ! [in/output],
      logical dip_l(*)               ! [output],

      integer i, ic, info
      double precision tmp(3,nch)


      tmp = DEF_MDIP
      dip_l(1:nch) = .false.
      do i = 1, nch
         if (len_trim(dip_s(i)) > 0) then
            if (trim(adjustl(dip_s(i))) == 'all') then
               dip_l(1:nch) = .true.
               tmp(:,1:nch) = dip(:,1:nch)
               exit
            else
               call wfat_read_orb_id(nch, orbid, spinid, dip_s(i), ic,
     &              info)
               if (info /= 0) then
                  write(luout, '(4x, "info = ", i3)') info
                  write(luout, '(4x, "Bad entry = ", a)')
     &                 trim(adjustl(dip_s(i)))
                  if (info == 5) then
                     call errquit('WFAT: The above entry of input '//
     &                    'directive mydipole contains an '//
     &                    'unspecified ionized channel.', 8, INPUT_ERR)
                  else
                     call errquit('WFAT: The above entry of input '//
     &                    'mydipole does not conform to the '//
     &                    'correct format.', 8, INPUT_ERR)
                  endif
               endif
               dip_l(ic) = .true.
               tmp(:,ic) = dip(:,i)
            endif
         endif
      enddo

      dip(:,1:nch) = tmp(:,1:nch)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_myerg(nchar, erg_s, nch, orbid, spinid,
     &     erg, erg_l)
c     DESCRIPTION:
c     Rearranges contents of erg so that they are filled based on the
c     channel indices contained in erg_s.
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "wfat_const.fh"
      integer nchar                  ! [input],
      character(len=nchar) erg_s(*)  ! [input],
      integer nch                    ! [input],
      integer orbid(nch)             ! [input],
      integer spinid(nch)            ! [input],
      double precision erg(*)        ! [in/output],
      logical erg_l(*)               ! [output],

      integer i, ic, info
      double precision tmp(nch)


      tmp = DEF_MERG
      erg_l(1:nch) = .false.
      do i = 1, nch
         if (len_trim(erg_s(i)) > 0) then
            if (trim(adjustl(erg_s(i))) == 'all') then
               erg_l(1:nch) = .true.
               tmp(1:nch) = erg(1:nch)
               exit
            else
               call wfat_read_orb_id(nch, orbid, spinid, erg_s(i), ic,
     &              info)
               if (info /= 0) then
                  write(luout, '(4x, "info = ", i3)') info
                  write(luout, '(4x, "Bad entry = ", a)')
     &                 trim(adjustl(erg_s(i)))
                  if (info == 5) then
                     call errquit('WFAT: The above entry of input '//
     &                    'directive myenergy contains an '//
     &                    'unspecified ionized channel.', 8, INPUT_ERR)
                  else
                     call errquit('WFAT: The above entry of input '//
     &                    'myenergy does not conform to the '//
     &                    'correct format.', 8, INPUT_ERR)
                  endif
               endif
               erg_l(ic) = .true.
               tmp(ic) = erg(i)
            endif
         endif
      enddo

      erg(1:nch) = tmp(1:nch)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_inp_mypos(rtdb, maxn)
c     DESCRIPTION:
c     Parse the user-specified center of the ionized channels.
      implicit none
#include "inp.fh"
#include "errquit.fh"
#include "wfat_const.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      integer maxn

      character*1 u
      integer i, l_s, k_s, l_p, k_p, n, i1, i2, k
      double precision tmp(3), b2a, ufactor
      parameter (b2a=0.52917721067d0)

      character(len=maxidlen) field
      

c     ==== Read unit ====
      if (inp_a(field)) then
         if (inp_compare(.false., 'angstrom', field)) then
            u = 'a'
         else if (inp_compare(.false., 'bohr', field)) then
            u = 'b'
         else
            u = 'a'
            call inp_prev_field()
         endif
      else
         u = 'a'   ! default unit is angstrom
      endif

      
c     ==== Unit's conversion factor ====
      ufactor = 1.0d0
      if (u == 'a') ufactor = 1.0d0/b2a
      

      if (.not.ma_push_get(mt_byte, maxidlen*maxn, 'chid', l_s, k_s))
     &     call errquit('WFAT: cannot allocate chid',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl, 3*maxn, 'pos', l_p, k_p))
     &     call errquit('WFAT: cannot allocate pos', 0, MA_ERR)

c     mychpos \
c       "sid1 oid1"  posx1 posy1 posz1 \
c       "sid2 oid2"  posx2 posy2 posz2 \
c       ...

      
c     ==== Parse channel positions ====
      byte_mb(k_s:k_s+maxidlen*maxn-1) = ''
      dbl_mb(k_p:k_p+3*maxn-1) = DEF_MPOS
      n = 0
      do i = 1, maxn            
         if (inp_a(field)) then
            if (.not. inp_f(tmp(1))) call errquit('WFAT: Cannot '//
     &     'read the x component of manual chpos.', 18, INPUT_ERR)
            if (.not. inp_f(tmp(2))) call errquit('WFAT: Cannot '//
     &     'read the y component of manual chpos.', 18, INPUT_ERR)
            if (.not. inp_f(tmp(3))) call errquit('WFAT: Cannot '//
     &     'read the z component of manual chpos.', 18, INPUT_ERR)
            n = n + 1
            i1 = k_s + (i-1)*maxidlen
            call wfat_mabyte_copy(maxidlen, field, i1)
            i1 = k_p + (i-1)*3 
            i2 = i1 + 3 - 1
            dbl_mb(i1:i2) = tmp * ufactor
         endif
      enddo
      if (n == 0) call errquit('WFAT: Input directive ''mychpos'''//
     &     'is given without entries.', 8, INPUT_ERR)

      
c     ==== Write to RTDB ====
      if (.not. rtdb_put(rtdb, 'wfat:inp:chpos_n', mt_int, 1, n))
     &     call errquit('WFAT: Failed writing wfat:inp:chpos_n to rtdb'
     &     , 1, INPUT_ERR)
      if (.not. rtdb_cput(rtdb, 'wfat:inp:chpos_id', n*maxidlen,
     &     byte_mb(k_s)))
     &     call errquit('WFAT: Failed writing wfat:inp:chpos_id to rtdb'
     &     , 1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:chpos', mt_dbl, 3*n,
     &     dbl_mb(k_p)))
     &     call errquit('WFAT: Failed writing wfat:inp:chpos to rtdb',
     &     1, INPUT_ERR)

      
      if (.not. MA_pop_stack(l_p))
     &     call errquit('WFAT: cannot deallocate pos', 0, MA_ERR)
      if (.not. MA_pop_stack(l_s))
     &     call errquit('WFAT: cannot deallocate chid', 0, MA_ERR)
      
      end
c=======================================================
      
      
c=======================================================
c=======================================================
      logical function wfat_inp_omgpos(au, omgpos)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "inp.fh"      
      logical au
      double precision omgpos(3)

      character*4 field
      logical l
      integer i


      wfat_inp_omgpos = .false.


      if (inp_f(omgpos(1))) then
         au = .false.
         do i = 2, 3
            if (.not. inp_f(omgpos(i))) then
               call errquit('WFAT: Input directive ''omgpos'' '//
     &           'expects either a string ''auto'' or a sequence of '//
     &           'three decimal numbers.', 1, INPUT_ERR)
            endif
         enddo
      else
         l = inp_a(field)
         if (inp_compare(.false., 'auto', field)) then
            au = .true.
            omgpos = 0.0d0
         else            
            call errquit('WFAT: Input directive ''omgpos'' expects '//
     &           'either a string ''auto'' or a sequence of three '//
     &           'decimal numbers.', 1, INPUT_ERR)
         endif
      endif
      

      wfat_inp_omgpos = .true.

      end
c=======================================================
      

c=======================================================
c=======================================================
      logical function wfat_inp_s2a1e(nq, q, d, mag, xzr, p, pparam, 
     &     nel, fedvr, nfe, fen, ndvr, nr, kein, kemult, maxit, nt)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "inp.fh"
#include "wfat_input.fh"
      integer nq
      double precision q(2)
      double precision d
      integer mag
      integer xzr
      character*3 p   ! 'cou' or 'eff'
      double precision pparam(2)
      integer nel
      logical fedvr
      integer nfe
      double precision fen(MAX_FEN)
      integer ndvr(MAX_FEN)
      integer nr(3)
      double precision kein
      double precision kemult
      integer maxit
      integer nt(3)

      character*50 field
      character*1 adum, xzr0
      logical sys_exist, sys_stat, mag_exist, mag_stat, p_exist, p_stat
      logical pparam_exist(2), pparam_stat(2), fedvr_exist, fedvr_stat
      logical fen_exist, fen_stat, ndvr_exist, ndvr_stat, lagdvr_exist
      logical lagdvr_stat, nr_exist(3), nr_stat(3), kein_exist
      logical kein_stat, kemult_exist, kemult_stat, maxit_exist
      logical maxit_stat, nt_exist(3), nt_stat(3), nel_exist, nel_stat
      logical xzr_exist, xzr_stat, using_fedvr, using_lagdvr
      integer i, n_fen, n_ndvr

      wfat_inp_s2a1e = .false.

      sys_exist = .false.
      sys_stat = .false.
      mag_exist = .false.
      mag_stat = .false.
      p_exist = .false.
      p_stat = .false.
      pparam_exist = .false.
      pparam_stat = .false.
      nel_exist = .false.
      nel_stat = .false.
      fedvr_exist = .false.
      fedvr_stat = .false.
      fen_exist = .false.
      fen_stat = .false.
      ndvr_exist = .false.
      ndvr_stat = .false.
      lagdvr_exist = .false.
      lagdvr_stat = .false.
      nr_exist = .false.
      nr_stat = .false.
      kein_exist = .false.
      kein_stat = .false.
      kemult_exist = .false.
      kemult_stat = .false.
      maxit_exist = .false.
      maxit_stat = .false.
      nt_exist = .false.
      nt_stat = .false.
      

      do while (inp_read())
         if (inp_a(field)) then

            if (inp_compare(.false., 'system', field)) then
               sys_exist = .true.
               q = -999.0d0
c              == q1 ==
               sys_stat = inp_f(q(1))
c              == internuc. distance ==
               if (inp_f(d)) then
c                 == q2 ==
                  if (.not. inp_f(q(2))) sys_stat = .false.
                  nq = 2
               else
                  nq = 1
               endif
c           ==== Magnetic quantum number ====
            else if (inp_compare(.false., 'magnetic', field)) then
               mag_exist = .true.
               mag_stat = inp_i(mag)
c           ==== xz-plane reflection symmetry ====
            else if (inp_compare(.false., 'xzreflect', field)) then
               xzr_exist = .true.
               xzr_stat = inp_a(xzr0)
c           ==== Potential type ====
            else if (inp_compare(.false., 'potential', field)) then
               p_exist = .true.
               p_stat = inp_a(p)
               if (inp_compare(.false., 'cou', p)) then
c                 Syntax: potential cou [<float soft_coulomb>]
                  if (inp_f(pparam(1))) then
                     pparam_exist(1) = .true.
                     pparam_stat(1) = .true.
                  endif
               else if (inp_compare(.false., 'eff', p)) then
c                 Syntax: potential eff <float param1> <float param2> [<integer nelec>]
                  pparam_stat(1) = inp_f(pparam(1))
                  pparam_stat(2) = inp_f(pparam(2))
                  if (pparam_stat(1) .and. pparam_stat(2)) then
                     pparam_exist(1:2) = .true.
                  endif
c                 == Number of electrons ==
                  nel_stat = inp_i(nel)           ! 3)
                  if (nel_stat) nel_exist = .true.
c                 3) By definition, when effective potential is chosen, the number
c                    of electrons can be any number N because the remaining N-1
c                    electrons are represented as an effective potential screening.
               else
                  call errquit('WFAT: The potential type of the '//
     &                 'input directive ''potential'' can only admit '//
     &                 '''cou'' (Coulomb) or ''eff'' (effective '//
     &                 'Coulomb)', 27, INPUT_ERR)
               endif
c           ==== FEDVR option ====
            else if (inp_compare(.false., 'fedvr', field)) then
               fedvr_exist = .true.
               fedvr_stat = .true.
c           ==== Finite element (FE) nodes ====
            else if (inp_compare(.false., 'fenodes', field)) then
               fen_exist = .true.
               fen = -999.0d0
               n_fen = 0
               do i = 1, MAX_FEN
                  if (inp_f(fen(i))) then
                     n_fen = n_fen + 1
                  else if (inp_a(adum)) then
                     call errquit('WFAT: The values provided for the '//
     &                    'input directive ''fenodes'' must be a '//
     &                    'decimal.', 40, INPUT_ERR)
                  else
                     fen_stat = .true.
                     exit
                  endif
               enddo
c           ==== Number of DVRs in FEs ====
            else if (inp_compare(.false., 'ndvr', field)) then
               ndvr_exist = .true.
               ndvr = -999
               n_ndvr = 0
               do i = 1, MAX_FEN
                  if (inp_i(ndvr(i))) then
                     n_ndvr = n_ndvr + 1
                  else if (inp_a(adum)) then
                     call errquit('WFAT: The values provided for the '//
     &                    'input directive ''ndvr'' must be an '//
     &                    'integer.', 40, INPUT_ERR)
                  else
                     ndvr_stat = .true.
                     exit
                  endif
               enddo
c           ==== Laguerre DVR option ====
            else if (inp_compare(.false., 'lagdvr', field)) then
               lagdvr_exist = .true.
               lagdvr_stat = .true.
c           ==== Number of Laguerre DVRs ====
            else if (inp_compare(.false., 'lagdvr_n', field)) then
               nr = -999
               nr_stat(1) = inp_i(nr(1))
               if (nr_stat(1)) nr_exist(1) = .true.
               nr_stat(2) = inp_i(nr(2))
               if (nr_stat(2)) nr_exist(2) = .true.
               nr_stat(3) = inp_i(nr(3))
               if (nr_stat(3)) nr_exist(3) = .true.
c           ==== Initial ke ====
            else if (inp_compare(.false., 'ke_init', field)) then
               kein_exist = .true.
               kein_stat = inp_f(kein)
c           ==== ke multipler ====
            else if (inp_compare(.false., 'ke_multiplier', field)) then
               kemult_exist = .true.
               kemult_stat = inp_f(kemult)
c           ==== Max. iteration to determine ke ====
            else if (inp_compare(.false., 'maxit', field)) then
               maxit_exist = .true.
               maxit_stat = inp_i(maxit)
c           ==== Number of assoc. Legendre DVRs ====
            else if (inp_compare(.false., 'legdvr_n', field)) then
               nt = -999
               nt_stat(1) = inp_i(nt(1))
               if (nt_stat(1)) nt_exist(1) = .true.
               nt_stat(2) = inp_i(nt(2))
               if (nt_stat(2)) nt_exist(2) = .true.
               nt_stat(3) = inp_i(nt(3))
               if (nt_stat(3)) nt_exist(3) = .true.
c           ==== End statement ====
            else if (inp_compare(.false., 'end', field)) then
               exit
            else
               call errquit('WFAT: Unknown input directive '''//
     &              trim(adjustl(field))//''' is found under the '//
     &              '''system_s2a1e'' subdirective', 1, INPUT_ERR)
            endif
         endif
      enddo


      using_fedvr = (fedvr_exist .or. fen_exist .or. ndvr_exist)
      using_lagdvr = (lagdvr_exist .or. nr_exist(1) .or. kein_exist
     &     .or. kemult_exist .or. maxit_exist)
      if (using_fedvr .and. using_lagdvr)
     &     call errquit('WFAT: The basis choices FEDVR and '//
     &     'Laguerre DVR are mutually exclusive. Meaning, the '//
     &     'appearance of any FEDVR-related directives (''fedvr'', '//
     &     '''nodes'', and ''ndvr'') implies the mandatory absence '//
     &     'of any Laguerre DVR-related directives (''lagdvr'', '//
     &     '''lagdvr_nr'', ''ke_init'', ''ke_multplier'', and '//
     &     '''maxit''), and vice versa.)', 50, INPUT_ERR)
            
      
c     ==== Mandatory inputs' existence checks ====
      if (.not. sys_exist)
     &     call errquit('INPUT ERROR: Cannot find input directive '//
     &     '''system'' needed for WFAT-2a1e calculation.', 0, INPUT_ERR)
      if (.not. mag_exist)
     &     call errquit('INPUT ERROR: Cannot find input directive '//
     &     '''magnetic'' needed for WFAT-2a1e calculation.',
     &     0, INPUT_ERR)
      if (.not. xzr_exist)
     &     call errquit('INPUT ERROR: Cannot find input directive '//
     &     '''xzreflect'' needed for WFAT-2a1e calculation.',
     &     0, INPUT_ERR)

c     ==== Inputs values checks ====
c     ==== xz-plane reflection symmetry ====
      if (xzr0 == 's') then
         xzr = 1
      else if (xzr0 == 'a') then
         xzr = 2
         if (mag == 0) call errquit('WFAT: States with a magnetic '//
     &        'quantum number of 0 can only be symmetric under '//
     &        'reflection onto xz-plane.', 14, INPUT_ERR)
      else
         call errquit('WFAT: Possible values for the input directive '//
     &        '''xzreflect'' are ''s'' (symmetric) or ''a'' (anti '//
     &        'symmetric).', 15, INPUT_ERR)
      endif
c     ==== Basis ====
      if (.not. using_fedvr .and. .not. using_lagdvr)       ! 1)
     &     using_fedvr = .true.
c     1) If no basis information is provided, then defaults to FEDVR.
      
      fedvr = using_fedvr       !1)
c     1) The choice of radial basis can be implied from the existence of
c        directives related to the basis. However, it is recommended to
c        explicitly state the choice by typing 'fedvr' for FEDVR or 'lagdvr'
c        for Laguerre DVR.

      if (fedvr) then
c        == FE nodes ==
         if (.not. fen_exist) then
            n_fen = DEF_NFENODES
            fen(1:n_fen) = (/ 1.0, 3.0, 6.0, 10.0, 14.0, 18.0, 22.0,
     &                        26.0, 30.0 /)
         endif
c        == FE number of DVRs ==
         if (.not. ndvr_exist) then
            n_ndvr = DEF_NFENODES
            ndvr(1:n_ndvr) = (/ 6, 25, 18, 16, 15, 15, 10, 10, 10 /)
         endif
c        == Number of FEs ==
         if (n_fen == n_ndvr) then
            nfe = n_fen
         else
            call errquit('WFAT: The number of entries provided for '//
     &           '''fenodes'' must be the same as that for ''ndvr''.',
     &           21, INPUT_ERR)
         endif
      else
c        == Number of radial Laguerre DVRs ==
         if (.not. nr_exist(1)) nr(1) = DEF_LAGDVR_NR
         if (.not. nr_exist(2)) nr(2) = nr(1)
         if (.not. nr_exist(3)) nr(3) = nr(1)
c        == Initial/guess ke ==
         if (.not. kein_exist) kein = DEF_KEIN
c        == ke multiplier ==
         if (.not. kemult_exist) kemult = DEF_KEMULT
c        == Max. iterations ==
         if (.not. maxit_exist) maxit = DEF_2A1E_MAXIT
      endif

c     == Number of angular assoc. Legendre DVRs ==
      if (.not. nt_exist(1)) nt(1) = DEF_LEGDVR_NT
      if (.not. nt_exist(2)) nt(2) = nt(1)
      if (.not. nt_exist(3)) nt(3) = nt(1)      
      
c     ==== Potential ====
      if (.not. p_exist) then
         if (nq == 1) then
            p = 'eff'
         else if (nq == 2) then
            p = 'cou'
         endif
      endif
      if (p == 'cou') then
         if (.not. pparam_exist(1)) pparam(1) = DEF_SOFT_COULOMB
      else if (p == 'eff') then
         if (.not. (pparam_exist(1) .and. pparam_exist(2)))
     &        call errquit('WFAT: When potential type is ''eff'', two'//
     &        ' potential parameters are required.', 14, INPUT_ERR)
      endif
      
c     ==== Number of electrons ====
      if (.not. nel_exist) then         !2)
         if (p == 'cou') then
            nel = 1
         else if (p == 'eff') then
            nel = nint( sum(q(1:nq)) )
         endif
      endif
c     2) The default value for the number of electrons depends on potential
c        type.

      wfat_inp_s2a1e = .true.
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_perturb(
     &                            type,   type_exist,   type_stat,
     &                            ord,    ord_exist,    ord_stat,
     &                            bt,     bt_exist,     bt_stat,
     &                            gm,     gm_exist,     gm_stat,
     &                            dgnthr, dgnthr_exist, dgnthr_stat,
     &                            offthr, offthr_exist, offthr_stat, 
     &                            e1thr,  e1thr_exist,  e1thr_stat, 
     &                            itrthr, itrthr_exist, itrthr_stat,  
     &                            maxitr, maxitr_exist, maxitr_stat,  
     &                            ordthr, ordthr_exist, ordthr_stat,
     &                            maxord, maxord_exist, maxord_stat)

      implicit none
#include "inp.fh"
#include "errquit.fh"      
      character*2 type                   ! [output]
      logical type_exist, type_stat      ! [output]
      integer ord                        ! [output]
      logical ord_exist, ord_stat        ! [output]
      double precision bt                ! [output]
      logical bt_exist, bt_stat          ! [output]
      double precision gm                ! [output]
      logical gm_exist, gm_stat          ! [output]
      double precision dgnthr            ! [output]
      logical dgnthr_exist, dgnthr_stat  ! [output]
      double precision offthr            ! [output]
      logical offthr_exist, offthr_stat  ! [output]
      double precision e1thr             ! [output]
      logical e1thr_exist, e1thr_stat    ! [output]
      double precision itrthr            ! [output]
      logical itrthr_exist, itrthr_stat  ! [output]
      integer maxitr                     ! [output]
      logical maxitr_exist, maxitr_stat  ! [output]
      double precision ordthr            ! [output]
      logical ordthr_exist, ordthr_stat  ! [output]
      integer maxord                     ! [output]
      logical maxord_exist, maxord_stat  ! [output]

      
      character*50 field
      integer nfields, i
      

      type_exist = .false.
      type_stat = .false.

      ord_exist = .false.
      ord_stat = .false.

      bt_exist = .false.
      bt_stat = .false.

      gm_exist = .false.
      gm_stat = .false.

      dgnthr_exist = .false.
      dgnthr_stat = .false. 

      offthr_exist = .false.
      offthr_stat = .false. 

      e1thr_exist = .false. 
      e1thr_stat = .false.

      itrthr_exist = .false.
      itrthr_stat = .false. 

      maxitr_exist = .false.
      maxitr_stat = .false. 

      ordthr_exist = .false.
      ordthr_stat = .false. 

      maxord_exist = .false.
      maxord_stat = .false. 


      nfields = inp_n_field()
      do i = 2, nfields
         if (inp_a(field)) then

            if (inp_compare(.false., 'bw', field)) then
               type = 'bw'
               type_exist = .true.
               type_stat = .true.
            else if (inp_compare(.false., 'rs', field)) then
               type = 'rs'
               type_exist = .true.
               type_stat = .true.
            else if (inp_compare(.false., 'order', field)) then
               ord_exist = .true.
               ord_stat = inp_i(ord)
            else if (inp_compare(.false., 'beta', field)) then
               bt_exist = .true.
               bt_stat = inp_f(bt)
            else if (inp_compare(.false., 'gamma', field)) then
               gm_exist = .true.
               gm_stat = inp_f(gm)

            else if (inp_compare(.false., 'dgnthr', field)) then
               dgnthr_exist = .true.
               dgnthr_stat = inp_f(dgnthr)
            else if (inp_compare(.false., 'offthr', field)) then
               offthr_exist = .true.
               offthr_stat = inp_f(offthr)
            else if (inp_compare(.false., 'e1thr', field)) then
               e1thr_exist = .true.
               e1thr_stat = inp_f(e1thr)
            else if (inp_compare(.false., 'itrthr', field)) then
               itrthr_exist = .true.
               itrthr_stat = inp_f(itrthr)
            else if (inp_compare(.false., 'maxitr', field)) then
               maxitr_exist = .true.
               maxitr_stat = inp_i(maxitr)
            else if (inp_compare(.false., 'ordthr', field)) then
               ordthr_exist = .true.
               ordthr_stat = inp_f(ordthr)
            else if (inp_compare(.false., 'maxord', field)) then
               maxord_exist = .true.
               maxord_stat = inp_i(maxord)
            
            else
               call errquit('WFAT: The specified keyword '''
     &              //trim(adjustl(field))//''' inside ''perturb'' '//
     &              'directive is unavailable.', 1, INPUT_ERR)
            endif
cOLD         else
cOLD            call errquit('WFAT: One parameter for directive '//
cOLD     &           '''perturb'' is not preceded by a string.', 
cOLD     &           1, INPUT_ERR)
         endif
      enddo


      if (type /= 'bw' .and. type /= 'rs') then        ! Just for safety.
         type_stat = .false.
      endif

      
      wfat_inp_perturb = .true.

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_parse_quadgrid(rtdb)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb

      character*50 field
      logical stat
      integer ozn, gs
      double precision cgrid_thr


      do while (inp_read())
         if (inp_a(field)) then
c           ==== gridset ====
            if (inp_compare(.false., 'gridset', field)) then
               stat = inp_i(gs)
               if (stat) then
                  if (.not. rtdb_put(rtdb, 'wfat:inp:gqual', mt_int, 1,
     &                 gs)) call errquit(
     &                 'WFAT: Failed writing inp:gqual to RTDB', 1,
     &                 RTDB_ERR)
               endif
c           ==== omegaz ====
            else if (inp_compare(.false., 'omegaz', field)) then
               stat = inp_i(ozn)
               if (stat) then
                  if (.not. rtdb_put(rtdb, 'wfat:inp:ozn', mt_int, 1,
     &                 ozn)) call errquit(
     &                 'WFAT: Failed writing inp:ozn to RTDB',
     &                 1, RTDB_ERR)
               endif
c           ==== cgridthr ====
            else if (inp_compare(.false., 'cgridthr', field)) then
               stat = inp_f(cgrid_thr)
               if (stat) then
                  if (.not. rtdb_put(rtdb, 'wfat:inp:cgrid_thr', mt_dbl,
     &                 1, cgrid_thr)) call errquit(
     &                 'WFAT: Failed writing inp:cgrid_thr to RTDB',
     &                 1, RTDB_ERR)
               endif
            else if (inp_compare(.false., 'end', field)) then
               exit
            else
               call errquit('WFAT: Unknown input directive '''//
     &              trim(adjustl(field))//''' is found under the '//
     &              '''quadgrid'' subdirective', 1, INPUT_ERR)
            endif

         endif
      enddo

      end
c=======================================================
      
      
cc=======================================================
cc=======================================================
c      subroutine wfat_process_input()
c
c      implicit none
c
c
cc     ==== beta ====
c      if (wfat_nbeta > 1) then
c         dbt = (wfat_beta_last-wfat_beta_start) / (wfat_nbeta-1)
c         do i = 1, wfat_nbeta
c            wfat_beta(i) = wfat_beta_start + (i-1)*dbt
c         enddo
c      else
c         wfat_beta(1) = wfat_beta_start
c      endif
c
cc     ==== gamma ====
c      if (wfat_ngamma > 1) then
c         dgm = (wfat_gamma_last-wfat_gamma_start) / (wfat_ngamma-1)
c         do i = 1, wfat_ngamma
c            wfat_gamma(i) = wfat_gamma_start + (i-1)*dgm
c         enddo
c      else
c         wfat_gamma(1) = wfat_gamma_start
c      endif
c
cc     ==== orbid ====
c      
c      
c      
c
c      end
cc=======================================================
     
