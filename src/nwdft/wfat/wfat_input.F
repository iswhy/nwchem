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
      write(luout, '(5x, "Reading WFAT inputs ...")')
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
#include "mafdecls.fh"
#include "rtdb.fh"      
      integer rtdb
      
      logical tt, ll1, ll2, ll3
      character*50 field, tail_unit
      integer i1, i2, k, nfblocks
      double precision ang2bohr

      parameter(ang2bohr=1.8897261255d0)
      
      logical wfat_inp_omgpos
      logical wfat_inp_perturb
      logical wfat_inp_s2a1e
      external wfat_inp_omgpos
      external wfat_inp_perturb
      external wfat_inp_s2a1e

      
c     ==== Initial status of WFAT inputs ====
      nfblocks = 0
      
      nktmd_exist = .false.
      nktmd_stat = .false.

      omgpos_exist = .false.
      omgpos_stat = .false.

      movec_exist = .false.
      movec_stat = .false.
      movec1_exist = .false.
      movec1_stat = .false.

      wfat1_exist = .false.
      wfat1_stat = .false.

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

      s2a1e_exist = .false.
      s2a1e_stat = .true.

      
c     ========================================================
c     ==================== INPUTS PARSING ====================
c     ========================================================

      do while (inp_read()) 
         if (inp_a(field)) then

cfordebug            write(luout, '("input field : ", a)') field

            if (inp_compare(.false., 'print', field)) then
               call wfat_parse_print(rtdb)
            else if (inp_compare(.false., 'field', field)) then
               nfblocks = nfblocks + 1
               call wfat_parse_field(rtdb, nfblocks)
            else if (inp_compare(.false., 'dt', field)) then
               call wfat_parse_dt(rtdb)
            else if (inp_compare(.false., 'tmax', field)) then
               call wfat_parse_tmax(rtdb)               
            else if (inp_compare(.false., 'beta', field)) then
               call wfat_parse_angle(rtdb, 'beta')
            else if (inp_compare(.false., 'gamma', field)) then
               call wfat_parse_angle(rtdb, 'gamma')
            else if (inp_compare(.false., 'kx', field)) then
               call wfat_parse_momentum(rtdb, 'x')
            else if (inp_compare(.false., 'ky', field)) then
               call wfat_parse_momentum(rtdb, 'y')
            else if (inp_compare(.false., 'movec_id', field)) then
               call wfat_parse_movec(rtdb)
            else if(inp_compare(.false., 'paraqn', field)) then
               call wfat_parse_paraqn(rtdb)
            else if (inp_compare(.false., 'lmax', field)) then
               call wfat_parse_inpi(rtdb, 'lmax')
            else if (inp_compare(.false., 'quadgrid', field)) then
               call wfat_parse_quadgrid(rtdb)
            else if (inp_compare(.false., 'si_type', field)) then
               call wfat_parse_inpa(rtdb, 'si_type')
            else if (inp_compare(.false., 'degenthr', field)) then
               call wfat_parse_inpf(rtdb, 'dgnthr')
            else if (inp_compare(.false., 'intgl', field)) then
               call wfat_parse_intgl(rtdb)
            else if (inp_compare(.false., 'vectors', field)) then   ! Needed for restarting MEWFAT task.
c              Loads MO coefficients from a file. Along with the restart directive, it can be
c              used to skip the SCF calculation that yield the orbitals.
c                Syntax:  vectors <movec file1> [movec file2]
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
               call wfat_parse_w1x(rtdb)
            else if (inp_compare(.false., 'mydipole', field)) then
               call wfat_parse_mydip(rtdb)
            else if (inp_compare(.false., 'myenergy', field)) then
               call wfat_parse_myerg(rtdb)
            else if (inp_compare(.false., 'mychpos', field)) then
               call wfat_inp_mypos(rtdb, MAX_N_ORBID)
            else if (inp_compare(.false., 'myffactor', field)) then
               call wfat_parse_ffc(rtdb)
            else if (inp_compare(.false., 'nlgl', field)) then
               call wfat_parse_inpi(rtdb, 'nlgl')
            else if (inp_compare(.false., 'diagthr', field)) then
               call wfat_parse_inpf(rtdb, 'dthr')
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
            else if (inp_compare(.false., 'system_2a1e', field)) then
               s2a1e_exist = .true.
               s2a1e_stat = wfat_inp_s2a1e(wfat2a_nc, wfat2a_q, 
     &              wfat2a_inter, wfat2a_mag, wfat2a_xzr, wfat2a_pot, 
     &              wfat2a_potparam, wfat2a_nel, wfat2a_fedvr, 
     &              wfat2a_nfe, wfat2a_fenodes, wfat2a_ndvr, wfat2a_nr, 
     &              wfat2a_kein, wfat2a_kemult, wfat2a_maxit, wfat2a_nt)
               s2a1e = s2a1e_stat
               if (.not. rtdb_put(rtdb, 'wfat:inp:s2a1e', mt_log, 1, 
     &              s2a1e)) call errquit('WFAT: Failed writing '//
     &              'inp:s2a1e to rtdb', 42, RTDB_ERR)
            else if (inp_compare(.false., 'omgpos', field)) then
               omgpos_exist = .true.
               omgpos_stat = wfat_inp_omgpos(autopos, omgpos)
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
      subroutine wfat_parse_print(rtdb)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      
      character*255 ss
      character*10 options(7)
      logical tt
      integer i1, i2, k

      
      options = (/
     &     'srfc      ',
     &     'rate      ',
     &     'yield     ',
     &     'avbeta    ', 
     &     'avgamma   ',
     &     'tmd       ',
     &     'tmd_atom  ' /)


c     ==== Read the entire line containing the print directive ====
      if (.not. inp_line(ss)) call errquit('WFAT: Error in reading '//
     &     'the ''print'' line.', 1, INPUT_ERR)
      i1 = 0       ! i1 must be initiated to zero before first call of inp_strtok
      tt = inp_strtok(ss, ' ', i1, i2)     !! The 'print' keyword

      
c     ==== Iterate over keywords ====
      do while (inp_strtok(ss, ' ', i1, i2))        !! The rest of 'print' line
         tt = inp_match(7, .false., ss(i1:i2), options, k)

         if (options(k) == 'srfc') then
            if (.not. rtdb_put(rtdb, 'wfat:inp:p_srfc', mt_log, 1,
     &           .true.)) call errquit('WFAT: Failed writing '//
     &           'inp:p_srfc to rtdb', 1, RTDB_ERR)

         else if (options(k) == 'rate') then
            if (.not. rtdb_put(rtdb, 'wfat:inp:p_rate', mt_log, 1,
     &           .true.)) call errquit('WFAT: Failed writing '//
     &           'inp:p_rate to rtdb', 1, RTDB_ERR)

         else if (options(k) == 'yield') then
            if (.not. rtdb_put(rtdb, 'wfat:inp:p_yild', mt_log, 1,
     &           .true.)) call errquit('WFAT: Failed writing '//
     &           'inp:p_yild to rtdb', 1, RTDB_ERR)

         else if (options(k) == 'avbeta') then
            if (.not. rtdb_put(rtdb, 'wfat:inp:p_avbeta', mt_log, 1,
     &           .true.)) call errquit('WFAT: Failed writing '//
     &           'inp:p_avbeta to rtdb', 1, RTDB_ERR)

         else if (options(k) == 'avgamma') then
            if (.not. rtdb_put(rtdb, 'wfat:inp:p_avgamma', mt_log, 1,
     &           .true.)) call errquit('WFAT: Failed writing '//
     &           'inp:p_avgamma to rtdb', 1, RTDB_ERR)

         else if (options(k) == 'tmd') then
            if (.not. rtdb_put(rtdb, 'wfat:inp:p_tmd', mt_log, 1,
     &           .true.)) call errquit('WFAT: Failed writing '//
     &           'inp:p_tmd to rtdb', 1, RTDB_ERR)

         else if (options(k) == 'tmd_atom') then
            if (.not. rtdb_put(rtdb, 'wfat:inp:p_tmd_atom', mt_log, 1,
     &           .true.)) call errquit('WFAT: Failed writing '//
     &           'inp:p_tmd_atom to rtdb', 1, RTDB_ERR)

         else
            call errquit('WFAT: Unknown keyword is specified for the '//
     &           'print directive', 64, INPUT_ERR)
cOLD            exit
         endif
      enddo
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_field(rtdb, i)

c     field
c        type < dc | gaussian | hann | external <"laser_file".laser> >
c        frequency <double> [ au | ev | thz | nm | bohr ]
c        center <double> [ fs | au ]
c        width <double> [ fs | au ]
c        max <double> [ au | tw/cm2 ]
c        cep <double> [ rad | deg ]
c     end
      implicit none
#include "wfat_field_param.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      integer i                 ! [input],

      character*50 key, tp0, fru, ceu, wiu, mxu, cpu
      character*2 max_str
      character*2 is
      character*4 notype
      character*2 tp, TP_OPTS(4)
      logical l
      integer k, typeg
      double precision mx, mx0, fr, fr0, ce, ce0, wi, wi0, cp, cp0

      double precision wfat_convert_freq, wfat_convert_time,
     &     wfat_convert_efield, wfat_convert_angle
      external wfat_convert_freq, wfat_convert_time,
     &     wfat_convert_efield, wfat_convert_angle

      parameter (notype='none')


      TP_OPTS = (/ 'dc', 'ga', 'ha', 'ex' /)
c     dc = DC/static field
c     ga = Gaussian envelope
c     ha = Hann envelope
c     ex = Field from an external file
      

c     ==== Ensure the number of field blocks is below limit ====
      write(max_str, '(i2)') max_nfield
      if (i <= max_nfield) then
         if (.not. rtdb_put(rtdb, 'wfat:field:n', mt_int, 1, i)) 
     &        call errquit('WFAT: Failed writing field:n to rtdb',
     &        1, RTDB_ERR)
      else
         call errquit('WFAT: No. of fields exceeds the maximum, '//
     &        'which is '//trim(adjustl(max_str)), 1, INPUT_ERR)
      endif

      
c     ==== Parse field parameters ====
      tp0 = notype
      do while (inp_read())
         if (inp_a(key)) then

c           >>>> Field type (dc, external, gaussian, etc) <<<<
            if (inp_compare(.false., 'type', key)) then
               l = inp_a(tp0)
               call wfat_abbrev_fieldtp(tp0, tp)
               if (tp == 'ex') call wfat_inp_externfield(rtdb)

c           >>>> Frequency <<<<
            else if (inp_compare(.false., 'frequency', key)) then
               l = inp_f(fr0)
               if (.not. inp_a(fru)) fru = 'au'
               
c           >>>> Center <<<<
            else if (inp_compare(.false., 'center', key)) then
               l = inp_f(ce0)
               if (.not. inp_a(ceu)) ceu = 'au'

c           >>>> Width <<<<
            else if (inp_compare(.false., 'width', key)) then
               l = inp_f(wi0)
               if (.not. inp_a(wiu)) wiu = 'au'
               
c           >>>> Amplitude <<<<
            else if (inp_compare(.false., 'max', key)) then
               l = inp_f(mx0)
               if (.not. inp_a(mxu)) mxu = 'au'
               
c           >>>> CEP <<<<
            else if (inp_compare(.false., 'cep', key)) then
               l = inp_f(cp0)
               if (.not. inp_a(cpu)) cpu = 'rad'

            else if (inp_compare(.false., 'end', key)) then
               exit
            else
               call errquit('WFAT: Unknown input directive '''//
     &              trim(adjustl(key))//''' is found under the '//
     &              '''field'' subdirective', 1, INPUT_ERR)
            endif

         endif
      enddo


      write(is, '(i2)') i
      is = adjustl(is)


c     ==== Write field parameters to RTDB ====
c     >>>> Field type <<<<
      if (inp_match(4, .false., tp, TP_OPTS, k)) then
         if (.not. rtdb_cput(rtdb, 'wfat:inp:field:type:'//is, 1, tp))
     &        call errquit('WFAT: Failed writing inp:field:type:'//
     &        is//' to rtdb', 1, RTDB_ERR)
         
c        ==== Identify type group ====
c        typeg = 1 ==> dc
c        typeg = 2 ==> external
c        typeg = 3 ==> everything else
         if (inp_compare(.false., 'dc', tp)) then
            typeg = 1
         else if (inp_compare(.false., 'ex', tp)) then
            typeg = 2
         else
            typeg = 3
         endif
         if (.not. rtdb_put(rtdb, 'wfat:field:typeg', mt_int, 1, typeg))
     &        call errquit('WFAT: Failed writing field:typeg to rtdb',
     &        1, RTDB_ERR)

c        ==== Check for the number of field blocks ====
c        ====  when dc or external field is used   ====
         if (i > 1) then
            if (typeg==1 .or. typeg==2) 
     &           call errquit('WFAT: Only one field block is allowed '//
     &           'when dc or external field is specified.',
     &           1, INPUT_ERR)
         endif      
      else
         call errquit('WFAT: The type directive with a proper value'//
     &        ' is missing in a field block.', 40, INPUT_ERR)
      endif
      
      if (typeg==1 .or. typeg==3) then
c        >>>> Amplitude <<<<
         mx = wfat_convert_efield(trim(adjustl(mxu)), mx0, 'au')
         if (.not. rtdb_put(rtdb, 'wfat:inp:field:max:'//is, mt_dbl, 1,
     &        mx)) call errquit('WFAT: Failed writing inp:field:max:'//
     &        is//' to rtdb', 1, RTDB_ERR)
      endif

      if (typeg==3) then
c        >>>> Frequency <<<<
         fr = wfat_convert_freq(trim(adjustl(fru)), fr0, 'au')
         if (.not. rtdb_put(rtdb, 'wfat:inp:field:freq:'//is, mt_dbl, 1,
     &        fr)) call errquit('WFAT: Failed writing inp:field:freq:'//
     &        is//' to rtdb', 1, RTDB_ERR)
c        >>>> Center <<<<
         ce = wfat_convert_time(trim(adjustl(ceu)), ce0, 'au')
         if (.not. rtdb_put(rtdb, 'wfat:inp:field:center:'//is, mt_dbl, 
     &        1, ce)) call errquit('WFAT: Failed writing '//
     &        'inp:field:center:'//is//' to rtdb', 1, RTDB_ERR)
c        >>>> Width <<<<
         wi = wfat_convert_time(trim(adjustl(wiu)), wi0, 'au')
         if (.not. rtdb_put(rtdb, 'wfat:inp:field:width:'//is, mt_dbl, 
     &        1, wi)) call errquit('WFAT: Failed writing '//
     &        'inp:field:width:'//is//' to rtdb', 1, RTDB_ERR)
c        >>>> CEP <<<<
         cp = wfat_convert_angle(trim(adjustl(cpu)), cp0, 'rad')
         if (.not. rtdb_put(rtdb, 'wfat:inp:field:cep:'//is, mt_dbl, 1, 
     &        cp)) call errquit('WFAT: Failed writing inp:field:cep:'//
     &        is//' to rtdb', 1, RTDB_ERR)
      endif
      
         
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_abbrev_fieldtp(tp0, tp)

      implicit none
#include "errquit.fh"
      character*(*) tp0
      character*2 tp

      if (trim(adjustl(tp0)) == 'dc') then
         tp = 'dc'
      else if (trim(adjustl(tp0)) == 'gaussian') then
         tp = 'ga'
      else if (trim(adjustl(tp0)) == 'hann') then
         tp = 'ha'
      else if (trim(adjustl(tp0)) == 'external') then
         tp = 'ex'
      else
         call errquit('WFAT: Unkown field type.', 90, INPUT_ERR)
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_expand_fieldtp(tp0, tp)

      implicit none
#include "errquit.fh"
      character*2 tp0
      character*(*) tp

      if (tp0 == 'dc') then
         tp = 'dc'
      else if (tp0 == 'ga') then
         tp = 'gaussian'
      else if (tp0 == 'ha') then
         tp = 'hann'
      else if (tp0 == 'ex') then
         tp = 'external'
      else
         call errquit('WFAT: Unkown field type.', 91, INPUT_ERR)
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_inp_externfield(rtdb)
c     <float dt>  <string t_unit>
c     <string f_unit>
c     data:
c     <float f(1)>
c     <float f(2)>
c     ...
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      
      character*2 t_unit
      character*6 f_unit
      character*2 TM_UNIT(2)
      character*6 FL_UNIT(2)
      character*5 ss
      character*1000 ext_fn, ext_fnb, dts, line
      integer k, h, bu, io, ntime
      double precision dt, fval

      
      h = 1000
      bu = h + 1
      TM_UNIT = (/ 'fs', 'au' /)
      FL_UNIT = (/ 'au    ', 'tw/cm2' /)
      
      if (inp_a(ext_fn)) then
c        ==== Open external field file ====
         open(h, file=ext_fn, action='read', status='old',
     &        iostat=io)

         if (io == 0) then
            if (.not. rtdb_cput(rtdb, 'wfat:field:ext_fn', 1,
     &           trim(adjustl(ext_fn)))) call errquit('WFAT: Failed '//
     &           'writing field:ext_fn to rtdb', 1, RTDB_ERR)
            
c           == open a binary file named <ext_fn>.dat ==
            ext_fnb = trim(adjustl(ext_fn))//'.dat'
            open(bu, file=ext_fnb, action='write', form='unformatted')
            if (.not. rtdb_cput(rtdb, 'wfat:field:ext_fnb', 1,
     &           trim(adjustl(ext_fnb)))) call errquit('WFAT: Failed '//
     &           'writing field:ext_fnb to rtdb', 1, RTDB_ERR)
         endif

c        ==== Read dt and its unit ====
         read(h, *) dt, t_unit
         if (inp_match(2, .false., t_unit, TM_UNIT, k)) then
            if (.not. rtdb_cput(rtdb, 'wfat:field:ext_tu', 1, t_unit)) 
     &        call errquit('WFAT: Failed writing field:ext_tu to rtdb',
     &        1, RTDB_ERR)
         else
            call errquit('WFAT: The unit choice for time data of the '//
     &           'external field is unavailable', 1, INPUT_ERR)
         endif

         if (.not. rtdb_put(rtdb, 'wfat:inp:field:dt', mt_dbl, 1, dt))
     &        call errquit('WFAT: Failed writing inp:field:dt to rtdb',
     &        1, RTDB_ERR)

c        ==== Read field data unit ====
         read(h, '(a)') f_unit
         if (inp_match(2, .false., f_unit, FL_UNIT, k)) then
            if (.not. rtdb_cput(rtdb, 'wfat:field:ext_fu', 1,
     &        trim(adjustl(f_unit)))) call errquit('WFAT: Failed '//
     &        'writing field:ext_fu to rtdb', 1, RTDB_ERR)
         else
            call errquit('WFAT: The unit choice for field data of '//
     &           'the external field is unavailable', 1, INPUT_ERR)
         endif
         
c        ==== Calculate ntime and tmax ====
         read(h, '(a)') ss
         if (ss == 'data:') then
            ntime = 0
            do
               read(h, '(a)', iostat=io) line
               if (io /= 0) exit  ! Exit the loop on error or end-of-file
            
               if (len( trim(adjustl(line)) ) > 0) then
                  ntime = ntime + 1
                  read(line, *) fval
                  write(bu) fval
               end if
            end do
            if (.not. rtdb_put(rtdb, 'wfat:field:ntime', mt_int, 1,
     &           ntime)) call errquit('WFAT: Failed writing '//
     &           'field:ntime to rtdb', 1, RTDB_ERR)
            if (.not. rtdb_put(rtdb, 'wfat:inp:field:tmax', mt_dbl, 1,
     &           (ntime-1)*dt)) call errquit('WFAT: Failed writing '//
     &           'inp:field:tmax to rtdb', 1, RTDB_ERR)
         else
            call errquit('WFAT: The keyword ''data:'' is not found '//
     &           'in the external field file.', 1, INPUT_ERR)
         endif
         
         close(h)
         close(bu)
      else
         call errquit('WFAT: The file containing the external field '//
     &        'is incorrectly formatted or does not exist.',
     &        1, INPUT_ERR)
      endif


      write(luout, '(/, 3x, "*** ATTENTION ***")')
      write(luout, 101)

 101  format(
     &     5x, 'The field has been taken from an external file.',
     &     'In this case, the user is ', /,
     &     5x, 'responsible in ensuring that the time data values ',
     &     'are equally spaced.')
      

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_parse_dt(rtdb)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer rtdb
      
      character*2 un
      character*2 DT_UNIT(2)
      integer k
      double precision dt, dt0
      double precision wfat_convert_time
      external wfat_convert_time

      
      DT_UNIT = (/ 'au', 'fs' /)

      if (.not. inp_f(dt0))
     &     call errquit('WFAT: Cannot parse dt', 88, INPUT_ERR)
      
      if (inp_a(un)) then
         if (.not. inp_match(2, .false., un, DT_UNIT, k)) 
     &        call errquit('WFAT: The unit choice for the directive '//
     &        '''dt'' is unavailable', 1, INPUT_ERR)
         dt = wfat_convert_time(un, dt0, 'au')
      else
         dt = dt0
      endif

      if (.not. rtdb_put(rtdb, 'wfat:inp:field:dt', mt_dbl, 1, dt))
     &     call errquit('WFAT: Cannot write inp:field:dt to RTDB.',
     &     21, RTDB_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_tmax(rtdb)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer rtdb

      character*2 un
      character*2 TMAX_UNIT(2)
      integer k
      double precision tmax, tmax0
      double precision wfat_convert_time
      external wfat_convert_time

      
      TMAX_UNIT = (/ 'au', 'fs' /)

      if (.not. inp_f(tmax0))
     &     call errquit('WFAT: Cannot parse tmax', 88, INPUT_ERR)
      
      if (inp_a(un)) then
         if (.not. inp_match(2, .false., un, TMAX_UNIT, k)) 
     &        call errquit('WFAT: The unit choice for the directive '//
     &        '''tmax'' is unavailable', 1, INPUT_ERR)
         tmax = wfat_convert_time(un, tmax0, 'au')
      else
         tmax = tmax0
      endif

      if (.not. rtdb_put(rtdb, 'wfat:inp:field:tmax', mt_dbl, 1, tmax))
     &     call errquit('WFAT: Cannot write :inp:field:tmax to RTDB.',
     &     21, RTDB_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_movec(rtdb)
c     (rtdb, nelm, orbid0, spinid0, nor0)
c     Input syntax:
c       a <list of 0 or negative integers> \
c       b <list of 0 or negative integers>
      implicit none
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "wfat_const.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer rtdb

      character*4 nelm_str
      character*1 afield
      integer i, j, ior, ior1, ior2, is, io, nor0, ifield
      integer orbid0(MAX_N_ORBID), spinid0(MAX_N_ORBID)
      integer orbid0_(MAX_N_ORBID), spinid0_(MAX_N_ORBID)

      logical wfat_check_invec
      external wfat_check_invec


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
     &              '(corresponding to occupied orbitals) are '//
     &              'supported.', 1, INPUT_ERR)
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
            if (wfat_check_invec(spinid0_(is), is-1, spinid0_(1:is-1))) 
     &           call errquit('WFAT: A particular spin channel can '//
     &           'only be specified once.', 1, INPUT_ERR)
         endif
         spinid0(ior1:ior2) = spinid0_(is)

c        ==== Check duplicate in orbital ====
         if (io > 1) then
            if (wfat_check_invec(orbid0_(io), io-1, orbid0_(1:io-1)))
     &           call errquit('WFAT: A particular spatial orbital in '//
     &           'a given spin channel can only be specified once.',
     &           1, INPUT_ERR)
         endif
         orbid0(ior1:ior2) = orbid0_(1:io)
         

         ior = ior + io
         if (ior > MAX_N_ORBID) then
            write(nelm_str, '(i4)')  MAX_N_ORBID
            call errquit('WFAT: The number of the specified '//
     &           'molecular orbitals exceeds limit, which '//
     &           'is equal to '//nelm_str//'.', 1, INPUT_ERR)
         endif

      enddo

      if (.not. rtdb_put(rtdb, 'wfat:inp:n_orbid', mt_int, 1, nor0))
     &     call errquit('WFAT: Failed writing wfat:inp:n_orbid to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:orbid', mt_int, nor0, orbid0))
     &     call errquit('WFAT: Failed writing inp:orbid to rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:spinid', mt_int, nor0, 
     &     spinid0)) call errquit('WFAT: Failed writing inp:spinid to'//
     &     ' rtdb', 1, INPUT_ERR)


      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_parse_paraqn(rtdb)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "wfat_const.fh"
      integer rtdb
      
      character*4 nelm_str, n_str
      character*5 afield
      character*6 NFIELD_FORM
      logical cond1, cond2
      integer i, j, iqn, iqn1, iqn2, i_n, i_m, ifield
      integer NM_LEN, MAX_NM_LEN, NFIELD_WIDTH
      integer n0(MAX_N_LEN*MAX_M_LEN), n0_(MAX_N_LEN*MAX_M_LEN)
      integer m0(MAX_N_LEN*MAX_M_LEN), m0_(MAX_N_LEN*MAX_M_LEN)
      
      parameter(NFIELD_WIDTH=4, NFIELD_FORM='(i3)')
      
      logical wfat_check_invec
      external wfat_check_invec

      
      MAX_NM_LEN = MAX_N_LEN * MAX_M_LEN
      iqn = 1
      i_n = 0
      
c     ==== Loop for reading n quantum numbers ====
      do i = 1, MAX_N_LEN
         
         if (.not. inp_a(afield)) then
            if (i_n .ne. 0) then
               exit
            else
               call errquit('WFAT: No n-parabolic quantum '//
     &              'number has been specified. Specify at '//
     &              'least one n-quantum number using the '//
     &              'format ''nx'' where ''x'' is zero or '//
     &              'positive integer.', 1, INPUT_ERR)
            endif
         endif

         cond1 = (afield(1:1) == 'n')
         cond2 = (verify(trim(adjustl(afield(2:NFIELD_WIDTH))),
     &           '0123456789') == 0)
         if (cond1 .and. cond2) then
            i_n = i_n + 1
            read(afield(2:NFIELD_WIDTH), NFIELD_FORM) n0_(i_n)

            i_m = 0
c           ==== Loop for reading m quantum numbers ====
            do j = 1, MAX_M_LEN
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
            NM_LEN = iqn2

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
            if (iqn > MAX_NM_LEN) then
               write(nelm_str, '(i4)')  MAX_NM_LEN
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


      if (.not. rtdb_put(rtdb, 'wfat:inp:n_npara', mt_int, 1, NM_LEN))
     &     call errquit('WFAT: Failed writing inp:n_npara to rtdb',
     &     32, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:npara', mt_int, NM_LEN, n0))
     &     call errquit('WFAT: Failed writing inp:npara to rtdb',
     &     33, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:mpara', mt_int, NM_LEN, m0))
     &     call errquit('WFAT: Failed writing inp:mpara to rtdb',
     &     34, INPUT_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_ffc(rtdb)   !(ftp, ergf)
c     Input syntax:
c       field_factor <float energy>
      implicit none
#include "wfat_const.fh"
#include "inp.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
      integer rtdb

      character(len=CHID_CHARLEN) ffc_s(MAX_N_LEN*MAX_M_LEN*MAX_N_ORBID)
      logical l1, l2
      integer nfields, n, i
      double precision tmp, ffc(MAX_N_LEN*MAX_M_LEN*MAX_N_ORBID)

      
      nfields = inp_n_field()
      if (nfields == 2) then
c        myffactor  <float value>
         l1 = inp_f(tmp)
         n = MAX_N_LEN * MAX_M_LEN * MAX_N_ORBID
         ffc_s(1:n) = 'all'
         ffc(1:n) = tmp
      else if (nfields > 2) then
c        myffactor \
c          "sid1 oid1 n1 m1"  <float value> \
c          "sid2 oid2 n2 m2"  <float value> \
c          ...
         if (mod(nfields,2) == 0) call errquit('WFAT: User-set field '//
     &        'factors are specified through pairs of channel IDs and'//
     &        ' the custom field factor corresponding to the channel.',
     &        5, INPUT_ERR)
         n = (nfields-1) / 2
         ffc_s = ''
         ffc = DEF_MFFC
         do i = 1, n
            l1 = inp_a(ffc_s(i))
            l2 = inp_f(ffc(i))
         enddo
      else
         call errquit('WFAT: Input directive myffactor has been given'//
     &        ' without any value.', 61, INPUT_ERR)
      endif
      
      if (.not. rtdb_put(rtdb, 'wfat:inp:n_myffc', mt_int, 1, n)) 
     &     call errquit('WFAT: Failed writing inp:n_myffc to RTDB',
     &     11, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'wfat:inp:myffc_s', n, ffc_s))
     &     call errquit('WFAT: Failed writing inp:myffc_s to RTDB',
     &     12, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:myffc', mt_dbl, n, ffc)) 
     &     call errquit('WFAT: Failed writing inp:myffc to RTDB',
     &     13, RTDB_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_w1x(rtdb)   !(nc, maxn, w1x, n)
c     Input syntax when using oe-wfat:
c       wfat1_exclude \
c          "sid1 oid1 n1 m1" \
c          "sid2 oid2 n2 m2" \
c          ...
c     where sid and oid are the spin and orbital IDs (integers) specified
c     through the movec_id directive, and n and m are parabolic quantum
c     numbers.
c     Input syntax when using me-wfat:
c       wfat1_exclude \
c          "n1 m1" \
c          "n2 m2" \
c          ...
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "wfat_const.fh"
      integer rtdb
      
      integer max_w1x, n, i
      character(len=CHID_CHARLEN) w1x(MAX_N_LEN*MAX_M_LEN*MAX_N_ORBID)
      
      max_w1x = MAX_N_LEN * MAX_M_LEN * MAX_N_ORBID
      w1x = ''
      n = 0
      do i = 1, max_w1x
         if (inp_a(w1x(i))) n = n + 1
      enddo

      if (n == 0) call errquit('WFAT: Input directive '//
     &     '''wfat1_exclude'' is given without entries.', 8, INPUT_ERR)

      if (.not. rtdb_put(rtdb, 'wfat:inp:n_wfat1_exclude', mt_int,
     &     1, n)) call errquit('WFAT: Failed writing '//
     &     'inp:n_wfat1_exclude to RTDB', 11, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'wfat:inp:wfat1_exclude', n, w1x))
     &     call errquit('WFAT: Failed writing inp:wfat1_exclude to '//
     &     'RTDB', 12, RTDB_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_mydip(rtdb)   !(nc, maxn, dip_s, dip, n)
c     DESCRIPTION:
c     Parse the user-specified unshifted moments of the ionized channels.
      implicit none
#include "inp.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "wfat_const.fh"
      integer rtdb
      
      character(len=ORBID_CHARLEN) dip_s(MAX_N_ORBID)
      integer n, i
      double precision tmp(3), dip(3,MAX_N_ORBID)

      
      if (inp_f(tmp(1))) then
c        mydipole  dipx dipy dipz
         n = MAX_N_ORBID
         if (.not. inp_f(tmp(2))) call errquit('WFAT: Cannot read '//
     &        'the y component of manual dipole.', 8, INPUT_ERR)
         if (.not. inp_f(tmp(3))) call errquit('WFAT: Cannot read '//
     &        'the z component of manual dipole.', 8, INPUT_ERR)
         do i = 1, MAX_N_ORBID
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
         do i = 1, MAX_N_ORBID            
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

      if (.not. rtdb_put(rtdb, 'wfat:inp:n_mydip', mt_int, 1, n)) 
     &     call errquit('WFAT: Failed writing inp:n_mydip to RTDB',
     &     11, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'wfat:inp:mydip_s', n, dip_s))
     &     call errquit('WFAT: Failed writing inp:mydip_s to RTDB',
     &     12, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:mydip', mt_dbl, 3*n, dip)) 
     &     call errquit('WFAT: Failed writing inp:mydip to RTDB',
     &     13, RTDB_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_myerg(rtdb)   ! (nc, maxn, erg_s, erg, n)

      implicit none
#include "inp.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "wfat_const.fh"
      integer rtdb

      character(len=ORBID_CHARLEN) erg_s(MAX_N_ORBID)
      integer n, i
      double precision tmp, erg(MAX_N_ORBID)
      

      if (inp_f(tmp)) then
c        myenergy  erg
         n = MAX_N_ORBID
         do i = 1, MAX_N_ORBID
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
         do i = 1, MAX_N_ORBID            
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

      if (.not. rtdb_put(rtdb, 'wfat:inp:n_myerg', mt_int, 1, n)) 
     &     call errquit('WFAT: Failed writing inp:n_myerg to RTDB',
     &     11, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'wfat:inp:myerg_s', n, erg_s))
     &     call errquit('WFAT: Failed writing inp:myerg_s to RTDB',
     &     12, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:inp:myerg', mt_dbl, 3*n, erg)) 
     &     call errquit('WFAT: Failed writing inp:myerg to RTDB',
     &     13, RTDB_ERR)
      
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


c=======================================================
c=======================================================
      subroutine wfat_parse_angle(rtdb, c)
c     Type 1:
c       beta <bt_start> <nbt> <bt_last>
c     bt_start, nbt, and bt_last are supplied by the user.
c     Type 2:
c       beta <bt_start>
c     This type is equivalent to bt <bt_start> 1 <bt_start>.
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"      
      integer rtdb
      character*(*) c

      integer nag
      double precision ag0, agn


      if (c /= 'beta' .and. c /= 'gamma') 
     &     call errquit('WFAT: Argument c of wfat_parse_angle '//
     &     'must be either beta or gamma.', 15, INPUT_ERR)

      
      if (inp_f(ag0)) then
         if (inp_i(nag)) then      ! type 1
            if (.not. inp_f(agn))
     &           call errquit('WFAT: A valid last angle is '//
     &           'expected.', 49, INPUT_ERR)
         else     ! type 2
            nag = 1
            agn = ag0
         endif
         
         if (.not. rtdb_put(rtdb, 'wfat:inp:'//c//'_start', mt_dbl, 1, 
     &        ag0)) call errquit('WFAT: Cannot write inp:'//c//
     &        '_start to rtdb', 53, INPUT_ERR)
         
         if (.not. rtdb_put(rtdb, 'wfat:inp:'//c//'_last', mt_dbl, 1,
     &        agn)) call errquit('WFAT: Cannot write inp:'//c//
     &        '_last to rtdb', 53, INPUT_ERR)
         
         if (.not. rtdb_put(rtdb, 'wfat:inp:n'//c, mt_int, 1, nag)) 
     &        call errquit('WFAT: Cannot write inp:n'//c//' to rtdb',
     &        53, INPUT_ERR)
      else
         call errquit('WFAT: Valid inputs for '//c//' directive is '//
     &        'expected.', 24, INPUT_ERR)
      endif

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_parse_momentum(rtdb, c)
c     Type 1:
c       kx <kx_start> <nkx> <kx_last>
c     kx_start, nkx, and kx_last are supplied by the user.
c     Type 2:
c       kx <kx_start>
c     This type is equivalent to kx <kx_start> 1 <kx_start>.
cOLD     Type 3:
cOLD       kx <nkx>
cOLD     This type is equivalent to kx -<kmax> <nkx> +<kmax>, where kmax is computed
cOLD     by NWChem according to the span of the broadest channel function.
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"      
      integer rtdb
      character*1 c

      integer nk
      double precision k0, kn


      if (c /= 'x' .and. c /= 'y') 
     &     call errquit('WFAT: Argument c of wfat_parse_momentum '//
     &     'must be either x or y.', 15, INPUT_ERR)

      
      if (inp_f(k0)) then
         if (inp_i(nk)) then      ! type 1
            if (.not. inp_f(kn))
     &           call errquit('WFAT: A valid last momentum point is '//
     &           'expected.', 49, INPUT_ERR)
         else     ! type 2
            nk = 1
            kn = k0
         endif

         if (.not. rtdb_put(rtdb, 'wfat:inp:k'//c//'_start', mt_dbl, 1, 
     &        k0)) call errquit('WFAT: Cannot write inp:k'//c//
     &        '_start to rtdb', 53, INPUT_ERR)
         
         if (.not. rtdb_put(rtdb, 'wfat:inp:k'//c//'_last', mt_dbl, 1, 
     &        kn)) call errquit('WFAT: Cannot write inp:k'//c//
     &        '_last to rtdb', 53, INPUT_ERR)
         
         if (.not. rtdb_put(rtdb, 'wfat:inp:nk'//c, mt_int, 1, nk)) 
     &        call errquit('WFAT: Cannot write inp:nk'//c//' to rtdb',
     &        53, INPUT_ERR)
      else
         call errquit('WFAT: Valid inputs for k'//c//' directive is '//
     &        'expected.', 24, INPUT_ERR)
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_intgl(rtdb)
c     DESCRIPTION:
c     Loads WFAT integrals from a file. Along with the restart directive, it can be
c     used to skip the calculation of WFAT integrals.
c     Syntax:
c       intgl load <file name to load> save <file name to save>
c     The "load  <file name to load>" and "save <file name to save>" are optional,
c     but one of them must exist when the "intgl" directive is given. When the
c     load file is not given, it will be ignored and the integrals will be calculated
c     by the program. If the save file is not given, an integral file with a default
c     name will be created in the same directory.
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer rtdb

      character*20 ss
      character*200 mode, fname
      logical lfound, sfound, l
      integer i

      lfound = .false.
      sfound = .false.
      
      do i = 1, 2
         l = inp_a(mode)
         l = inp_a(fname)
         if (inp_compare(.false., 'load', mode)) then
            if (lfound) then
               call errquit('WFAT: the load keywoard of the intgl '//
     &              'directive can only appear once.', 66, INPUT_ERR)
            else
               ss = 'wfat:inp:intgl_lname'
               lfound= .true.
               if (.not. rtdb_cput(rtdb, ss, 1, trim(adjustl(fname))))
     &              call errquit('WFAT: Failed writing '//ss//
     &              'to RTDB.', 1, INPUT_ERR)
            endif
         else if (inp_compare(.false., 'save', mode)) then
            if (sfound) then
               call errquit('WFAT: the save keywoard of the intgl '//
     &              'directive can only appear once.', 77, INPUT_ERR)
            else
               ss = 'wfat:inp:intgl_sname'
               sfound = .true.
               if (.not. rtdb_cput(rtdb, ss, 1, trim(adjustl(fname))))
     &              call errquit('WFAT: Failed writing '//ss//
     &              'to RTDB.', 1, INPUT_ERR)
            endif
         endif
      enddo

      call inp_prev_field()
      call inp_prev_field()

      if (.not. lfound .and. .not. sfound)
     &     call errquit('WFAT: When the intgl keyword is given, at '//
     &     'least one of the save or load files must be specified.',
     &     45, INPUT_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_inpa(rtdb, key)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"      
#include "mafdecls.fh"
      integer rtdb
      character*(*) key

      character*1000 val

      if (inp_a(val)) then
         if (.not. rtdb_cput(rtdb, 'wfat:inp:'//key, 1,
     &        trim(adjustl(val))) ) call errquit('WFAT: Cannot write '//
     &        'wfat:inp:'//key//' to RTDB.', 21, RTDB_ERR)
      else
         call errquit ("WFAT: "//key//"takes a string.", 20, INPUT_ERR)
      endif
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_parse_inpi(rtdb, key)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"      
#include "mafdecls.fh"
      integer rtdb
      character*(*) key

      integer val

      if (inp_i(val)) then
         if (.not. rtdb_put(rtdb, 'wfat:inp:'//key, mt_int, 1, val))
     &        call errquit('WFAT: Cannot write wfat:inp:'//key//
     &        ' to RTDB.', 21, RTDB_ERR)
      else
         call errquit ("WFAT: "//key//"takes an integer.",
     &        20, INPUT_ERR)
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_parse_inpf(rtdb, key)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"      
#include "mafdecls.fh"
      integer rtdb
      character*(*) key

      double precision val

      if (inp_f(val)) then
         if (.not.rtdb_put(rtdb, 'wfat:inp:'//key, mt_dbl, 1, val))
     &        call errquit('WFAT: Cannot write wfat:inp:'//key//
     &        ' to RTDB.', 21, RTDB_ERR)
      else
         call errquit ("WFAT: "//key//"takes a decimal value.", 
     &        20, INPUT_ERR)
      endif

      end
c=======================================================
      
      
      
      
     
