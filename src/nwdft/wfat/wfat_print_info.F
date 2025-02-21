c=======================================================
c=======================================================
      subroutine wfat_module_title(un)

      implicit none
      integer un

      write(un, '(4/)')
      call util_print_centered(
     &     un, '-------------------------------------------------',
     &     40, .false.)
      call util_print_centered(
     &     un, '<<<<<< Weak-Field Asymptotic Theory Module >>>>>>',
     &     40, .false.)
      call util_print_centered(
     &     un, '-------------------------------------------------',
     &     40, .false.)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_print_info(rtdb)

      implicit none
#include "global.fh"      
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "wfat_qm_param.fh"
#include "geom.fh"
#include "bas.fh"      
      integer rtdb

      logical rtdb_stat
      integer ndid

      
      ndid = ga_nodeid()
      if (ndid /= 0) return

      rtdb_stat = rtdb_parallel(.false.)
      
      
c     ==== Channel info ====
cOLD      call wfat_print_channel(rtdb, nel)
      

c     ==== Printing info ====
      call wfat_print_printsw(rtdb)      
      

c     ==== Orientation angles ====
      call wfat_print_angles(rtdb)
      
            
c     ==== Remaining WFAT info ====
      call wfat_print_others(rtdb)

      
      rtdb_stat = rtdb_parallel(.true.)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_level(rtdb)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb
      
      character*2 wlevel

      if (.not. rtdb_cget(rtdb,'wfat:level', 1, wlevel)) 
     &     call errquit('WFAT: Failed reading wfat level',0, RTDB_ERR)
      write(luout, '(2/)')
      call util_print_centered(luout, 'Level of calculation',
     &     40, .true.)
      if (wlevel == 'oe') then
         call util_print_centered(luout, 'One-electron WFAT',
     &        40, .false.)
      else if (wlevel == 'me') then
         call util_print_centered(luout, 'Many-electron WFAT',
     &        40, .false.)
      else
         call errquit('WFAT: WFAT level unknown.', 1, UERR)
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_qmsys(rtdb)

      implicit none
#include "wfat_qm_param.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "global.fh"
      integer rtdb

      logical status
      character*255 movecs, scftype, basname, title
      integer nmo(2), nb, ns

      logical movecs_read_header
      external movecs_read_header
      
      
      if (.not. rtdb_cget(rtdb, 'wfat:movec_name', 1, movecs))
     &     call errquit('WFAT: Cannot read wfat:movec_name from RTDB.',
     &     1, RTDB_ERR)

      status = movecs_read_header(movecs, title, basname, scftype,
     &         nb, ns, nmo, 2)

      
      if (ga_nodeid() == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Previous SCF info', 
     &        20, .true.)
         write(luout, '("")')
      endif

      if (ga_nodeid() == 0)
     &     write(luout, 900) 'initial species (neutral)',
     &     trim(movecs), trim(title), trim(scftype), trim(basname),
     &     nbas0, nspin0, nmo, nel, multip, charge

 900  format(
     &     3x, 'Information about the ', a, ' :', /,
     &     5x, '* MO file = ', a, /,
     &     5x, '* Job title = ', a, /,
     &     5x, '* SCF type = ', a, /,
     &     5x, '* Basis set name = ', a, /,
     &     5x, '* No. of bases = ', i5, /,
     &     5x, '* No. of spin channels = ', i3, /,
     &     5x, '* No. of MO''s in each channel = ', 2i5, /,
     &     5x, '* No. electrons in each channel = ', 2i4, /,
     &     5x, '* Spin multiplicity = ', i3, /,
     &     5x, '* Charge state = ', f7.3)
            
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_channels(nch, spinid, orbid)

      implicit none
#include "stdio.fh"
      integer nch, spinid(nch), orbid(nch)
      
      character*5 spinnm(2)
      character*7 io_str
      character*12 hl
      integer i, is, io, max_in_row

      parameter(max_in_row=4)

      
      spinnm = (/ 'alpha', ' beta' /)      
      write(luout, '(5x, "No. of channels = ", i6)') nch
      write(luout, '(5x, "Ionized channels = ")', advance='no')

      do i = 1, nch
         is = spinid(i)
         io = orbid(i)
         if (io == 0) then
            hl = 'HOMO'
         else if (io < 0) then
            write(io_str, '(i7)') abs(io)
            hl = "HOMO-"//adjustl(io_str)
         else if (io == 1) then
            hl = 'LUMO'
         else if (io > 1) then
            write(io_str, '(i7)') abs(io-1)
            hl = "LUMO+"//adjustl(io_str)
         endif

         write(luout, '(a, "/", a)', advance='no')
     &        trim(adjustl(spinnm(is))), trim(adjustl(hl))
         if (i < nch) then
            write(luout, '(", ")', advance='no')
         else if (i == nch) then
            write(luout, '("")')
         endif

         if (mod(i,max_in_row) == 0) then
            write(luout, '("")')
            write(luout, '(5x, 19x, "")', advance='no')
         endif
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_paraqn(nqn, nid, mid)

      implicit none
#include "stdio.fh"
      integer nqn, nid(nqn), mid(nqn)

      character*6 n_str, m_str
      integer i, max_in_row

      parameter(max_in_row=6)


      write(luout, '(5x, "No. of parabolic quantum numbers = ", i6)')
     &     nqn
      write(luout, '(5x, "Parabolic quantum numbers = ")',
     &     advance='no')
      do i = 1, nqn
         write(n_str, '(i5)') nid(i)
         write(m_str, '(sp, i5)') mid(i)
         
         write(luout, '(a)', advance='no')
     &        '('//trim(adjustl(n_str))//','//trim(adjustl(m_str))//
     &        ')'
         write(luout, '("; ")', advance='no')

         if (mod(i,max_in_row) == 0) then
            write(luout, '("")')
            write(luout, '(5x, 28x, "")', advance='no')
         endif
      enddo
      write(luout, '("")')

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_printsw(rtdb)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      
      character*3 ispr      
      logical pr_srfc, pr_rate, pr_yild, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom
      

      write(luout, '(2/)')
      call util_print_centered(luout, 'Printing Information', 20,
     &     .true.)
      write(luout, '(1x)')

      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)

      call wfat_tf_to_yn(pr_srfc, .true., ispr)
      write(luout, '(5x, "Print structure factor : ", a3)') ispr
      call wfat_tf_to_yn(pr_rate, .true., ispr)
      write(luout, '(5x, "Print ion rate : ", a3)') ispr
      call wfat_tf_to_yn(pr_yild, .true., ispr)
      write(luout, '(5x, "Print ion yield : ", a3)') ispr
      call wfat_tf_to_yn(pr_avgamma, .true., ispr)
      write(luout, '(5x, "Print gamma-averaged ion yield : ", a3)') ispr
      call wfat_tf_to_yn(pr_avbeta, .true., ispr)
      write(luout, '(5x, "Print beta-averaged ion yield : ", a3)') ispr
      call wfat_tf_to_yn(pr_tmd, .true., ispr)
      write(luout, '(5x, "Print transverse momentum distribution '//
     &     '(TMD) : ", a3)') ispr
      call wfat_tf_to_yn(pr_tmd_atom, .true., ispr)
      write(luout, '(5x, "Print atom-resolved TMD : ", a3)') ispr

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_angles(rtdb)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      
      logical pr_srfc, pr_rate, pr_yild, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom
      integer nbt, ngm
      double precision bt_start, bt_last, gm_start, gm_last
      double precision bt_tmd, gm_tmd
      

      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)


      write(luout, '(2/)')
      call util_print_centered(luout, 'Orientation angle '//
     &     'information', 20, .true.)
      write(luout, '(1x)')
      if (pr_srfc .or. pr_rate) then
c        ==== Beta ====
         if (.not. rtdb_get(rtdb, 'wfat:beta_start', mt_dbl, 1,
     &        bt_start))
     &        call errquit('WFAT: Failed reading beta_start from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "First beta (deg) ", bt_start
         if (.not. rtdb_get(rtdb, 'wfat:beta_last', mt_dbl, 1,
     &        bt_last))
     &        call errquit('WFAT: Failed reading beta_last from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "Last beta (deg) ", bt_last
         if (.not. rtdb_get(rtdb, 'wfat:nbeta', mt_int, 1, nbt))
     &        call errquit('WFAT: Failed reading nbeta from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", i7)')
     &        "Number of beta points  ", nbt

c        ==== Gamma ====
         if (.not. rtdb_get(rtdb, 'wfat:gamma_start', mt_dbl, 1,
     &        gm_start))
     &        call errquit('WFAT: Failed reading gamma_start from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "First gamma (deg) ", gm_start
         if (.not. rtdb_get(rtdb, 'wfat:gamma_last', mt_dbl, 1,
     &        gm_last))
     &        call errquit('WFAT: Failed reading gamma_last from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "Last gamma (deg) ", gm_last
         if (.not. rtdb_get(rtdb, 'wfat:ngamma', mt_int, 1, ngm))
     &        call errquit('WFAT: Failed reading ngamma from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", i7)')
     &        "Number of gamma points ", ngm
      endif
      
      if (pr_tmd .or. pr_tmd_atom) then         
         if (.not. rtdb_get(rtdb, 'wfat:beta_tmd', mt_dbl, 1, bt_tmd))
     &        call errquit('WFAT: Failed reading bt_tmd from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "Beta (deg) for TMD ", bt_tmd
         if (.not. rtdb_get(rtdb, 'wfat:gamma_tmd', mt_dbl, 1, gm_tmd))
     &        call errquit('WFAT: Failed reading gm_tmd from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "gamma (deg) for TMD ", gm_tmd
      endif

      write(luout, '(/, 5x, "Note :")')
      write(luout, 601)
 601  format(7x,
     &     "To get an idea of how the molecule is oriented with "     ,
     &     "the above values of orientation ",   /,
     &     7x, 
     &     "angles, imagine the following rotations. First imagine "  ,
     &     "a pair of coinciding axes, the ", /,
     &     7x, 
     &     "(space-fixed) lab axes and the (body-fixed) internal "    ,
     &     "axes. Then rotate the molecule around ", /,
     &     7x, 
     &     "the lab z axis (or internal z, whichever is fine at "     ,
     &     "this point) counterclockwise as seen ",        /,
     &     7x, 
     &     "from lab +z axis by an ammount gamma. Thus, at this "  ,
     &     "point, gamma is the angle subtended ",  /,
     &     7x, 
     &     "by the lab and internal y axes. Next, rotate the ",
     &     "molecule around the lab y axis counter ",/,
     &     7x, 
     &     "clockwise as seen from lab +y axis by an ammount beta ",
     &     "so that beta is the angle between ", /,
     &     7x,
     &     "the lab and internal z axes. These rotations yield the ",
     &     "molecule's final orientation.")

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_others(rtdb)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      
      character*255 fname
      logical pr_srfc, pr_rate, pr_yild, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom, autopos
      integer nktmd, lmax, gqual
      double precision fieldmag, omgpos(3)
      

      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)


      write(luout, '(2/)')
      call util_print_centered(luout, 'Other WFAT inputs', 30, .true.)
      write(luout, '(1x)')

cOLDC     ==== Field magnitude ====
cOLD      if (pr_rate .or. pr_tmd .or. pr_tmd_atom) then
cOLD         if (.not. rtdb_get(rtdb, 'wfat:fieldmag', mt_dbl, 1, fieldmag))
cOLD     &        call errquit('WFAT: Failed reading fieldmag from rtdb',
cOLD     &        1, INPUT_ERR)
cOLD         write(luout, '(5x, a28, " = ", es10.3, " a.u.")')
cOLD     &        "Field strength ", fieldmag
cOLD      endif


c     ==== No. of points for TMD ====
      if (pr_tmd .or. pr_tmd_atom) then
         if (.not. rtdb_get(rtdb, 'wfat:nktmd', mt_int, 1, nktmd))
     &        call errquit('WFAT: Failed reading nktmd from rtdb',
     &        1, RTDB_ERR)
         write(luout, '(5x, a28, " = ", i7)')
     &        "No. of x momentum points ", nktmd
         write(luout, '(5x, a28, " = ", i7)')
     &        "No. of y momentum points ", nktmd
      endif

c     ==== Max. angular momentum ====
      if (.not. rtdb_get(rtdb, 'wfat:lmax', mt_int, 1, lmax))
     &     call errquit('WFAT: Failed reading lmax from rtdb',
     &     1, INPUT_ERR)
      write(luout, '(5x, a28, " = ", i3)')
     &     "Max. angular momentum ", lmax

c     ==== Grid quality ====
      if (.not. rtdb_get(rtdb, 'wfat:gqual', mt_int, 1, gqual))
     &     call errquit('WFAT: Failed reading gqual from rtdb',
     &     1, INPUT_ERR)
      write(luout, '(5x, a28, " = ", i3)') "Grid set ", gqual

c     ==== Integral file ====
      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('WFAT: Cannot read wfat:intgl_name from RTDB',
     &     1, RTDB_ERR)
      write(luout, '(5x, a28, " = ", a)')
     &     "Integral file ", trim(adjustl(fname))

c     ==== Auto/manual positioning of origin ====
      if (.not. rtdb_get(rtdb, 'wfat:autopos', mt_log, 1, autopos))
     &     call errquit('WFAT: Cannot read wfat:autopos from RTDB.',
     &     1, RTDB_ERR)
      if (autopos) then
         write(luout, '(5x, a28, " = ", a)') "Channel pos. ", 'auto'
      else
         if (.not. rtdb_get(rtdb, 'wfat:omgpos', mt_dbl, 3, omgpos))
     &     call errquit('WFAT: Cannot read wfat:omgpos from RTDB.',
     &     1, RTDB_ERR)
         write(luout, '(5x, a28, " = ", 3es14.5)')
     &        "Channel pos. ", omgpos(1:3)

         write(luout, '(/, 7x, "*** ATTENTION ***")')
         write(luout, 602)
      endif
      
 602  format(9x,
     &     "Manual specification for the input 'omgpos' such "      ,
     &     "as in the present run should be used with caution. ",    /,
     &     9x, 
     &     "This feature is only intended for research purposes. "    ,
     &     "Improper choice of the values for 'omgpos' can lead ", /,
     &     9x, 
     &     "to a behavior of the results that is not in agreement, "  ,
     &     "even qualitatively, with results from other methods. ",  /,
     &     9x, 
     &     "To let NWChem pick the right values for the user, "       ,
     &     "use 'auto' as the value for 'omgpos'. ")

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_print_mask(l, r)

      implicit none
#include "stdio.fh"
      double precision l    ! Smoothing length in bohr
      double precision r    ! Radius in bohr

      double precision ang2bohr
      parameter (ang2bohr=1.8897261255d0)
      
      write(luout, '(2/)')
      call util_print_centered(luout, 'Mask Information', 38, .true.)
      write(luout, '(1x)')
      write(luout, '(5x, "Masking is active. The integrand of the '//
     &     'WFAT integrals will be masked with ", /, '//
     &     '5x, "the following masking parameters:")')
      write(luout, '(7x, "Smoothing length = ", f8.4, " bohr", '//
     &   '" = ", f8.4, " angstrom")') l, l/ang2bohr
      write(luout, '(7x, "Radius = ", f8.4, " bohr", " = ", '//
     &   'f8.4, " angstrom")') r, r/ang2bohr
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_polb_tensors(nch, pb)

      implicit none
#include "stdio.fh"
      integer nch
      double precision pb(3,3,nch)

      character*1 c(3)
      integer i, j, k
      
      c = (/ 'x', 'y', 'z'/)
      
      write(luout, '(2/)')
      call util_print_centered(luout, 'Polarizability tensors', 27,
     &     .true.)
      write(luout, '("")')
      do i = 1, nch
         write(luout, '(5x, "Channel : ", i4)') i

         write(luout, '(7x, 1x)', advance='no')
         do k = 1, 3
            write(luout, '(2x, a13)', advance='no') c(k)
         enddo
         write(luout, '("")')
         
         do j = 1, 3
            write(luout, '(7x, a)', advance='no') c(j)
            do k = 1, 3
               write(luout, '(2x, es13.6)', advance='no')
     &              pb(j,k,i)
            enddo
            write(luout, '("")')
         enddo

         write(luout, '("")')
      enddo

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_print_o(nch, nqn, nid, mid, o)

      implicit none
#include "stdio.fh"
      integer nch
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      double precision o(2,nqn,nch)

      character*6 n_str, m_str
      integer i, j

      write(luout, '(2/)')
      call util_print_centered(luout, 'O coefficients', 23,
     &     .true.)
      write(luout, '("")')

      do i = 1, nch
         write(luout, '(5x, "Channel : ", i4)') i
         write(luout, '(5x, 2a4, 2a19)')
     &        'n', 'm', 'Constant part', 'Orientation part'
         call wfat_print_repeat(luout, 46, 5, '-')
         write(luout, '("")')
         do j = 1, nqn
            write(luout, '(5x, 2i4, 2es19.8)')
     &           nid(j), mid(j), o(1,j,i), o(2,j,i)
         enddo
         write(luout, '("")')
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_w1x(nch, oid, sid, nqn, nid, mid, w1x)

      implicit none
#include "stdio.fh"
      integer nch
      integer oid(nch)
      integer sid(nch)
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      logical w1x(nqn,nch)

      character*12 i_s, n_s, m_s, dt_s
      integer i, j


      write(luout, '(2/)')
      call wfat_info_banner(luout, .true.)
      if (any(w1x))
     &     write(luout, '(5x, a)') 'Some channels and parabolic '//
     &     'quantum numbers are treated with WFAT(0), they are:'

      do i = 1, nch
         do j = 1, nqn
            if (w1x(j,i)) then
               write(i_s, '(i12)') i
               write(n_s, '(i12)') nid(j)
               write(m_s, '(i12)') mid(j)
               
               write(luout, '(7x, "Channel: ", a, ",  n = ", '//
     &              'a, ",  m = ", a)')
     &              trim(adjustl(i_s)), trim(adjustl(n_s)),
     &              trim(adjustl(m_s))
            endif
         enddo
      enddo

      end
c=======================================================
      
      
