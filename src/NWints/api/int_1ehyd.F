c     =============== IMAM ===================
      
      subroutine int_1ehyd(i_basis,ish,j_basis,jsh,rcent,crg,lscr,scr,
cc AJL/Begin
     &    lvhd,vhd)
cc AJL/End
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geom.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
#include "util.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_hnd_1e
      external cando_nw_1e
      external cando_nw
c::passed
c:tex-\begin{verbatim}
      integer i_basis              !< [Input] basis set handle for ish
      integer ish                  !< [Input] i shell/contraction
      integer j_basis              !< [Input] basis set handle for jsh
      integer jsh                  !< [Input] j shell/contraction
      double precision rcent(1:3)  !< [Input] coordinate of potential center
      double precision crg         !< [Input] charge of potential center
      integer lscr                 !< [Input] length of scratch array
      integer lvhd                 !< [Input] length of integral buffer
      double precision scr(lscr)   !< [Scratch] scratch array
      double precision vhd(lvhd)   !< [Output] potential energy integrals
c:tex-\end{verbatim}
c::local

      logical doV
      logical ohnd_ok, onw_ok, onw1e_ok
      logical shells_ok, orel, oirel, ojrel, oNR, canAB
      integer i_geom, j_geom, ibas, jbas, ucont, uconts
      integer lbas, sbas, abas, isbas, jsbas
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_icfpS
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_icfpS
      integer ncents(1)
      double precision invexp(1), crg_(1)
c
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
c
      integer WarnP
      save WarnP
      data WarnP /0/
cc AJL/Begin
cc Soft-coded the value of msg
      character*20 msg
cc AJL/End
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('int_1eov'))
     &       call errquit('int_1eov: int_init was not called' ,0,
     &       INT_ERR)
c

      doV = .true.
      crg_(1) = crg
      ncents(1) = 1
      invexp(1) = 0.0d0

      
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1eov: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      oNR = .true.
      oirel = .false.
      ojrel = .false.
      orel = .false.
      canAB = .false.
c
      if (dyall_mod_dir) then
c
c     get basis set handles; relativistic integral option only valid
c     if both ibas and jbas are the ao basis.
c
        lbas = lc_bsh + BASIS_HANDLE_OFFSET
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        orel = ibas .eq. abas .and. jbas .eq. abas
      end if
c
c   i shell
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,ibas) .ne. 0)) then
        oirel = .true.
        isbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1estv: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          ibas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        isbas = ibas
      end if
c
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
      i_icfpS = infbs_cont(CONT_ICFP ,uconts,isbas)
c
c   j shell
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,jbas) .ne. 0)) then
        ojrel = .true.
        jsbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1estv: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          jbas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        jsbas = jbas
      end if
c
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
      j_icfpS = infbs_cont(CONT_ICFP ,uconts,jsbas)
c
      oNR = .not.(oirel.and.ojrel)
      orel = oirel.or.ojrel
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_1eov: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif

      ohnd_ok = cando_hnd_1e(i_basis,ish,0)
     &    .and. cando_hnd_1e(j_basis,jsh,0)
     &    .and. (.not.geom_any_finuc (i_geom))
     &    .and. (.not.geom_any_finuc (j_geom))
      onw_ok = cando_nw(i_basis,ish,0) .and. cando_nw(j_basis,jsh,0)
      onw1e_ok = cando_nw_1e(i_basis,ish,0)
     &     .and. cando_nw_1e(j_basis,jsh,0)

      
      if (orel) then
        call rel_onel (
     &      coords(1,i_cent,i_geom),
     &      dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfpS,isbas)),i_prim,i_gen,Li,
     &      coords(1,j_cent,j_geom),
     &      dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfpS,jsbas)),j_prim,j_gen,Lj,
     &      rcent(1),crg_(1),
     &      invexp(1),ncents(1),
     &      scr,scr,vhd,lvhd,.false.,.false.,doV,canAB,
     &      onw_ok,ohnd_ok,oNR,.false.,
     &      scr,lscr,rel_dbg,rel_typ)
      else if (ohnd_ok) then
        call hnd_stvint(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      j_prim, j_gen, Lj,
     &      rcent(1),crg_(1),ncents(1),
     &      scr,scr,Vhd,lvhd,.false.,.false.,doV,scr,lscr)
c
      elseif (onw1e_ok) then
         msg = 'random'
        call int_hf1sp(
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li, i_cent,
     &        coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfp,jbas)),
     &        j_prim, j_gen, Lj, j_cent,
     &        rcent(1),crg_(1),
     &        invexp(1),ncents(1),
     &        scr,scr,vhd,lvhd,.false.,.false.,doV,canAB,.false.,
cc AJL/Begin
c     &        scr,lscr,'int_1eov')
     &        scr,lscr,msg)
cc AJL/End
      elseif (onw_ok) then
        call hf1(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)), i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)), j_prim, j_gen, Lj,
     &      rcent(1),crg_(1),
     &      invexp(1),ncents(1),
     &      scr,scr,vhd,lvhd,.false.,.false.,doV,canAB,.false.,
     &      scr,lscr)
      else
        call errquit('int_1eov: could not do hnd, sp or nw integrals',
     &                0, INT_ERR)
      endif
c
*     We now have the cartesian integral block(s)  (jlo:jhi,ilo:ihi)
*      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4. Since int_nbf_* arrays are set to the appropriate size.
c
      if (li.eq.-1) i_gen = 1
      if (lj.eq.-1) j_gen = 1
c
      if (bas_spherical(ibas).and.bas_spherical(jbas)) then
*... transform both i and j integrals
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
c
        if (doV) call spcart_tran1e(vhd,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(ibas)) then
*.. transform on i component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = j_nbf_x
        if (doV) call spcart_tran1e(vhd,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(jbas)) then
*.. transform on j component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = i_nbf_x
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
        if (doV) call spcart_tran1e(vhd,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
      else
        call errquit(
     &      'int_1ehyd: should never reach transform blocked else',911,
     &       INT_ERR)
      endif
      return
c
      end
c     =============== IMAM ===================


c=======================================================
c=======================================================
      subroutine int_1ehyd_all(ibas, nbas, rcent, crg, vhyd)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"      
      integer ibas                 ! [input]
      integer nbas                 ! [input]
      double precision rcent(1:3)  ! [Input] coordinate of potential center
      double precision crg         ! [Input] charge of potential center
      double precision vhyd(nbas, nbas)   ! [output]

      integer l_s, l_hb
      integer k_s, k_hb
      integer nbas0, max1e, mem1, ncont, nib, njb, nblock
      integer i, j, ib1, ib2, jb1, jb2



c     ==== Initial checks ====
      if (.not. bas_numbf(ibas, nbas0)) then
         call errquit('int_1ehyd_all: Getting nbas0 failed',
     &        53, BASIS_ERR)
      else
         if (nbas0 /= nbas) then
            call errquit('int_1ehyd_all: Argument nbas is '//
     &           'inconsistent with handle ibas', 53,
     &           BASIS_ERR)
         endif
      endif
      
c     ==== Some allocations ====
      call int_mem_1e(max1e, mem1)
      if (.not. MA_alloc_get(MT_DBL, mem1, 'scr00', l_s, k_s))
     &  call errquit('int_1ehyd_all: 1e scratch allocation failed',
     &     31, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, max1e, 'vhyd_bl', l_hb, k_hb))
     &  call errquit('int_1ehyd_all: vhyd_bl allocation failed',
     &     32, MA_ERR)

c     ==== Begin filling vhyd ====
      if (.not. bas_numcont(ibas, ncont))
     &     call errquit('int_1ehyd_all: Getting ncont failed',
     &          54, BASIS_ERR)

      vhyd = 0.0d0
      do j = 1, ncont
         if (.not. bas_cn2bfr(ibas, j, jb1, jb2))
     &        call errquit('int_1ehyd_all: a cn2bfr call failed',
     &        57, BASIS_ERR)
         njb = jb2 - jb1 + 1
         
         do i = 1, j
            if (.not. bas_cn2bfr(ibas, i, ib1, ib2))
     &           call errquit('int_1ehyd_all: a cn2bfr call failed',
     &           57, BASIS_ERR)
            nib = ib2 - ib1 + 1

c           ==== Calculate a block of vhyd ====
            dbl_mb(k_s : k_s+mem1-1) = 0.0d0            
cOLD            call int_1ehyd(ibas, i, ibas, j, rcent, crg,
cOLD     &           mem1, dbl_mb(k_s), max1e, dbl_mb(k_hb))
            call int_1ehyd(ibas, j, ibas, i, rcent, crg,
     &           mem1, dbl_mb(k_s), max1e, dbl_mb(k_hb))

c           ==== Fill the corresponding block of vhyd ====
            nblock = nib * njb
            call int_1ehyd_fill(nib, njb, dbl_mb(k_hb:k_hb+nblock-1),
     &           vhyd(ib1:ib2, jb1:jb2))
         enddo
      enddo

c     ==== Apply the Hermitian property ====
      do j = 1, nbas
         do i = 1, j
            vhyd(j,i) = vhyd(i,j)
         enddo
      enddo
      vhyd = -vhyd     ! 1)
                       
                       
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_s))
     &     call errquit('int_1ehyd_all: scratch deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_hb))
     &     call errquit('int_1ehyd_all: vhyd_bl deallocation failed',
     &     4, MA_ERR)

      end

c     NOTES :
c       1) This sign shift is needed because the routine int_1ehyd is basically  
c          scrapped from the routine that computes the matrix elements of nuclear
c          attraction operator which bears a minus sign. And the minus sign there
c          is apparently absorbed into the matrix elements. So here, the sign    
c          shift negates that absorbed minus sign.                               
c=======================================================
      

c=======================================================
c=======================================================
      subroutine int_1ehyd_fill(m, n, vhbl, vh)

      implicit none
      integer m                   ! [input]
      integer n                   ! [input]
      double precision vhbl(m,n)  ! [input] 
      double precision vh(m,n)    ! [output]

      
      vh = vhbl
      
      end
c=======================================================      
