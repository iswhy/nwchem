
      subroutine imam_test()

      implicit none

#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"

c      when using this, replace all instances of 'ao_bas_han' with 'basis'.
#include "cscf.fh"
      
c     errquit.fh contains GEOM_ERR, BASIS_ERR
#include "errquit.fh"
#include "geom.fh"

c     cdft.fh contains geom, AO_bas_han
c#include "cdft.fh"   when using this, replace all instances of 'basis' with 'ao_bas_han'.

      character*50 fr
      character*16 tag
      character*3 npt_str
      integer natom, nbas, ncont, nbl, ang, nxp, ncfs, casp
      integer mem1, max1e, l_pe, k_pe, l_s, k_s, l_hy, k_hy
      integer i, j, k, ico, icent, ico1, ico2, ibf1, ibf2, icol
      integer ico_a, ico_b, ibf_a1, ibf_a2, ibf_b1, ibf_b2
      integer npoint, lgr_xyz, igr_xyz, latomcd, iatomcd
      integer latomch, iatomch, ltags, itags, lrq, irq, ipt
      integer lniz, iniz, ldocset, idocset
      integer lbas_cent_info, ibas_cent_info
      integer lbas_cset_info, ibas_cset_info
      integer lchi_ao, ichi_ao, lmocf, kmocf, lmo_pt, imo_pt, ibf
      integer imo, i1, i2
      integer nbas_a, nbas_b
      integer ndid
      double precision znuc, crg
      double precision c(3), rcent(3)
      double precision dx, dy, dz, mox



      

      
c     ==== Get node ID ====
      ndid = ga_nodeid()
cforcheck      if (ndid == 0) then
cforcheck         write(LuOut, '("====== CHECKING NO. OF ELECTRONS: =====")') 
cforcheck         write(LuOut, '("noc(1) = ", i4)') noc(1)
cforcheck         write(LuOut, '("noc(2) = ", i4)') noc(2)
cforcheck         write(LuOut, '("ipol = ", i4)') ipol
cforcheck         write(LuOut, '(" ")')
cforcheck      endif
      
      if (ndid == 0) then
         write(LuOut, '("====== CHECKING NO. OF ELECTRONS: =====")') 
         write(LuOut, '("nalpha = ", i4)') nalpha
         write(LuOut, '("nbeta = ", i4)') nbeta
         write(LuOut, '(" ")')
      endif
      

      
c     ==== Get no. of atoms ====
      if (.not. geom_ncent(geom, natom))
     &     call errquit('imam: testing natom failed', 51, GEOM_ERR)

      
c     ==== Loop over atoms ====
      do icent = 1,natom
         if (.not. geom_cent_get(geom,icent,tag,c(1),znuc))
     &        call errquit('imam: testing znuc failed', 52, GEOM_ERR)
         if (.not. bas_ce2bfr(basis, icent, ibf1, ibf2))
     &        call errquit('imam: testing ce2bfr', 60, BASIS_ERR)
         if (.not. bas_ce2cnr(basis, icent, ico1, ico2))
     &        call errquit('imam: testing ce2cnr', 61, BASIS_ERR)
         
         if (ndid.eq.0) write(LuOut, 400)
     &        icent, tag, znuc, c(1), c(2), c(3), ibf1, ibf2, ico1, ico2
      enddo
 400  format(4x, 'Center ID =  ', i2/,
     &       6x, 'Element = ', a2/,
     &       6x, 'Nuc. charge = ', f5.2/,
     &       6x, 'Coordinate = ', 3f8.4/,
     &       6x, 'First basis = ', i3/,
     &       6x, 'Last basis = ', i3/,
     &       6x, 'First contraction = ', i3/,
     &       6x, 'Last contraction = ', i3)


c     ==== Get no. of basis functions and contractions ====
      if (.not. bas_numbf(basis, nbas))
     &     call errquit('imam: testing nbas failed', 53, BASIS_ERR)
      if (.not. bas_numcont(basis, ncont))
     &     call errquit('imam: testing ncont failed', 54, BASIS_ERR)

      if (ndid.eq.0) write(LuOut, 401) nbas, ncont
 401  format(/,
     &       6x, 'No. of basis functions = ', i3/,
     &       6x, 'No. of contractions = ', i4/)
      

c     ==== Loop over contractions ====
      do ico = 1, ncont
         if (.not. bas_continfo(basis, ico, ang, nxp, ncfs, casp))
     &        call errquit('imam: testing continfo failed',
     &        55, BASIS_ERR)
         if (.not. bas_cn2ce(basis, ico, icent))
     &        call errquit('imam: testing cn2ce failed',
     &        56, BASIS_ERR)
         if (.not. bas_cn2bfr(basis, ico, ibf1, ibf2))
     &        call errquit('imam: testing cn2bfr failed',
     &        57, BASIS_ERR)
         
         if (ndid.eq.0) write(LuOut, 402) ico, ang, nxp, ncfs, casp,
     &        icent, ibf1, ibf2
      enddo
 402  format(4x, 'Contraction ID = ', i2/,
     &       6x, 'Ang. momentum = ', i2/,
     &       6x, 'No. of exponents = ', i3/,
     &       6x, 'No. of radial functions = ', i3/,
     &       6x, 'Cartesian/spherical ? ', i2/,
     &       6x, 'Center ID = ', i3/,
     &       6x, 'First basis = ', i3/,
     &       6x, 'Last basis = ', i3)



c========== Block of shell doublet size ========      
      ico_a = 6
      ico_b = 6
      if (.not. bas_cn2bfr(basis, ico_a, ibf_a1, ibf_a2))
     &        call errquit('imam: 2nd testing cn2bfr failed',
     &     57, BASIS_ERR)
      if (.not. bas_cn2bfr(basis, ico_b, ibf_b1, ibf_b2))
     &        call errquit('imam: 2nd testing cn2bfr failed',
     &        57, BASIS_ERR)

      call int_mem_1e(max1e, mem1)
      nbas_a = ibf_a2 - ibf_a1 + 1
      nbas_b = ibf_b2 - ibf_b1 + 1
      nbl = nbas_a * nbas_b
      write(luout,'("IMAM: shell block size = ", i5)') nbl
      write(luout,'("IMAM: max buffer size = ", i5)') max1e
      write(luout,'("IMAM: 1e scratch size = ", i5)') mem1
      if (nbl > max1e) 
     &     call errquit('imam: pe buffer allocation failed',
     &     58, BASIS_ERR)
         

c============== Scratch array preparation ==============
      if (.not. MA_push_get(MT_DBL, mem1, 'scr', l_s, k_s))
     &  call errquit('imam: 1e scratch allocation failed',
     &     59, MA_ERR)

      
c================ Pe ===============      
      if (.not. MA_push_get(MT_DBL, max1e, 'pe', l_pe, k_pe))
     &  call errquit('imam: pe buffer allocation failed',
     &     58, MA_ERR)      
      call int_1epe(basis, ico_a, basis, ico_b,
     &     mem1, dbl_mb(k_s), max1e, dbl_mb(k_pe))

      write(luout,'("IMAM : Checking PE")')
      write(fr,'(a, i3, 15a)') '(', nbas_b, '(2x, e11.5))'
      k = k_pe
      do icol = 1, nbas_a
         if (ndid.eq.0) write(luout,FMT=fr) dbl_mb(k : k+nbas_b-1)
         k = k + nbas_b
      enddo

c================ Vhyd ===============
      rcent(1:3) = 0.0d0
      rcent(1) = 2.0d0
      rcent(2) = 1.0d0
      rcent(3) = -3.0d0
      crg = 1.0d0
      if (.not. MA_push_get(MT_DBL, max1e, 'vhyd', l_hy, k_hy))
     &  call errquit('imam: pe buffer allocation failed',
     &     58, MA_ERR)
      call int_1ehyd(basis, ico_a, basis, ico_b,
     &     rcent, crg, mem1, dbl_mb(k_s), max1e, dbl_mb(k_hy))

      write(luout,'("IMAM : Checking Vhyd")')
      k = k_hy
      do icol = 1, nbas_a
         if (ndid.eq.0) write(luout,FMT=fr) dbl_mb(k : k+nbas_b-1)
         k = k + nbas_b
      enddo      

c================MO at grids ===============
c     ==== Calculating gr_xyz ====
      npoint = 4
      if (.not. ma_push_get (mt_dbl, 3*npoint, "gr_xyz",
     &     lgr_xyz, igr_xyz))
     $     call errquit ("imam: cannot alloc grid", 0, MA_ERR)
c     == 1st point ==
      dbl_mb(igr_xyz+0) = 1.5d0
      dbl_mb(igr_xyz+1) = 3.6d0
      dbl_mb(igr_xyz+2) = 4.1d0
c     == 2nd point ==
      dbl_mb(igr_xyz+3) = -2.8d0
      dbl_mb(igr_xyz+4) = -1.3d0
      dbl_mb(igr_xyz+5) = -0.4d0
c     == 3rd point ==
      dbl_mb(igr_xyz+6) = -2.8d0
      dbl_mb(igr_xyz+7) =  1.3d0
      dbl_mb(igr_xyz+8) = -0.4d0
c     == 4th point ==
      dbl_mb(igr_xyz+9)  =  1.5d0
      dbl_mb(igr_xyz+10) = -3.6d0
      dbl_mb(igr_xyz+11) =  4.1d0

      
c     ==== Calculating atomcd ====
      if (.not.ma_push_get(mt_dbl,3*natom,'atomcd', latomcd, iatomcd))
     &     call errquit('imam: cannot allocate xyz_atom',0, MA_ERR)
      
      if (.not.ma_push_get(mt_dbl,natom,'atomch',latomch,iatomch))
     &     call errquit('imam: cannot allocate charge',0, MA_ERR)
      
      if (.not.ma_push_get(mt_Byte,natom*16,'tags',ltags,itags))
     &     call errquit('imam: cannot allocate tags',0, MA_ERR)
      
      if (.not. geom_cart_get(geom, natom, byte_mb(itags),
     &     dbl_mb(iatomcd), dbl_mb(iatomch)))
     &     call errquit('imam: geom_cart_get failed', 0, GEOM_ERR)

      
c     ==== Calculating rq ====
      if(.not.ma_push_get(mt_dbl, npoint*natom, 'rq',lrq,irq))
     &     call errquit("imam: rq",0, MA_ERR)

      do icent = 1, natom
         do ipt = 1, npoint
            dx = dbl_mb(iatomcd+3*(icent-1)+0) -
     &           dbl_mb(igr_xyz+3*(ipt-1)+0)
            dy = dbl_mb(iatomcd+3*(icent-1)+1) -
     &           dbl_mb(igr_xyz+3*(ipt-1)+1)
            dz = dbl_mb(iatomcd+3*(icent-1)+2) -
     &           dbl_mb(igr_xyz+3*(ipt-1)+2)
            dbl_mb(irq + npoint*(icent-1) + ipt-1) =
     &           dx*dx + dy*dy + dz*dz
         enddo
      enddo

      
c     ==== Calculating niz ====
      if(.not.ma_push_get(MT_int, natom, 'niz', lniz, iniz))
     &     call errquit("imam: niz",0, MA_ERR)
      do icent = 1, natom
         int_mb(iniz + icent - 1) = 1
      enddo

      
c     ==== Calculating docset ====
      if (.not.ma_push_get(mt_log, ncont, 'docset', ldocset, idocset))
     &     call errquit('imam: cannot allocate ccdocset', ncont, MA_ERR)
      do ico = 1, ncont
         log_mb(idocset + ico - 1) = .true.
      enddo

      
c     ==== Calculating bas_cent_info and bas_cset_info ====
      if (.not.ma_push_get(mt_int, 3*natom, "bas_cent_info",
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit("imam: cannot allocate bas_cent_info",0,
     &     MA_ERR)
      
      if (.not.ma_push_get(mt_int, 6*ncont, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit("imam: cannot allocate bas_cset_info",0,
     &     MA_ERR)
      
      call xc_make_basis_info(basis, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natom)


c     ==== Allocating chi_ao ====
      if (.not.ma_push_get(mt_dbl, npoint*nbas,
     &     'chi_ao', lchi_ao, ichi_ao))
     &     call errquit('imam: cannot allocate chi_ao',0, MA_ERR)

      
c     ==== Calculating chi_ao, basis at points ====
      call xc_eval_basis(basis, 0, dbl_mb(ichi_ao),
     &        0d0, 0d0, 0d0, dbl_mb(irq),
     &        dbl_mb(igr_xyz), dbl_mb(iatomcd), npoint, natom,
     &        int_mb(iniz), log_mb(idocset),
     &        int_mb(ibas_cent_info), int_mb(ibas_cset_info))

c     ==== Calculating mocf ====
      if (.not. ma_push_get(mt_dbl,nbas,'mocf',lmocf,kmocf))
     $     call errquit('movecs_anal: ma 1 failed', nbas, MA_ERR)

c     ==== Allocating mo_pt ====
      if (.not.ma_push_get(mt_dbl, npoint*nbas,
     &     'mo_pt', lmo_pt, imo_pt))
     &     call errquit('imam: cannot allocate mo_pt',0, MA_ERR)

      write(luout,'("IMAM : Checking MOs at points")')
c     == Loop over orbitals ==
      do imo = 1, nbas
         call ga_get(g_movecs, 1, nbas, imo, imo, dbl_mb(kmocf), nbas)
c        == Loop over points ==
         do ipt = 1, npoint
            mox = 0.0d0
c           == Loop over bases ==
            do ibf = 1, nbas
               mox = mox +
     &               dbl_mb(kmocf + ibf-1) *
     &               dbl_mb(ichi_ao + npoint*(ibf-1) + ipt-1)
            enddo
            
            dbl_mb(imo_pt + npoint*(imo-1) + ipt-1) = mox
         enddo

         i1 = imo_pt + npoint*(imo-1)
         i2 = imo_pt + npoint*(imo-1) + npoint-1
         write(npt_str, '(i3)') npoint
         fr = '(i3, '//npt_str//'(2x, e11.5))'
c         write(fr,'(5a, i3, 15a)') '(i3, ', npoint, '(2x, e11.5))'
c         print *, npoint, fr
         if (ndid.eq.0) write(LuOut, FMT=fr) imo, dbl_mb(i1 : i2)
      enddo
      

c================ Deallocation ===============
      if (.not. MA_pop_stack(l_s))
     &     call errquit('imam: 1e scratch deallocation failed',
     &     60, GA_ERR)
      if (.not. MA_pop_stack(l_pe))
     &     call errquit('imam: pe buffer deallocation failed',
     &     61, GA_ERR)
      if (.not. MA_pop_stack(l_hy))
     &     call errquit('imam: vhyd buffer deallocation failed',
     &     62, GA_ERR)
      if (.not. MA_pop_stack(lgr_xyz))
     $     call errquit("imam: cannot deallocate grid", 0, MA_ERR)
      if (.not. MA_pop_stack(latomcd))
     &     call errquit('imam: cannot deallocate xyz_atom',0, MA_ERR)
      if (.not. MA_pop_stack(latomch))
     &     call errquit('imam: cannot deallocate charge',0, MA_ERR)
      if (.not. MA_pop_stack(ltags))
     &     call errquit('imam: cannot deallocate tags',0, MA_ERR)
      if(.not. MA_pop_stack(lrq))
     &     call errquit("imam: cannot deallocate rq",0, MA_ERR)
      if(.not. MA_pop_stack(lniz))
     &     call errquit("imam: cannot deallocate niz",0, MA_ERR)
      if (.not. MA_pop_stack(ldocset))
     &     call errquit('imam: cannot deallocate ccdocset',0, MA_ERR)
      if (.not. MA_pop_stack(lbas_cent_info))
     &     call errquit("imam: cannot deallocate bas_cent_info",0,
     &     MA_ERR)
      if (.not. MA_pop_stack(lbas_cset_info))
     &     call errquit("imam: cannot deallocate bas_cset_info",0,
     &     MA_ERR)
      if (.not. MA_pop_stack(lchi_ao))
     &     call errquit('imam: cannot deallocate chi_ao',0, MA_ERR)
      
      end
