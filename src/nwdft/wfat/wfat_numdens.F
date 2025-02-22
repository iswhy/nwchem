c=======================================================
c=======================================================
      subroutine wfat_numdens(bashan, geomhan, nsp, gq, numden)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "bas.fh"
#include "errquit.fh"
#include "geom.fh"
#include "wfat.fh"
#include "global.fh"      
      integer bashan
      integer geomhan
      integer nsp
      integer gq
      double precision numden(1:*)

      character*16 tag
      character*13 spinname
      logical t1
      integer natom, natomp, nbas, nrd, max_npt, npt, nptp
      integer nb0, nb1, nb2, nb3, nd0, nd1, nd2, nbasden
      integer pb0, pb1, pb2, pb3, pd0, pd1, pd2
      integer icent, ic, ipt, jpt, i, ipt_i, ipt_f
      integer l_rd, l_pt, l_wt, l_sc, l_icp, l_ipp
      integer k_rd, k_pt, k_wt, k_sc, k_icp, k_ipp
      integer ipt0, imf, ndid, g_de
      double precision atp(3), znuc, v_c(3), v_m(3), celf
      
      parameter (imf = 0)


      call ga_sync()
      ndid = ga_nodeid()
      

c     ==== Get number of atoms ====
      if (.not. geom_ncent(geomhan, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
c     ==== Get number of bases ====
      if (.not. bas_numbf(bashan, nbas))
     &     call errquit('WFAT: Getting nbas failed in ', 53, UERR)

      
c     ==== Allocate grid-related arrays ====
      call wfat_get_gridmaxpt(gq, max_npt)
c     == Radial points ==
      if (.not. MA_alloc_get(MT_DBL, 3*max_npt, 'gpt', l_pt, k_pt))
     &     call errquit('WFAT: gpt allocation failed', 59, MA_ERR)
c     == Angular points ==
      if (.not. MA_alloc_get(MT_DBL, max_npt, 'gwt', l_wt, k_wt))
     &     call errquit('WFAT: gwt allocation failed', 59, MA_ERR)

      
c     ==== Get the reduced density matrix ====
      nrd = nbas * nbas * nsp
      if (.not. MA_alloc_get(MT_DBL, nrd, 'rdm', l_rd, k_rd))
     &     call errquit('WFAT: rdm allocation failed', 59, MA_ERR)
      call wfat_get_rdm(dbl_mb(k_rd:k_rd+nrd-1))

      
c     ==== Allocate scratch array for basis and densities ====
c     == basis ==
      nb0 = nbas
      nb1 =  3 * nbas
      nb2 =  6 * nbas
      nb3 = 10 * nbas
c     == density ==      
      nd0 = nsp
      nd1 = 3 * nsp
      nd2 = 6 * nsp
      nbasden = nb0 + nb1 + nb2 + nb3 + nd0 + nd1 + nd2
      if (.not. MA_alloc_get(MT_DBL, nbasden, 'scr', l_sc, k_sc))
     &     call errquit('WFAT: scr allocation failed', 1, MA_ERR)
c     == basis ==
      pb0 = k_sc
      pb1 = pb0 + nb0
      pb2 = pb1 + nb1
      pb3 = pb2 + nb2
c     == density ==
      pd0 = pb3 + nb3
      pd1 = pd0 + nd0
      pd2 = pd1 + nd1


      if (.not. ga_create(mt_dbl, nsp, 1, 'dens', 0, 0, g_de))
     &     call errquit('WFAT: Cannot create dens global array',
     &     1, GA_ERR)
      call ga_zero(g_de)
      

c     ==== Obtaining grid ranges for the current process ====
      if (.not. ma_alloc_get(MT_INT, natom, 'icent_p', l_icp, k_icp))
     &     call errquit('WFAT: icent_p allocation failed.', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_INT, 2*natom, 'ipt_p', l_ipp, k_ipp))
     &     call errquit('WFAT: ipt_p allocation failed.', 1, MA_ERR)
      call wfat_get_gridrange(geomhan, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natomp, nptp)
      

c     ==== BEGIN QUADRATURE SUMMATION ====
c     ==== Loop over atoms ====
      jpt = 0
      numden(1:nsp) = 0.0d0
      do ic = 1, natomp
         icent = int_mb(k_icp + ic - 1)
         if (.not. geom_cent_get(geomhan, icent, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)

         call wfat_get_gridpt(int(znuc), gq, max_npt, npt,
     &                        dbl_mb(k_pt : k_pt+3*max_npt-1),
     &                        dbl_mb(k_wt : k_wt+max_npt-1))

         ipt_i = int_mb(k_ipp + (ic-1)*2)
         ipt_f = int_mb(k_ipp + 1 + (ic-1)*2)
c        ==== Loop over all quadrature points around the current atom ====
         do ipt = ipt_i, ipt_f
            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            call geom_cent2cent(geomhan, icent, v_c, imf, v_m, t1)

c           ==== Calculate basis at v_m ====
            call wfat_eval_basis_1pt(geomhan, bashan, 0, v_m,
     &           dbl_mb(pb0 : pb0+nb0-1),
     &           dbl_mb(pb1 : pb1+nb1-1),
     &           dbl_mb(pb2 : pb2+nb2-1),
     &           dbl_mb(pb3 : pb3+nb3-1))
            
            
c           ==== Calculate density at v_m ====
            call wfat_dens_1ptbas(nbas, nsp,
     &           dbl_mb(k_rd : k_rd+nrd-1), 0,
     &           dbl_mb(pb0 : pb0+nb0-1),
     &           dbl_mb(pb1 : pb1+nb1-1),
     &           dbl_mb(pb2 : pb2+nb2-1),
     &           dbl_mb(pd0 : pd0+nd0-1),
     &           dbl_mb(pd1 : pd1+nd1-1),
     &           dbl_mb(pd2 : pd2+nd2-1))
            
c           ==== Quadrature summation ====
            call wfat_celf_1pt(geomhan, icent, v_c, wf_msm, celf)
            if (abs(celf) > 1.0d-14) then
               numden(1:nsp) = numden(1:nsp) +
     &                         dbl_mb(pd0 : pd0+nd0-1) *
     &                         celf * dbl_mb(k_wt+ipt-1)
            endif

cc            jpt = jpt + 1
cc            write(luout, '(i7, " Coord. = ", 3e14.5)')   jpt, v_m
cc            write(luout, '("         sum density = ", e12.5)')
cc     &           sum(dbl_mb(pd0:pd0+nd0-1))
cc            write(luout, '("         sum basis = ", e12.5)')
cc     &           sum(dbl_mb(pb0:pb0+nb0-1))
cc            write(luout, '("         weight = ", e12.5)')
cc     &            dbl_mb(k_wt+ipt-1)
cc            write(luout, '("         cell = ", e12.5)')    celf
         enddo
      enddo


c     ==== Accumulate ====
      call ga_acc(g_de, 1, nsp, 1, 1, numden(1:nsp), nsp, 1.0d0)
      call ga_sync()
      call ga_get(g_de, 1, nsp, 1, 1, numden(1:nsp), nsp)

      if (ndid == 0) then
         write(luout, '(/, '//
     &        '3x, "*** INFO ***", /,'//
     &        '5x, "Numerically integrated density using '//
     &        'the chosen grid quality level in WFAT :")')
         do i = 1, nsp
            if (i == 1) spinname = 'alpha density'
            if (i == 2) spinname = 'beta density'
            if (i == 1 .and. nsp == 1) spinname = 'total density'
            write(luout, '(7x, a13, " = ", f12.6)') spinname, numden(i)
         enddo
      endif


      call ga_sync()
      

c     ==== Deallocations ====
      if (.not. ga_destroy(g_de))
     &     call errquit('WFAT: Cannot destroy dens.', 1, GA_ERR)
      if (.not. MA_free_heap(l_icp))
     &     call errquit('WFAT: icent_p deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ipp))
     &     call errquit('WFAT: ipt_p deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_pt))
     &     call errquit('WFAT: gpt deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_wt))
     &     call errquit('WFAT: gwt deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_rd))
     &     call errquit('WFAT: rdm deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)

      end
c=======================================================


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_numdens(bashan, geomhan, nsp, gq, numden)
cOLD
cOLD      implicit none
cOLD#include "stdio.fh"      
cOLD#include "mafdecls.fh"
cOLD#include "bas.fh"
cOLD#include "errquit.fh"
cOLD#include "geom.fh"
cOLD#include "wfat.fh"
cOLD      integer bashan
cOLD      integer geomhan
cOLD      integer nsp
cOLD      integer gq
cOLD      double precision numden(1:*)
cOLD
cOLD      character*16 tag
cOLD      character*13 spinname
cOLD      logical t1
cOLD      integer natom, nbas, nrd, max_npt, npt
cOLD      integer nb0, nb1, nb2, nb3, nd0, nd1, nd2, nbasden
cOLD      integer pb0, pb1, pb2, pb3, pd0, pd1, pd2
cOLD      integer icent, ipt, jpt, i
cOLD      integer l_rd, l_pt, l_wt, l_sc
cOLD      integer k_rd, k_pt, k_wt, k_sc
cOLD      integer ipt0, imf
cOLD      double precision atp(3), znuc, v_c(3), v_m(3), celf
cOLD      
cOLD      parameter (imf = 0)
cOLD
cOLD
cOLD
cOLDc     ==== Get number of atoms ====
cOLD      if (.not. geom_ncent(geomhan, natom))
cOLD     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
cOLDc     ==== Get number of bases ====
cOLD      if (.not. bas_numbf(bashan, nbas))
cOLD     &     call errquit('WFAT: Getting nbas failed in ', 53, UERR)
cOLD
cOLD      
cOLDc     ==== Allocate grid-related arrays ====
cOLD      call wfat_get_gridmaxpt(gq, max_npt)
cOLDc     == Radial points ==
cOLD      if (.not. MA_alloc_get(MT_DBL, 3*max_npt, 'gpt', l_pt, k_pt))
cOLD     &     call errquit('WFAT: gpt allocation failed', 59, MA_ERR)
cOLDc     == Angular points ==
cOLD      if (.not. MA_alloc_get(MT_DBL, max_npt, 'gwt', l_wt, k_wt))
cOLD     &     call errquit('WFAT: gwt allocation failed', 59, MA_ERR)
cOLD
cOLD      
cOLDc     ==== Get the reduced density matrix ====
cOLD      nrd = nbas * nbas * nsp
cOLD      if (.not. MA_alloc_get(MT_DBL, nrd, 'rdm', l_rd, k_rd))
cOLD     &     call errquit('WFAT: rdm allocation failed', 59, MA_ERR)
cOLD      call wfat_get_rdm(dbl_mb(k_rd:k_rd+nrd-1))
cOLD
cOLD      
cOLDc     ==== Allocate scratch array for basis and densities ====
cOLDc     == basis ==
cOLD      nb0 = nbas
cOLD      nb1 =  3 * nbas
cOLD      nb2 =  6 * nbas
cOLD      nb3 = 10 * nbas
cOLDc     == density ==      
cOLD      nd0 = nsp
cOLD      nd1 = 3 * nsp
cOLD      nd2 = 6 * nsp
cOLD      nbasden = nb0 + nb1 + nb2 + nb3 + nd0 + nd1 + nd2
cOLD      if (.not. MA_alloc_get(MT_DBL, nbasden, 'scr', l_sc, k_sc))
cOLD     &     call errquit('WFAT: scr allocation failed', 1, MA_ERR)
cOLDc     == basis ==
cOLD      pb0 = k_sc
cOLD      pb1 = pb0 + nb0
cOLD      pb2 = pb1 + nb1
cOLD      pb3 = pb2 + nb2
cOLDc     == density ==
cOLD      pd0 = pb3 + nb3
cOLD      pd1 = pd0 + nd0
cOLD      pd2 = pd1 + nd1
cOLD
cOLD
cOLDc     ==== BEGIN QUADRATURE SUMMATION ====
cOLDc     ==== Loop over atoms ====
cOLD      jpt = 0
cOLD      numden(1:nsp) = 0.0d0
cOLD      do icent = 1, natom
cOLD         if (.not. geom_cent_get(geomhan, icent, tag, atp, znuc))
cOLD     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
cOLD
cOLD         call wfat_get_gridpt(int(znuc), gq, max_npt, npt,
cOLD     &                        dbl_mb(k_pt : k_pt+3*max_npt-1),
cOLD     &                        dbl_mb(k_wt : k_wt+max_npt-1))
cOLD
cOLDc        ==== Loop over all quadrature points around the current atom ====
cOLD         do ipt = 1, npt
cOLD            ipt0 = k_pt + 3*(ipt-1)
cOLD            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
cOLD            call geom_cent2cent(geomhan, icent, v_c, imf, v_m, t1)
cOLD
cOLDc           ==== Calculate basis at v_m ====
cOLD            call wfat_eval_basis_1pt(geomhan, bashan, 0, v_m,
cOLD     &           dbl_mb(pb0 : pb0+nb0-1),
cOLD     &           dbl_mb(pb1 : pb1+nb1-1),
cOLD     &           dbl_mb(pb2 : pb2+nb2-1),
cOLD     &           dbl_mb(pb3 : pb3+nb3-1))
cOLD            
cOLD            
cOLDc           ==== Calculate density at v_m ====
cOLD            call wfat_dens_1ptbas(nbas, nsp,
cOLD     &           dbl_mb(k_rd : k_rd+nrd-1), 0,
cOLD     &           dbl_mb(pb0 : pb0+nb0-1),
cOLD     &           dbl_mb(pb1 : pb1+nb1-1),
cOLD     &           dbl_mb(pb2 : pb2+nb2-1),
cOLD     &           dbl_mb(pd0 : pd0+nd0-1),
cOLD     &           dbl_mb(pd1 : pd1+nd1-1),
cOLD     &           dbl_mb(pd2 : pd2+nd2-1))
cOLD            
cOLDc           ==== Quadrature summation ====
cOLD            call wfat_celf_1pt(geomhan, icent, v_c, wf_msm, celf)
cOLD            if (abs(celf) > 1.0d-14) then
cOLD               numden(1:nsp) = numden(1:nsp) +
cOLD     &                         dbl_mb(pd0 : pd0+nd0-1) *
cOLD     &                         celf * dbl_mb(k_wt+ipt-1)
cOLD            endif
cOLD
cOLDcc            jpt = jpt + 1
cOLDcc            write(luout, '(i7, " Coord. = ", 3e14.5)')   jpt, v_m
cOLDcc            write(luout, '("         sum density = ", e12.5)')
cOLDcc     &           sum(dbl_mb(pd0:pd0+nd0-1))
cOLDcc            write(luout, '("         sum basis = ", e12.5)')
cOLDcc     &           sum(dbl_mb(pb0:pb0+nb0-1))
cOLDcc            write(luout, '("         weight = ", e12.5)')
cOLDcc     &            dbl_mb(k_wt+ipt-1)
cOLDcc            write(luout, '("         cell = ", e12.5)')    celf
cOLD         enddo
cOLD      enddo
cOLD
cOLD      
cOLD      write(luout, '(/, '//
cOLD     &     '3x, "*** INFO ***", /,'//
cOLD     &     '5x, "Numerically integrated density using '//
cOLD     &     'the chosen grid quality level in WFAT :")')
cOLD      do i = 1, nsp
cOLD         if (i == 1) spinname = 'alpha density'
cOLD         if (i == 2) spinname = 'beta density'
cOLD         if (i == 1 .and. nsp == 1) spinname = 'total density'
cOLD         write(luout, '(7x, a13, " = ", f12.6)') spinname, numden(i)
cOLD      enddo
cOLD      
cOLD
cOLDc     ==== Deallocations ====
cOLD      if (.not. MA_free_heap(l_pt))
cOLD     &     call errquit('WFAT: gpt deallocation failed', 61, MA_ERR)
cOLD      if (.not. MA_free_heap(l_wt))
cOLD     &     call errquit('WFAT: gwt deallocation failed', 61, MA_ERR)
cOLD      if (.not. MA_free_heap(l_rd))
cOLD     &     call errquit('WFAT: rdm deallocation failed', 61, MA_ERR)
cOLD      if (.not. MA_free_heap(l_sc))
cOLD     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)
cOLD
cOLD      end
cOLDc=======================================================      
      
