c=======================================================
c=======================================================
c==== Imam ====
      subroutine wfat_eval_basis_radf(ghan, bhan, maxder, vpt,
     &     radf, bval, d1bval, d2bval, d3bval)

      implicit none
#include "wfat_qm_param.fh"      
#include "stdio.fh"      
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "errquit.fh"      
      integer ghan                  ! [input] geometry handle
      integer bhan                  ! [input] basis handle
      integer maxder                ! [input] highest derivative order requested
      double precision vpt(3)       ! [input] Cartesian coordinate of the evaluation point
      double precision radf(1:*)    ! [output] the value of the radial function of each basis at vpt
      double precision bval(1:*)    ! [output] the value of basis at vpt
      double precision d1bval(3,*)  ! [output] the value of the 1st derivs of basis at vpt
      double precision d2bval(6,*)  ! [output] the value of the 2nd derivs of basis at vpt
      double precision d3bval(10,*) ! [output] the value of the 3rd derivs of basis at vpt
      
      integer ncont, natom, nbas, npoint
      integer lds, lac, laq, ltag, lrq, lnz, lbcent_info, lbcset_info
      integer ids, iac, iaq, itag, irq, inz, ibcent_info, ibcset_info
      integer ic
      double precision dx, dy, dz

      logical slas
      double precision rsl
      
      parameter (npoint = 1)


c     ==== Filling docset ====
      if (.not. bas_numcont(bhan, ncont))
     &     call errquit('WFAT: Getting ncont failed', 32, BASIS_ERR)
      if (.not. MA_alloc_get(mt_log, ncont, 'docset', lds, ids))
     &     call errquit('WFAT: Cannot allocate docset', ncont, MA_ERR)
      log_mb(ids : ids+ncont-1) = .true.
      
c     ==== Get no. of atoms ====
      if (.not. geom_ncent(ghan, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)

c     ==== Get atomic coordinates and charges ====
      if (.not. MA_alloc_get(mt_dbl, 3*natom, 'a_coord', lac, iac))
     &     call errquit('WFAT: cannot allocate atom coord.',0, MA_ERR)
      if (.not. MA_alloc_get(mt_dbl, natom, 'a_charge', laq, iaq))
     &     call errquit('WFAT: cannot allocate atom charge',0, MA_ERR)
      if (.not. MA_alloc_get(mt_Byte, natom*16, 'a_tag', ltag, itag))
     &     call errquit('WFAT: cannot allocate atom tag',0, MA_ERR)
      if (.not. geom_cart_get(ghan, natom, byte_mb(itag),
     &     dbl_mb(iac), dbl_mb(iaq)))
     &     call errquit('WFAT: Fail to get atom coord and charges',
     &                  0, GEOM_ERR)

      
c     ==== Get the square of point-atom distances ====
      if(.not. MA_alloc_get(mt_dbl, natom, 'rq', lrq, irq))
     &     call errquit("WFAT: Cannot allocate rq",0, MA_ERR)
      do ic = 1, natom
         dx = dbl_mb(iac+3*(ic-1)+0) - vpt(1)
         dy = dbl_mb(iac+3*(ic-1)+1) - vpt(2)
         dz = dbl_mb(iac+3*(ic-1)+2) - vpt(3)
         dbl_mb(irq + npoint*(ic-1)) = dx*dx + dy*dy + dz*dz
      enddo      

c     ==== Filling niz ====
      if(.not. MA_alloc_get(MT_int, natom, 'niz', lnz, inz))
     &     call errquit("WFAT: Cannot allocate niz",0, MA_ERR)
      int_mb(inz : inz+natom-1) = 1

c     ==== Filling bas_cent_info and bas_cset_info ====
      if (.not. MA_alloc_get(mt_int, 3*natom, "bas_cent_info",
     &     lbcent_info, ibcent_info))
     &     call errquit("WFAT: cannot allocate bas_cent_info",0,
     &     MA_ERR)
      if (.not. MA_alloc_get(mt_int, 6*ncont, 'bas_cset_info',
     &     lbcset_info, ibcset_info))
     &     call errquit("WFAT: cannot allocate bas_cset_info",0,
     &     MA_ERR)
      call xc_make_basis_info(bhan, int_mb(ibcent_info),
     &     int_mb(ibcset_info), natom)      

      
c     ==== Calculating basis at the input point ====
      if (.not. bas_numbf(bhan, nbas))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)
      call wfat_eval_basis(bhan, maxder, radf(1:nbas),
     &                   bval(1:nbas),
     &                   d1bval(1:3, 1:nbas),
     &                   d2bval(1:6, 1:nbas),
     &                   d3bval(1:10, 1:nbas),      
     &                   dbl_mb(irq),
     &                   vpt, dbl_mb(iac), npoint, natom,
     &                   int_mb(inz), log_mb(ids),
     &                   int_mb(ibcent_info),
     &                   int_mb(ibcset_info))

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(lds))
     &     call errquit('WFAT: Cannot deallocate docset', ncont, MA_ERR)
      if (.not. MA_free_heap(lac))
     &     call errquit('WFAT: cannot deallocate atom coord.',0, MA_ERR)
      if (.not. MA_free_heap(laq))
     &     call errquit('WFAT: cannot deallocate atom charge',0, MA_ERR)
      if (.not. MA_free_heap(ltag))
     &     call errquit('WFAT: cannot deallocate atom tag',0, MA_ERR)
      if (.not. MA_free_heap(lrq))
     &     call errquit("WFAT: cannot deallocate rq",0, MA_ERR)
      if (.not. MA_free_heap(lnz))
     &     call errquit("WFAT: cannot deallocate niz",0, MA_ERR)
      if (.not. MA_free_heap(lbcent_info))
     &     call errquit("WFAT: cannot deallocate bcent_info",0, MA_ERR)
      if (.not. MA_free_heap(lbcset_info))
     &     call errquit("WFAT: cannot deallocate bcset_info",0, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
c==== Imam ====
      subroutine wfat_eval_basis_1pt(ghan, bhan, maxder, vpt,
     &                bval, d1bval, d2bval, d3bval)

      implicit none
#include "wfat_qm_param.fh"      
#include "stdio.fh"      
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "errquit.fh"      
      integer ghan                  ! [input] geometry handle
      integer bhan                  ! [input] basis handle
      integer maxder                ! [input] highest derivative order requested
      double precision vpt(3)       ! [input] Cartesian coordinate of the evaluation point
      double precision bval(1:*)    ! [output] the value of basis at vpt
      double precision d1bval(3,*)  ! [output] the value of the 1st derivs of basis at vpt
      double precision d2bval(6,*)  ! [output] the value of the 2nd derivs of basis at vpt
      double precision d3bval(10,*) ! [output] the value of the 3rd derivs of basis at vpt
      
      integer ncont, natom, nbas, npoint
      integer lds, lac, laq, ltag, lrq, lnz, lbcent_info, lbcset_info
      integer ids, iac, iaq, itag, irq, inz, ibcent_info, ibcset_info
      integer ic
      double precision dx, dy, dz

      logical slas
      double precision rsl, kp
      
      parameter (npoint = 1)


c     ==== Filling docset ====
      if (.not. bas_numcont(bhan, ncont))
     &     call errquit('WFAT: Getting ncont failed', 9, BASIS_ERR)
      if (.not. MA_alloc_get(mt_log, ncont, 'docset', lds, ids))
     &     call errquit('WFAT: Cannot allocate docset', ncont, MA_ERR)
      log_mb(ids : ids+ncont-1) = .true.
      
c     ==== Get no. of atoms ====
      if (.not. geom_ncent(ghan, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)

c     ==== Get atomic coordinates and charges ====
      if (.not. MA_alloc_get(mt_dbl, 3*natom, 'a_coord', lac, iac))
     &     call errquit('WFAT: cannot allocate atom coord.',0, MA_ERR)
      if (.not. MA_alloc_get(mt_dbl, natom, 'a_charge', laq, iaq))
     &     call errquit('WFAT: cannot allocate atom charge',0, MA_ERR)
      if (.not. MA_alloc_get(mt_Byte, natom*16, 'a_tag', ltag, itag))
     &     call errquit('WFAT: cannot allocate atom tag',0, MA_ERR)
      if (.not. geom_cart_get(ghan, natom, byte_mb(itag),
     &     dbl_mb(iac), dbl_mb(iaq)))
     &     call errquit('WFAT: Fail to get atom coord and charges',
     &                  0, GEOM_ERR)

      
c     ==== Get the square of point-atom distances ====
      if(.not. MA_alloc_get(mt_dbl, natom, 'rq', lrq, irq))
     &     call errquit("WFAT: Cannot allocate rq",0, MA_ERR)
      do ic = 1, natom
         dx = dbl_mb(iac+3*(ic-1)+0) - vpt(1)
         dy = dbl_mb(iac+3*(ic-1)+1) - vpt(2)
         dz = dbl_mb(iac+3*(ic-1)+2) - vpt(3)
         dbl_mb(irq + npoint*(ic-1)) = dx*dx + dy*dy + dz*dz
      enddo      

c     ==== Filling niz ====
      if(.not. MA_alloc_get(MT_int, natom, 'niz', lnz, inz))
     &     call errquit("WFAT: Cannot allocate niz",0, MA_ERR)
      int_mb(inz : inz+natom-1) = 1

c     ==== Filling bas_cent_info and bas_cset_info ====
      if (.not. MA_alloc_get(mt_int, 3*natom, "bas_cent_info",
     &     lbcent_info, ibcent_info))
     &     call errquit("WFAT: cannot allocate bas_cent_info",0,
     &     MA_ERR)
      if (.not. MA_alloc_get(mt_int, 6*ncont, 'bas_cset_info',
     &     lbcset_info, ibcset_info))
     &     call errquit("WFAT: cannot allocate bas_cset_info",0,
     &     MA_ERR)
      call xc_make_basis_info(bhan, int_mb(ibcent_info),
     &     int_mb(ibcset_info), natom)      

      
c     ==== Calculating basis at the input point ====
      if (.not. bas_numbf(bhan, nbas))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)
      call xc_eval_basis(bhan, maxder, 
     &                   bval(1:nbas),
     &                   d1bval(1:3, 1:nbas),
     &                   d2bval(1:6, 1:nbas),
     &                   d3bval(1:10, 1:nbas),      
     &                   dbl_mb(irq),
     &                   vpt, dbl_mb(iac), npoint, natom,
     &                   int_mb(inz), log_mb(ids),
     &                   int_mb(ibcent_info),
     &                   int_mb(ibcset_info))
      

cdebug      slas = .false.
cdebug      rsl = 5.0d0
cdebug      kp = sqrt(2*abs(dbl_mb(k_moerg)))
cdebug      call wfat_eval_basis_slater(bhan, maxder, 
cdebug     &                   bval(1:nbas),
cdebug     &                   d1bval(1:3, 1:nbas),
cdebug     &                   d2bval(1:6, 1:nbas),
cdebug     &                   d3bval(1:10, 1:nbas),      
cdebug     &                   dbl_mb(irq),
cdebug     &                   vpt, dbl_mb(iac), npoint, natom,
cdebug     &                   int_mb(inz), log_mb(ids),
cdebug     &                   int_mb(ibcent_info),
cdebug     &                   int_mb(ibcset_info), slas, kp, rsl)

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(lds))
     &     call errquit('WFAT: Cannot deallocate docset', ncont, MA_ERR)
      if (.not. MA_free_heap(lac))
     &     call errquit('WFAT: cannot deallocate atom coord.',0, MA_ERR)
      if (.not. MA_free_heap(laq))
     &     call errquit('WFAT: cannot deallocate atom charge',0, MA_ERR)
      if (.not. MA_free_heap(ltag))
     &     call errquit('WFAT: cannot deallocate atom tag',0, MA_ERR)
      if (.not. MA_free_heap(lrq))
     &     call errquit("WFAT: cannot deallocate rq",0, MA_ERR)
      if (.not. MA_free_heap(lnz))
     &     call errquit("WFAT: cannot deallocate niz",0, MA_ERR)
      if (.not. MA_free_heap(lbcent_info))
     &     call errquit("WFAT: cannot deallocate bcent_info",0, MA_ERR)
      if (.not. MA_free_heap(lbcset_info))
     &     call errquit("WFAT: cannot deallocate bcset_info",0, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
c==== Imam ====
      subroutine wfat_eval_basis_npt(ghan, bhan, maxder, npt, vpt,
     &                bval, d1bval, d2bval, d3bval)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "errquit.fh"      
      integer ghan                      ! [input] geometry handle
      integer bhan                      ! [input] basis handle
      integer maxder                    ! [input] highest derivative order requested
      integer npt                       ! [input] no. of points
      double precision vpt(3,npt)       ! [input] Cartesian coordinate of the evaluation point
      double precision bval(npt,*)      ! [output] the value of basis at vpt
      double precision d1bval(npt,3,*)  ! [output] the value of the 1st derivs of basis at vpt
      double precision d2bval(npt,6,*)  ! [output] the value of the 2nd derivs of basis at vpt
      double precision d3bval(npt,10,*) ! [output] the value of the 3rd derivs of basis at vpt

      character*16 atag
      integer ncont, natom, nbas
      integer lds, lac, lrq, lnz, lbcent_info, lbcset_info
      integer ids, iac, irq, inz, ibcent_info, ibcset_info
      integer ic, iac1, iac2, i
      double precision dr(3), acord(3), ach
      


c     ==== Filling docset ====
      if (.not. bas_numcont(bhan, ncont))
     &     call errquit('WFAT: Getting ncont failed', 54, BASIS_ERR)
      if (.not. MA_alloc_get(mt_log, ncont, 'docset', lds, ids))
     &     call errquit('WFAT: Cannot allocate docset', ncont, MA_ERR)
      log_mb(ids : ids+ncont-1) = .true.

      
c     ==== Get no. of atoms ====
      if (.not. geom_ncent(ghan, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)

      
c     ==== Get the square of point-atom distances ====
      if (.not. MA_alloc_get(mt_dbl, 3*natom, 'a_coord', lac, iac))
     &     call errquit('WFAT: cannot allocate atom coord.',0, MA_ERR)
      if (.not. MA_alloc_get(mt_dbl, npt*natom, 'rq', lrq, irq))
     &     call errquit("WFAT: Cannot allocate rq",0, MA_ERR)
      do ic = 1, natom
         if (.not. geom_cent_get(ghan, ic, atag, acord, ach))
     &        call errquit('WFAT: Cannot get center info', 1, GEOM_ERR)
         iac1 = iac + (ic-1)*3
         iac2 = iac1 + 2
         dbl_mb(iac1:iac2) = acord
         do i = 1, npt
            dr = acord - vpt(1:3, i)
            dbl_mb(irq + npt*(ic-1) + i - 1) = sum(dr*dr)
         enddo
      enddo
      

c     ==== Filling niz ====
      if(.not. MA_alloc_get(MT_int, natom, 'niz', lnz, inz))
     &     call errquit("WFAT: Cannot allocate niz",0, MA_ERR)
      int_mb(inz : inz+natom-1) = 1

      
c     ==== Filling bas_cent_info and bas_cset_info ====
      if (.not. MA_alloc_get(mt_int, 3*natom, "bas_cent_info",
     &     lbcent_info, ibcent_info))
     &     call errquit("WFAT: cannot allocate bas_cent_info",0,
     &     MA_ERR)
      if (.not. MA_alloc_get(mt_int, 6*ncont, 'bas_cset_info',
     &     lbcset_info, ibcset_info))
     &     call errquit("WFAT: cannot allocate bas_cset_info",0,
     &     MA_ERR)
      call xc_make_basis_info(bhan, int_mb(ibcent_info),
     &     int_mb(ibcset_info), natom)      

      
c     ==== Calculating basis at the input point ====
      if (.not. bas_numbf(bhan, nbas))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)
      call xc_eval_basis(bhan, maxder, 
     &                   bval(1:npt, 1:nbas),
     &                   d1bval(1:npt, 1:3, 1:nbas),
     &                   d2bval(1:npt, 1:6, 1:nbas),
     &                   d3bval(1:npt, 1:10, 1:nbas),      
     &                   dbl_mb(irq),
     &                   vpt, dbl_mb(iac), npt, natom,
     &                   int_mb(inz), log_mb(ids),
     &                   int_mb(ibcent_info),
     &                   int_mb(ibcset_info))

c      call xc_eval_basis(bhan, maxder, 
c     &                   bval,
c     &                   d1bval,
c     &                   d2bval,
c     &                   d3bval,      
c     &                   dbl_mb(irq),
c     &                   vpt, dbl_mb(iac), npt, natom,
c     &                   int_mb(inz), log_mb(ids),
c     &                   int_mb(ibcent_info),
c     &                   int_mb(ibcset_info))

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(lds))
     &     call errquit('WFAT: Cannot deallocate docset', ncont, MA_ERR)
      if (.not. MA_free_heap(lac))
     &     call errquit('WFAT: cannot deallocate atom coord.',0, MA_ERR)
      if (.not. MA_free_heap(lrq))
     &     call errquit("WFAT: cannot deallocate rq",0, MA_ERR)
      if (.not. MA_free_heap(lnz))
     &     call errquit("WFAT: cannot deallocate niz",0, MA_ERR)
      if (.not. MA_free_heap(lbcent_info))
     &     call errquit("WFAT: cannot deallocate bcent_info",0, MA_ERR)
      if (.not. MA_free_heap(lbcset_info))
     &     call errquit("WFAT: cannot deallocate bcset_info",0, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval_basis(basis_hand, maxder, radf, chi, d1chi, 
     &     d2chi, d3chi, rq, qxyz, ccoord, nq, natoms, iniz, docset,
     &     bas_cent_info, bas_cset_info)

      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
c!!! Will we still need this?
#include "bas.fh"
#include "stdio.fh"  
c
      integer basis_hand        !< [Input] basis set handle
      integer maxder            !< [Input] max derivative order to compute
      integer nq                !< [Input] number of quadrature points
      integer natoms            !< [Input] number of atoms
      double precision radf(nq,*) !<[Output] radial part of the basis functions
      double precision chi(nq,*) !< [output] basis function values
      double precision d1chi(nq,3,*) !< [output] basis function 1st derivs
      double precision d2chi(nq,6,*) !< [output] basis function 2nd derivs
      double precision d3chi(nq,10,*) !< [output] basis function 3rd derivs
      double precision rq(nq,natoms) !< [Input] squared pt-ctr distances
      double precision qxyz(3,nq) !< [Input] grid point coordinates
      double precision ccoord(3,natoms) !< [Input] nuclear coordinates
      integer iniz(*)           !< [Input] iniz(i).eq.0 means the basis 
                                !< functions on atom \f$i\f$ are not evaluated,
                                !< only used in derivatives where atoms
                                !< positions may be frozen
      logical docset(*)         !< [Input] whether cset set is significant
      integer bas_cent_info(3,natoms) !< [Input] Basis set info for each
                                      !< center
                                      !< - 1. The maximum angular momentum
                                      !< - 2. The lowest rank shell on the
                                      !<   atom
                                      !< - 3. The highest rank shell on the
                                      !<   atom
      integer bas_cset_info(6,*) !< [Input] Basis set info for each shell
                                 !< - 1. The lowest rank basis function in the
                                 !<   shell
                                 !< - 2. The highest rank basis function in the
                                 !<   shell
                                 !< - 3. The angular momentum of the shell
                                 !< - 4. The number of primitives in the shell
                                 !< - 5. The number of contractions in the shell
                                 !< - 6. Whether the shell is Cartesian or
                                 !<   spherical harmonic
c     
      external nang_cart, nang_cart_tot
      integer nang_cart, nang_cart_tot
c
c     Local declarations
c
      integer gen_l
      integer MXprim,MXang,MXcont
      integer izprim,iccoef,iccoef1
      integer max_elem_ang_scr,max_der_component,max_ang_for_der
      integer iang,id1ang,id2ang,id3ang,il2k,idrd,irad,itop
      integer iscr,lscr
      integer ictr,icset,ichi,ib
      integer ifirst,ilast,ltype,nprimo,ncontr,nshbf(2)
      integer l,nang,iang_beg,iang_end
      integer icontr,l0(2)
      logical spherical,LResult
cOLD #define Ax 2
cOLD #define Ay 3
cOLD #define Az 4
cOLD #define Axx 5
cOLD #define Axy 6
cOLD #define Axz 7
cOLD #define Ayy 8
cOLD #define Ayz 9
cOLD #define Azz 10
cOLD #define Axxx 11
cOLD #define Axxy 12
cOLD #define Axxz 13
cOLD #define Axyy 14
cOLD #define Axyz 15
cOLD #define Axzz 16
cOLD #define Ayyy 17
cOLD #define Ayyz 18
cOLD #define Ayzz 19
cOLD #define Azzz 20
cOLD c
cOLD #define Ad1 5
cOLD #define Ad2 6
cOLD #define Ad3 7
cOLD #define Ad4 8
cOLD #define Ad5 9
cOLD #define Af1 10
cOLD #define Af2 11
cOLD #define Af3 12
cOLD #define Af4 13
cOLD #define Af5 14
cOLD #define Af6 15
cOLD #define Af7 16
cOLD c
cOLD #define dzero 1
cOLD #define dx 2
cOLD #define dy 3
cOLD #define dz 4
cOLD #define dxx 5
cOLD #define dxy 6
cOLD #define dxz 7
cOLD #define dyy 8
cOLD #define dyz 9
cOLD #define dzz 10
cOLD #define dxxx 11
cOLD #define dxxy 12
cOLD #define dxxz 13
cOLD #define dxyy 14
cOLD #define dxyz 15
cOLD #define dxzz 16
cOLD #define dyyy 17
cOLD #define dyyz 18
cOLD #define dyzz 19
cOLD #define dzzz 20
c
c     Determine when to use general code
c
      spherical = bas_is_spherical(basis_hand)
      if (spherical) then
         if (maxder .eq. 0) then
            gen_l = 4
         elseif (maxder .eq. 1) then
            gen_l = 3
         else
            gen_l = 2
         endif
      else
         if (maxder .le. 1) then
            gen_l = 4
         else
            gen_l = 3
         endif
      endif
c
c     Get the max number of contractions, number of primitives
c     and angular momentum for any given shell of this basis set.
c
      if( .not. bas_ncontr_cn_max(basis_hand,MXcont) )then
         call errquit('Exiting in xc_eval_basis.',1, BASIS_ERR)
      end if
      if( .not. bas_nprim_cn_max(basis_hand,MXprim) )then
         call errquit('Exiting in xc_eval_basis.',2, BASIS_ERR)
      end if
      if( .not. bas_high_angular(basis_hand,MXang) )then
         call errquit('Exiting in xc_eval_basis.',3, BASIS_ERR)
      end if
c
c     Allocate scratch space for primitives for a single shell
c
      
      izprim = 0
      iccoef = izprim + MXprim
      itop = iccoef + MXprim*MXcont
c
c     Allocate scratch space for angular and radial functions
c     Allocate cartesian no matter what since in the general case
c     the cartesian functions are evaluated first and then transformed
c
      max_elem_ang_scr = nang_cart_tot(MXang)
      max_der_component = nang_cart_tot(maxder)
c
c     If basis function derivatives are needed, always evaluate x, y and z
c     angular functions for derivative purposes no matter what the angular
c     momentum is
c
      max_ang_for_der = min(maxder,1)
      max_elem_ang_scr = max(max_elem_ang_scr,
     &                       nang_cart_tot(max_ang_for_der))
c
      iang = itop
      itop = iang + nq*max_elem_ang_scr
      id1ang = iang
      id2ang = iang
      id3ang = iang
      il2k = iang
      if (MXang .ge. gen_l) then
c!!! Could subtract off everything below gen_l here and change the
c!!! dimensioning of the d*ang arrays to save the space if needed
         if (maxder .ge. 1) then
            id1ang = itop
            il2k = id1ang + nq*3*max_elem_ang_scr
            itop = il2k + (MXang+1)*(MXang+1)*(MXang+1)
         endif
         if (maxder .ge. 2) then
            id2ang = itop
            itop = id2ang + nq*6*max_elem_ang_scr
         endif
         if (maxder .ge. 3) then
            id3ang = itop
            itop = id3ang + nq*10*max_elem_ang_scr
         endif
      endif
      idrd = itop
      irad = idrd + nq*max_der_component*MXcont
      itop = irad + nq*(maxder+1)*MXcont
      if (.not. MA_Push_Get(mt_dbl, itop, 'xc_eval_basis', lscr, iscr))
     &     call errquit('Could not allocate scratch in xc_eval_basis',
     &     itop, MA_ERR)
      izprim = iscr
      iccoef = iccoef + iscr
      iang = iang + iscr
      id1ang = id1ang + iscr
      id2ang = id2ang + iscr
      id3ang = id3ang + iscr
      il2k = il2k + iscr
      idrd = idrd + iscr
      irad = irad + iscr
c
c     Offset for basis functions in chi arrays
c
      ichi = 1
c
      do 100 ictr = 1,natoms
         if (iniz(ictr) .eq. 0) goto 100
c
c     Evaluate the angular functions up through the max angular momentum
c     on the current center, so that the angular function values can be
c     for all functions on this center - also make sure that the angular
c     functions needed for derivatives get evaluated
c
         call eval_angular(dbl_mb(iang), dbl_mb(id1ang), dbl_mb(id2ang),
     &        dbl_mb(id3ang), qxyz, ccoord(1,ictr), nq, dbl_mb(il2k),
     &        max(bas_cent_info(1,ictr),max_ang_for_der), maxder, gen_l,
     &        spherical)         
c
c     Loop over all contraction sets on this center
c
         do 200 icset = bas_cent_info(2,ictr),bas_cent_info(3,ictr)
c
c     Check if we have neglected this contraction
c
            if (.not. docset(icset)) goto 200
c
c     Get info about current contraction set - first and last basis function,
c     angular momentum, number of primitives, number of contractions and
c     whether spherical harmomic
c
            ifirst = bas_cset_info(1,icset)   ! The first basis among the m=-l,...,+l functions.
            ilast  = bas_cset_info(2,icset)   ! The last basis among to m=-l,...,+l functions.
            ltype  = bas_cset_info(3,icset)   ! Angular momentum type: 0(s), 1(p), 2(d), etc.
            nprimo = bas_cset_info(4,icset)
            ncontr = bas_cset_info(5,icset)
            spherical = bas_cset_info(6,icset) .ne. 0
cfordebug            write(luout, '("idx: ", 4i5)') ifirst, ilast, ltype, ncontr
cfordebug            if (ltype > 3) call errquit('Stop inside basis', 9, UERR)
c
c     Number of basis functions counting all contractions
c
            nshbf(1) = ilast-ifirst+1
c
c     Angular momentum
c
            if (ltype .lt. 0) then
               l0(1)=0
               l0(2)=1
               nshbf(1)=1
               nshbf(2)=3
            else
               l0(1) = ltype
cold               l = ltype
            endif
c
c     Get exponents and contraction coefficients for this contraction set
c
c!!! These calls could be causing too much overhead - check
            LResult = bas_get_exponent(basis_hand,icset,dbl_mb(izprim))
            LResult = bas_get_coeff(basis_hand,icset,dbl_mb(iccoef))
            do icontr=1,ncontr
               l = l0(icontr)
c
c     Cartesian/spherical harmonic
c
            if (spherical) then
               nang = 2*l+1
               iang_end = (l+1)*(l+1)
               iang_beg = iang_end - nang + 1
            else
               nang = nang_cart(l)
               iang_end = nang_cart_tot(l)
               iang_beg = iang_end - nang + 1
            endif
c
c     Evaluate the radial components of the contracted Gaussian functions
c
            iccoef1=iccoef+(icontr-1)*nprimo
            if(maxder.eq.0) then
               call eval_radial0(dbl_mb(irad), rq(1,ictr), nq, 
     &              dbl_mb(izprim), dbl_mb(iccoef1), nprimo, ncontr)
            else
               call eval_radial(dbl_mb(irad), rq(1,ictr), nq, maxder,
     &              dbl_mb(izprim), dbl_mb(iccoef1), nprimo, ncontr)
            endif            
c
c     Compute and store the derivatives of the radial function
c

c           ==== Apply mask ====
cc            call wfat_mask_exp(nq, qxyz, orig, zi, kp, amsk, bmsk,
cc     &           dbl_mb(irad:irad+nq*(maxder+1)-1:maxder+1))
cc2            write(luout, '("q = ", 3es15.6)') qxyz(:,1)
            do ib = ifirst, ilast
               radf(:,ib) = dbl_mb(irad:irad+nq*(maxder+1)-1:maxder+1)
cc2               write(luout, '("radf : ", i4, es21.12)') ib, radf(1,ib)
            enddo
            

            
            call eval_der_radial(dbl_mb(idrd), dbl_mb(irad),
     &           dbl_mb(iang), nq, maxder, max_der_component)

c     Evaluate basis functions and begin evaluation of basis function
c     derivative by multiplying angular function by the radial function
c     and its derivatives
            call eval_ang_x_drad(chi(1,ichi), d1chi(1,1,ichi),
     &           d2chi(1,1,ichi), d3chi(1,1,ichi), dbl_mb(idrd),
     &           dbl_mb(iang), nq, iang_beg, iang_end,
     &           max_der_component)
c
c     Finish off with terms involving angular function derivatives
c
            if (l .gt. 0 .and. maxder .gt. 0) then
               call eval_dang_x_rad(chi(1,ichi), d1chi(1,1,ichi),
     &              d2chi(1,1,ichi), d3chi(1,1,ichi), dbl_mb(idrd),
     &              dbl_mb(iang), dbl_mb(id1ang), dbl_mb(id2ang),
     &              dbl_mb(id3ang), nq, iang_beg, iang_end,
     &              max_der_component, gen_l)
            endif
c
c     Advance offset into chi arrays by number of basis functions
c     just processed
c
            ichi = ichi + nshbf(icontr)
            enddo
c
 200     continue
 100  continue
c
      if (.not. MA_Pop_Stack(lscr))
     &     call errquit('xc_eval_basis: could not pop stack',0, MA_ERR)
c
      return

      end
c=======================================================
      
      
      
