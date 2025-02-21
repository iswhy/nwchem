c=======================================================
c=======================================================
      subroutine wfat_eval_basis_slater(basis_hand, maxder, chi, d1chi, 
     &     d2chi, d3chi, rq, qxyz, ccoord, nq, natoms, iniz, docset,
     &     bas_cent_info, bas_cset_info, slas, kp, rsl)
c
c     $Id$
c
c     Evaluate basis functions and their derivatives on a grid
c
c     BGJ - 7/00
c
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
      logical slas               !< [Input]
      double precision kp
      double precision rsl
      
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
      integer ictr,icset,ichi
      integer ifirst,ilast,ltype,nprimo,ncontr,nshbf(2)
      integer l,nang,iang_beg,iang_end
      integer icontr,l0(2)
      logical spherical,LResult
      double precision rad_sl
#define Ax 2
#define Ay 3
#define Az 4
#define Axx 5
#define Axy 6
#define Axz 7
#define Ayy 8
#define Ayz 9
#define Azz 10
#define Axxx 11
#define Axxy 12
#define Axxz 13
#define Axyy 14
#define Axyz 15
#define Axzz 16
#define Ayyy 17
#define Ayyz 18
#define Ayzz 19
#define Azzz 20
c
#define Ad1 5
#define Ad2 6
#define Ad3 7
#define Ad4 8
#define Ad5 9
#define Af1 10
#define Af2 11
#define Af3 12
#define Af4 13
#define Af5 14
#define Af6 15
#define Af7 16
c
#define dzero 1
#define dx 2
#define dy 3
#define dz 4
#define dxx 5
#define dxy 6
#define dxz 7
#define dyy 8
#define dyz 9
#define dzz 10
#define dxxx 11
#define dxxy 12
#define dxxz 13
#define dxyy 14
#define dxyz 15
#define dxzz 16
#define dyyy 17
#define dyyz 18
#define dyzz 19
#define dzzz 20
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
            ifirst = bas_cset_info(1,icset)
            ilast  = bas_cset_info(2,icset)
            ltype  = bas_cset_info(3,icset)
            nprimo = bas_cset_info(4,icset)
            ncontr = bas_cset_info(5,icset)
            spherical = bas_cset_info(6,icset) .ne. 0
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

            if (slas) then
               call eval_radial(rad_sl, rsl**2, 1, 0,
     &              dbl_mb(izprim), dbl_mb(iccoef1), nprimo, ncontr)
               call wfat_slater_radial(kp, rsl, rad_sl, nq,
     &              maxder, rq(1,ictr), dbl_mb(irad))
            endif
c
c     Compute and store the derivatives of the radial function
c
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


c=======================================================
c=======================================================
      subroutine wfat_slater_radial(kp, rsl, rad_sl, nq, maxd, r2,
     &     crad)

      implicit none
#include "errquit.fh"
      double precision kp                ! [input],
      double precision rsl               ! [input],
      double precision rad_sl            ! [input],
      integer nq                         ! [input],
      integer maxd                       ! [input],
      double precision r2(nq)            ! [input],
      double precision crad(0:maxd,nq)   ! [output],

      integer i
      double precision amp, r


      if (kp < 0.0d0) call errquit('WFAT: The argument kp of '//
     &     'wfat_slater_radial must be positive or zero', 83, UERR)
      
      do i = 1, nq
         r = sqrt(r2(i))
         if (r > rsl) then
            amp = rad_sl * exp(kp * rsl)
            crad(0,i) = amp * exp(-kp * r)
         endif
      enddo

      end
c=======================================================
      
