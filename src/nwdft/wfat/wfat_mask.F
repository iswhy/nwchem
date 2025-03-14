c=====================================================
c=====================================================
      subroutine wfat_get_masks(rtdb, erg, ismask_i, mask_param, 
     &     ismask_e, ergmask, amsk, bmsk)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      double precision erg(*)
      logical ismask_i
      double precision mask_param(2)
      logical ismask_e
      double precision ergmask
      double precision amsk
      double precision bmsk

      logical using_mask
      double precision abmask(2)

      if (.not. rtdb_get(rtdb, 'wfat:mask', mt_log, 1, using_mask))
     &     call errquit('WFAT: Failed reading mask from RTDB', 1, 
     &     RTDB_ERR)
      
      ismask_i = .false.
      if (using_mask .and. rtdb_get(rtdb, 'wfat:mask_g', mt_dbl, 2,
     &     mask_param)) then
         ismask_i = .true.
      else
         mask_param = (/ -1.0d6, -1.0d6 /)
      endif

      ismask_e = .false.
      if (using_mask .and. rtdb_get(rtdb, 'wfat:mask_gjk', mt_dbl, 2,
     &     abmask)) then
         ismask_e = .true.
         ergmask = erg(1)
         amsk = abmask(1)
         bmsk = abmask(2)
      else
         ergmask = 0.0d0
         amsk = -1.0d6
         bmsk = -1.0d6
      endif

      end
c=====================================================
      

c=======================================================
c=======================================================
      subroutine wfat_mask_exp(nq, vpt, orig, zi, kp, a, b, rf, mask)

      implicit none
#include "stdio.fh"
      integer nq
      double precision vpt(3,nq)
      double precision orig(3)
      double precision zi
      double precision kp
      double precision a
      double precision b
      double precision rf(nq)
      double precision mask(nq)

      integer i
      double precision p, x, dv(3), r, a_, b_

      a_ = a / 2
      b_ = b + a/2
      do i = 1, nq
         dv = vpt(:,i) - orig
         r = sqrt(sum(dv*dv))
         x = kp * r
         if (abs(rf(i)) > 1.0d-20) then
            p = 1 / rf(i) * (2*x)**(zi/kp-1) * exp(-x)
            mask(i) = p - (p-1) / 2 * erfc((r-b_)/a_)
         else
            mask(i) = 0.0d0
         endif
cc2         if ( isnan(mask(i)) )
cc2     &        write(luout, '("maskf : ", 4es21.12)') r, x, p, rf(i)
      enddo

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_mask_intd(geom, natom, nch, nqn, nsh, v_m,
     &     mask, intd)

      implicit none
      integer geom
      integer natom                           ! [input], the total number of atoms that corresponds to the geometry handle 'geom'.
      integer nch
      integer nqn
      integer nsh
      double precision v_m(3,nch)             ! [input], the Cartesian coordinate of the evaluation point in molecular-frame.
      double precision mask(2)
      double complex intd(nsh, nqn, nch)

      logical t1
      integer i, j, imf
      double precision v_c(3), r(natom), r_, f, a, b
      parameter( imf=0 )

      a = mask(1) / 2
      b = mask(2) + mask(1)/2
      do i = 1, nch
         r = 0.0d0
         do j = 1, natom
            call geom_cent2cent(geom, imf, v_m(:,i), j, v_c, t1)
            r(j) = sqrt( sum(v_c**2) )
         enddo
         r_ = minval(r)         
cOLD         f = 0.5d0 * erfc(mask(1) * (r_-mask(2)))
         f = 0.5d0 * erfc((r_-b)/a)

         intd(:,:,i) = intd(:,:,i) * f
      enddo

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_mask_intdj(geom, natom, nch, nqn, nsh,
     &     v_m, mask, intd_j)

      implicit none
      integer geom
      integer natom                           ! [input], the total number of atoms that corresponds to the geometry handle 'geom'.
      integer nch
      integer nqn
      integer nsh
      double precision v_m(3,nch)             ! [input], the Cartesian coordinate of the evaluation point in molecular-frame.
      double precision mask(2)
      double complex intd_j(nsh, nqn, 3, nch)

      integer i, j
      double precision v(3,3)

      do i = 1, nch
         do j = 1 ,3
            v(:,j) = v_m(:,i)
         enddo
         call wfat_mask_intd(geom, natom, 3, nqn, nsh, v, mask,
     &        intd_j(:,:,:,i))
      enddo
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_collect_maskf(geom, natom, nch, cgrid, nptp, 
     &     nptop, chpos, v_m, v_oc, mask_param, maskf, maskfo)

      implicit none
#include "errquit.fh"
      integer geom
      integer natom
      integer nch
      logical cgrid(nch)
      integer nptp
      integer nptop
      double precision chpos(3,nch)
      double precision v_m(3,nptp)      ! [input], Cartesian coordinate of the nucleus-grid evaluation points in molecular frame.
      double precision v_oc(3,nptop)    ! [input], Cartesian coordinate of the channel-origin-grid evaluation points in channel-origin frame.
      double precision mask_param(2)
      double precision maskf(nptp)
      double precision maskfo(nptop,nch)

      logical t
      integer i, j
      double complex maskf_, v_t(3)


      call errquit('WFAT: Make sure all instances of call '//
     &     'wfat_collect_maskf have used the new list of '//
     &     'arguments', 54, UERR)

c     ==== Evaluating mask function at nucleus grids ====
      do i = 1, nptp
         maskf_ = cmplx(1.0d0, 0.0d0, kind=kind(1.0d0))
         call wfat_mask_intd(geom, natom, 1, 1, 1, v_m(:,i), 
     &        mask_param, maskf_)
         maskf(i) = real(maskf_)
      enddo

c     ==== Evaluating mask function at channel-origin grids ====
      do j = 1, nch
         if (cgrid(j)) then
            do i = 1, nptop
               call geom_arb2cent(-1000, chpos(:,j), v_oc(:,i), 0, v_t,
     &              t)          ! 1)
c              v_t = Channel-origin grid point in molecular frame.
c              1) The argument geom is supplied with an arbitrary integer 1000 because
c              geom_arb2cent will not verify it if icent2=0, which it is in the call above.
               maskf_ = cmplx(1.0d0, 0.0d0, kind=kind(1.0d0))
               call wfat_mask_intd(geom, natom, 1, 1, 1, v_t, 
     &              mask_param, maskf_)
               maskfo(i,j) = real(maskf_)
            enddo
         else
            maskfo(:,j) = 0.0d0
         endif
      enddo

      end
c=======================================================
    
