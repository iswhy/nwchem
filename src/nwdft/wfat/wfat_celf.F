c=======================================================
c=======================================================
      subroutine wfat_celf_1pt_c(geom, nch, geomc, icent, cgrid, chpos,
     &     v_m, msm, celf)
c     DESCRIPTION:
c     geomc is an array of geometry objects. Each element of geomc consists
c     of atomic coordinates such that the second up to the last atoms and
c     their coordinates are identical to those in geom. The first atom in the
c     elements of geomc is the parabolic channel center of each channel.
      implicit none
      integer geom
      integer nch
      integer geomc(nch)
      integer icent
      logical cgrid(nch)
      double precision chpos(3,nch)
      double precision v_m(3,nch)
      integer msm
      double precision celf(nch)

      logical t
      integer ic, i, k
      double precision v_c(3)


      ic = icent + 1

c     ==== Calculate cell function for a point belonging to a ====
c     ====  nucleus and when the channel does not have center ====
      do i = 1, nch
         if (.not. cgrid(i) .and. icent /= -1) then
            k = i
            call geom_cent2cent(geom, 0, v_m(:,k), icent, v_c, t)
            call wfat_celf_1pt(geom, icent, v_c, msm, celf(k))
            exit
         endif
      enddo

      
c     ==== Cell functions for all channels ====
      if (icent /= -1) then
c        == When the center is a nucleus ==
         do i = 1, nch
            if (.not. cgrid(i)) then
               celf(i) = celf(k)
            else
               call geom_cent2cent(geomc(i), 0, v_m(:,i), ic, v_c, t)
               call wfat_celf_1pt(geomc(i), ic, v_c, msm, celf(i))
            endif
         enddo
      else
c        == When the center is a channel center ==
         do i = 1, nch
            if (.not. cgrid(i)) then
               celf(i) = 0.0d0
            else
               call geom_cent2arb(geomc(i), 0, v_m(:,i), chpos(:,i),
     &              v_c, t)
               call wfat_celf_1pt(geomc(i), 1, v_c, msm, celf(i))
c              The second argument above is 1 because the channel center
c              is the first center of the geometry object geomc(i).
            endif
         enddo
      endif
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_celf_1pt(geom, ia1, d1, msm, celf)
c     DESCRIPTION:
c     This subroutine evaluates the function w_A defined in eq. (23) in
c     Ref. [1].
      
      implicit none
#include "errquit.fh"
#include "geom.fh"
      integer geom, ia1, msm
      double precision d1(3)
      double precision celf

      logical check
      integer ia2, natom
      double precision d2(3), tt, parf, tot_parf

      
      if (.not. geom_ncent(geom, natom))
     &    call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
      if (.not. geom_check_cent(geom, 'wfat_celf_1pt', ia1))
     &    call errquit('WFAT: Argument ia1 of wfat_celf_1pt refers '//
     &                 'to a non-existent center ID in geom.', 51,
     &                 GEOM_ERR)
      
c     ==== Loop over all atoms to be summed to  ====
c     ==== obtain the summed partition function ====
      tot_parf = 0.0d0
      do ia2 = 1, natom
         call geom_cent2cent(geom, ia1, d1, ia2, d2, check)
         
c        ==== Calculate one-atom partition function ====
c        ====   function at the current grid point  ====
         call wfat_parf_1pt(geom, ia2, d2, msm, tt)
         if (ia2 == ia1) parf = tt 

c        ====   Calculate the summed partition   ====
c        ==== function at the current grid point ====
         tot_parf = tot_parf + tt
      enddo

      celf = parf / tot_parf       ! Evaluates eq. (23) of Ref. [1].

      end
      
c  References :
c    [1] "Quadrature schemes for integrals of density functional theory",
c        Mol. Phys., 78 (4) 1993.
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_parf_1pt(geom, ia1, d1, msm, parf)

c     DESCRIPTION:
c       This subroutine evaluates the partition function, parttf,
c       defined in eq. (22) in Ref. [1] at a spatial point d1(1:3).
c       This evaluation point is in the frame of reference
c       of atom ia1.

      implicit none
#include "errquit.fh"
#include "geom.fh"      
      integer geom, ia1, msm
      double precision d1(3)
      double precision parf

      character*16 tg1, tg2
      logical check
      integer ia2, natom
      double precision rr1, va1(3), zn1, atr1
      double precision rr2, va2(3), zn2, atr2, d2(3)
      double precision v12(3), dist
      double precision mu, chi, u, a, nu, smof
      

      if (.not. geom_ncent(geom, natom))
     &    call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
      if (.not. geom_check_cent(geom, 'wfat_parf_1pt', ia1))
     &    call errquit('WFAT: Argument ia1 of wfat_parf_1pt refers '//
     &                 'to a non-existent center ID in geom.', 51,
     &                 GEOM_ERR)

      rr1 = sqrt(sum(d1*d1))
      
c     ==== Get atomic coordinate and nuclear charge for ia1 ====
      if (.not. geom_cent_get(geom, ia1, tg1, va1, zn1))
     &     call errquit('WFAT: Getting zn1 failed', 51, GEOM_ERR)
      call wfat_get_atomrad(int(zn1), atr1)      
      
      parf = 1.0d0
      do ia2 = 1, natom
         
         if (ia2 == ia1) cycle
         
c        ==== Get atomic coordinate and nuclear charge for ia2 ====   
         if (.not. geom_cent_get(geom, ia2, tg2, va2, zn2))
     &        call errquit('WFAT: Getting zn2 failed', 51, GEOM_ERR)
         v12 = va2 - va1
         dist = sqrt(sum(v12*v12))
      
c        ==== Calculate the distance of the input ====
c        ====  Cartesian coordinate from atom ia2 ====
         call geom_cent2cent(geom, ia1, d1, ia2, d2, check)
         rr2 = sqrt(sum(d2*d2))
      
c        ==== Calculate mu, chi, u, a, and nu ====
         mu = (rr1 - rr2)/dist
         call wfat_get_atomrad(int(zn2), atr2)
         chi = atr1 / atr2
         u = (chi-1) / (chi+1)
         a = u / (u**2 - 1)
         nu = mu + a * (1-mu**2)
      
c        ==== Evaluate the mask function ====
         call wfat_eval_smoothf(nu, msm, smof)
         parf = parf * smof
      
      enddo

      end
c  References :
c    [1] "Quadrature schemes for integrals of density functional theory",
c        Mol. Phys., 78 (4) 1993.      
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval_smoothf(x, msm, smof)

c     DESCRIPTION:
c       This subroutine evaluates s(mu) whose first derivative is defined
c       in eq. (24) of Ref. [1] at mu=x.
      
      implicit none
      integer msm
      double precision x
      double precision smof

      integer i
      double precision smof1, gg, am
      double precision xlm_combin
      external xlm_combin
      
      
c     ==== Do calculation ====
      smof1 = 0.0d0
      smof = 0.0d0
      do i = 0, msm           !! 1)
         gg = xlm_combin(msm,i) * (-1)**i / dble(2*i+1)
         smof1 = smof1 + gg      !! 2)
         smof = smof + gg * (x**(2*i+1)+1)
      enddo
      
      am = -1.0d0 / 2 / smof1
      smof = am * smof + 1

      end
c  Notes:
c    1) This loop runs over all terms in the binomial expansion expression of s(mu) defined
c       in eq. (24) of Ref. [1]. Note that the expression of s(mu), which is what this loop
c       calculates, is obtained by expanding eq. (24) and then integrating over mu.
c    2) smof1 is for calculating the constant prefactor A_m in eq. (24) of Ref. [1]. Here,
c       the same binomial expansion as that for s(mu) is used except that mu is substituted
c       with 1.
c
c     References :
c    [1] "Quadrature schemes for integrals of density functional theory",
c        Mol. Phys., 78 (4) 1993.      
c=======================================================
