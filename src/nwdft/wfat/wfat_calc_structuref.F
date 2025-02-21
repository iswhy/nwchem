c=======================================================
c=======================================================
      subroutine wfat_calc_srfc0(nat, nqn, bt, gm, erg, chmmt, gn)

      implicit none
      integer nat, nqn
      double precision bt, gm                    ! [input], orientation angle in degree.
      double precision erg, chmmt(3)                 ! [input], 
      double complex gn(nat, nqn)                    ! [in/output], 

      integer iqn, ic
      double precision protmat(3,3), mu(3), kp
      
      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)

      
c     ==== Begin calculating structure factor ====
      mu = matmul(transpose(protmat), -chmmt)
      kp = sqrt(2 * abs(erg))
      
      do iqn = 1, nqn      
         do ic = 1, nat
            gn(ic,iqn) = exp(-kp * mu(3)) * gn(ic,iqn)
         enddo
      enddo
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_srfc1(nat, nqn, bt, gm, erg, chmmt, jn, kn)

      implicit none
      integer nat, nqn
      double precision bt, gm                    ! [input], orientation angle in degree.
      double precision erg, chmmt(3)                 ! [input], 
      double complex jn(nat, nqn)                    ! [in/output], 
      double complex kn(nat, nqn)                    ! [in/output], 

      integer iqn, ic
      double precision protmat(3,3), mu(3), kp
      
      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)

      
c     ==== Begin calculating structure factor ====
      mu = matmul(transpose(protmat), -chmmt)
      kp = sqrt(2 * abs(erg))
      
      do iqn = 1, nqn      
         do ic = 1, nat
            jn(ic,iqn) = exp(-kp * mu(3)) * jn(ic,iqn)
            kn(ic,iqn) = exp(-kp * mu(3)) * kn(ic,iqn)
         enddo
      enddo
      
      end
c=======================================================

      
c=======================================================
c=======================================================
c==== Imam ====
      subroutine wfat_get_srfc0(geom, bt, gm, nch, erg, nqn, mid, 
     &                          lmax_, chmmt, intgl, nopw, srfc)
      
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "stdio.fh"
      integer geom
      integer nch, nqn, mid(nqn), lmax_
      double precision bt, gm                    ! [input], orientation angle in degree.
      double precision erg(nch), chmmt(3,nch)
      double complex intgl(1:*)    ! 1)
      logical nopw                 ! 3)
      double complex srfc(1:*)     ! 2)
      
      integer ich, iqn, l, m, i_lm, icent, i, j, nat, nsh, l0
      double precision btr, gmr, wgf, mu(3), kappa, protmat(3,3)
      double precision deg2rad
      double complex g_nu0, g_nu, czero, imnum
      parameter ( deg2rad = 0.01745329251994d0 )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )
      
c     3) When nopw=.true., this subroutine simply multiplies the input asymptotic
c        coefficients with the dipole-moment containing exponential factor to yield
c        the structure factor. Also, when nopw=.true., lmax_ must be zero.
      

c     NOTE:
C     Consider replacing the executable portion of this subroutine by a call
c     to wfat_calc_gn followed by a call to wfat_calc_srfc0.
      
      if (nopw .and. lmax_ /= 0)
     &     call errquit('WFAT: If nopw is .true., then lmax_ must be '//
     &     'equal to zero (wfat_get_srfc0).', 43, UERR)
      
      
c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 77, GEOM_ERR)

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      btr = bt * deg2rad
      gmr = gm * deg2rad

      
c     ==== Begin calculating structure factor ====
      nsh = (1 + lmax_)**2

      
      do ich = 1, nch
         mu = matmul(transpose(protmat), -chmmt(:, ich))
         kappa = sqrt(2 * abs(erg(ich)))
         do iqn = 1, nqn
      
            do icent = 1, nat
               g_nu = czero
               if (nopw) then
                  l0 = 0
               else
                  l0 = abs(mid(iqn))
               endif
               
               do l = l0, lmax_
                  do m = -l, l
                     if (nopw) then
                        wgf = 1.0d0
                     else
                        call wfat_wignerf(l, mid(iqn), m, btr, wgf)
                     endif
c check                  call wfat_wignerf(l, mid(iqn), m, -btr, wgf)
c check                  if (abs(btr) < 1.0d-13) write(*, '("wigner:", e15.6, i3, i3, i3, e15.6)') btr, l, mid(iqn), m, wgf
                     
                     i_lm = l**2 + l + m + 1
                     i = nsh*nat*nqn*(ich-1) + nsh*nat*(iqn-1) +
     &                   nsh*(icent-1) + i_lm
                     g_nu0 = intgl(i) * wgf * exp(-imnum * m * gmr)
c check                  g_nu0 = intgl(i_lm, icent, iqn, ich) * wgf * exp(-imnum * m * -gmr)
                     g_nu = g_nu + g_nu0
                  enddo
               enddo
               j = 1 + nat*nqn*(ich-1) + nat*(iqn-1) + (icent-1)
               srfc(j) = exp(-kappa * mu(3)) * g_nu
            enddo
            
         enddo
      enddo
         
      end
c     Notes:
c       1) intgl is conceptually a 4D array with size nsh * natom * nqn * nch where
c          nsh=(1+lmax_)^2. That is, the first dimension is for spherical harmoncis,
c          second dimension for atoms, third dimension for parabolic quantum numbers, and
c          the last dimension for ionized channels.
c       2) srfc, the structure factor, is conceptually a 3D array with size natom*nqn*nch.
c          That is, the first dimension is for atoms, second dimension for parabolic 
c          quantum numbers, and last dimension for ionized channels.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_passrot(bt0, gm0, protmat)

      implicit none
      double precision bt0, gm0            ! [input], orientation angles in degrees.
      double precision protmat(1:3, 1:3)

      double precision bt, gm              ! orientation angles in radians.
      double precision prot_y(3,3), prot_z(3,3), deg2rad
      parameter (deg2rad=0.01745329251994d0)
      
      
      bt = bt0 * deg2rad
      gm = gm0 * deg2rad
      
      prot_y(1:3, 1) = (/  cos(bt),    0.0d0,  sin(bt) /)    ! See eq. (1) in Ref. [1] and the discussion around it.
      prot_y(1:3, 2) = (/    0.0d0,    1.0d0,    0.0d0 /)
      prot_y(1:3, 3) = (/ -sin(bt),    0.0d0,  cos(bt) /)
      
      prot_z(1:3, 1) = (/  cos(gm), -sin(gm),    0.0d0 /)
      prot_z(1:3, 2) = (/  sin(gm),  cos(gm),    0.0d0 /)
      prot_z(1:3, 3) = (/    0.0d0,    0.0d0,    1.0d0 /)  
      
c      protmat = matmul(prot_y, prot_z)
      protmat = matmul(prot_z, prot_y)


cOLD      bt = bt0 * deg2rad
cOLD      gm = gm0 * deg2rad
cOLD      
cOLD      prot_y(1:3, 1) = (/  cos(bt),   0.0d0, -sin(bt) /)    ! See eq. (1) in Ref. [1] and the discussion around it.
cOLD      prot_y(1:3, 2) = (/    0.0d0,   1.0d0,    0.0d0 /)
cOLD      prot_y(1:3, 3) = (/  sin(bt),   0.0d0,  cos(bt) /)
cOLD      
cOLD      prot_z(1:3, 1) = (/  cos(gm), sin(gm),    0.0d0 /)
cOLD      prot_z(1:3, 2) = (/ -sin(gm), cos(gm),    0.0d0 /)
cOLD      prot_z(1:3, 3) = (/    0.0d0,   0.0d0,    1.0d0 /)  
cOLD      
cOLD      protmat = matmul(prot_z, prot_y)

      end
c     References:
c     [1] "Structure factors for tunneling ionization rates of diatomic molecules", http://dx.doi.org/10.1016/j.adt.2015.02.001
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_wignerf(j, m1, m2, beta, wgf)

      !!  d_{m1,m2}^j (beta)
      !!  beta is in radian
      implicit none
#include "errquit.fh"      
      integer j, m1, m2
      double precision beta
      double precision wgf
      
      integer k, kmax
      double precision nm1, nm2, nm3, nm4, dm1, dm2, dm3 ,dm4
      double precision cc, ss, gg, factorial
      external factorial
      
      
      if (abs(m1) > j .or. abs(m2) > j) then
         call errquit('WFAT: Arguments m1 and m2 of wfat_wignerf '//
     &        'must be such that abs(m1)>j and abs(m2)>j', 1,
     &        UERR)
      endif
      
      wgf = 0.0d0
      kmax = max(j-m1, j+m2, m2-m1)
      do k = 0, kmax
         if (j-m1 < k .or. j+m2 < k .or. k < m2-m1) then
            cycle
         endif
         
         nm1 = factorial(j + m2)
         nm2 = factorial(j - m2)
         nm3 = factorial(j + m1)
         nm4 = factorial(j - m1)
      
         dm1 = factorial(j + m2 - k)
         dm2 = factorial(k)
         dm3 = factorial(j - m1 - k)
         dm4 = factorial(k + m1 - m2)
      
         cc = (cos(beta/2))**(2*j - 2*k + m2 - m1)
         ss = (sin(beta/2))**(2*k - m2 + m1)
         
         gg = (-1)**(k-m2+m1) * sqrt(nm1 * nm2 * nm3 * nm4) /
     &        (dm1 * dm2 * dm3 * dm4) * cc * ss
      
         wgf = wgf + gg
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_gn(nat, nqn, mid, lmax, bt, gm, intgl, nopw,
     &     gn)

      implicit none
#include "errquit.fh"
      integer nat, nqn, mid(nqn), lmax
      double precision bt, gm                    ! [input], orientation angle in degree.
      double complex intgl((1+lmax)**2, nat, nqn)    !! 1) nsh * natom * nqn * nch
      logical nopw
      double complex gn(nat, nqn)                    !! 2)

      integer iqn, ic, l, m, i_lm, l0
      double precision btr, gmr, deg2rad, wgf
      double complex gn0, czero, imnum
      parameter ( deg2rad = 0.01745329251994d0 )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )


      if (nopw .and. lmax /= 0)
     &     call errquit('WFAT: If nopw is .true., then lmax must be '//
     &     'equal to zero (wfat_calc_gn).', 43, UERR)
      
      
      btr = bt * deg2rad
      gmr = gm * deg2rad
      
c     ==== Begin calculating asymptotic coefficients ====
      do iqn = 1, nqn

         do ic = 1, nat
            gn(ic, iqn) = czero
            if (nopw) then
               l0 = 0
            else
               l0 = abs(mid(iqn))
            endif
            
            do l = l0, lmax
               do m = -l, l
                  if (nopw) then
                     wgf = 1.0d0
                  else
                     call wfat_wignerf(l, mid(iqn), m, btr, wgf)
                  endif
                  
                  i_lm = l**2 + l + m + 1
                  gn0 = intgl(i_lm,ic,iqn) * wgf *
     &                  exp(-imnum*m*gmr)
                  gn(ic, iqn) = gn(ic, iqn) + gn0
               enddo
            enddo

         enddo
      enddo
         
      end
c=======================================================
