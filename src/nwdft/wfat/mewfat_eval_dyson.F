c=======================================================
c=======================================================
      subroutine mewfat_eval_dyson(sx, n0, pcf, rcf, phi, dy)

      implicit none
      integer sx                    ! [input], the spin of electron removed from the neutral.
      integer n0(2)                 ! [input], no. of electrons in neutral.
      double precision pcf(n0(sx))  ! [input], P det. coefficients.
      double precision rcf          ! [input], R det. coefficients.
      double precision phi(n0(sx))  ! [input], occupied MO's of the neutral in the sx channel evaluated at the evaluation point.
      double precision dy           ! [output], the Dyson orbital at the evaluation point.

      integer j, dl, ntot0
      double precision a

      
      dy = 0.0d0
      do j = 1, n0(sx)
         a = (-1)**j * pcf(j) * phi(j)
         dy = dy + a
      enddo

      dl = 0
      if (sx==2) dl = 1
      ntot0 = sum(n0)
      dy = (-1)**(ntot0+dl*n0(1)) / sqrt(dble(ntot0)) * rcf * dy
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_norm_dyson(sx, n0, pcf, rcf, norm)

      implicit none
      integer sx                    ! [input], the spin of electron removed from the neutral.
      integer n0(2)                 ! [input], no. of electrons in neutral.
      double precision pcf(n0(sx))  ! [input], P det. coefficients.
      double precision rcf          ! [input], R det. coefficients.
      double precision norm         ! [output], the Dyson orbital norm.

      integer j, ntot0
      double precision norm0

      
      norm = 0.0d0
      do j = 1, n0(sx)
         norm0 = pcf(j)**2
         norm = norm + norm0
      enddo

      ntot0 = sum(n0)
      norm = norm / dble(ntot0) * rcf**2
      
      end
c=======================================================
      
      
