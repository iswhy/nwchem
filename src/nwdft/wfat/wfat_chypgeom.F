c=======================================================
c=======================================================
      subroutine wfat_chgm_1pt(r, lmin_, lmax_, a_, b_, lmb, chgm)

      implicit none
#include "stdio.fh"      
      integer lmin_, lmax_
      double precision r, a_, b_, lmb
      double precision chgm(lmin_:lmax_)
      
      integer l
      double precision a, aa, b, bb, b0, z, chgm_1, chgm_2, chgm_3
      double complex zz, uu
      
      
c     ==== Prepare initial values for recurrence calculation ====
      z = lmb * r
      zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
      a = a_ + (lmax_ + 1)
      b = b_ + 2*(lmax_ + 1)
      
      b0 = b - 0
c      call wfat_chgm(a, b0, z, chgm_3)
      call wfat_cchg(a, b0, zz, uu)           ! 1)
      chgm_3 = real(uu)
      
      b0 = b - 1
c      call wfat_chgm(a, b0, z, chgm_2)
      call wfat_cchg(a, b0, zz, uu)           ! 1)
      chgm_2 = real(uu)
      
      b0 = b - 2
c      call wfat_chgm(a, b0, z, chgm_1)
      call wfat_cchg(a, b0, zz, uu)           ! 1)
      chgm_1 = real(uu)
      
      
c     ==== Calculate the confluent hypergeometric function by recurrence ====
cc      if (abs(chgm_3) < 5.0d2) then
cc         call wfat_chgm_1pt_recur(lmin_, lmax_, a_, b_, z,
cc     &                            chgm_1, chgm_2, chgm_3, chgm)
cc      else
         do l = lmin_, lmax_
            aa = l + a_
            bb = 2*l + b_
c            call wfat_chgm(aa, bb, z, chgm(l))
            call wfat_cchg(aa, bb, zz, uu)           ! 1)
            chgm(l) = real(uu)
         enddo
cc      endif
         
      
c     ==== Scan for too big values ====
      do l = lmin_, lmax_
c         if (abs(chgm(l)) > 1.0d+20 .or. isnan(chgm(l)) ) then
         if ( isnan(chgm(l)) ) then        
c              print *, "rmax = ", r0, l, chgm(l)
            chgm(l) = 0.0d0
         endif
      enddo
      
      end
c     NOTES:
c     1) wfat cchg is used instead of wfat_chgm that is more specifically
c        used for real arguments because for some reason the latter
c        produced NaN in calculating the structure factors.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chgm_1pt_recur(m, ll, a0, b0, z, chgm_1_, chgm_2_,
     &                               chgm_3_, chgm)
      
c      !! DESCRIPTION:
c      !!   chgm_1_ = M(ll+1+a0, 2ll+0+b0, z)
c      !!   chgm_2_ = M(ll+1+a0, 2ll+1+b0, z)
c      !!   chgm_3_ = M(ll+1+a0, 2ll+2+b0, z)
c      !!   chgm    = M(l+a0,    2l+b0,    z)  with l = m, m+1, ..., ll
c      !!   m >= 0
      
      implicit none
#include "errquit.fh"
      integer m, ll
      double precision a0, b0, z
      double precision chgm_1_, chgm_2_, chgm_3_
      double precision chgm(m:ll)
      
      integer l
      double precision a, b, chgm_1, chgm_2, chgm_3, chgm_4


      if (m < 0) then
         call errquit('WFAT: Argument m of wfat_chgm_1pt_recur '//
     &                'cannot be negative.', 1, UERR)
      endif
      
            
      do l = ll, m, -1
         a = a0 + l
         b = b0 + 2*l + 1       
      
         !!==== chgm_4 ====!! Use eq. (13.4.4) of Ref. [1]
         chgm_4 = (b * chgm_2_ - z * chgm_3_) / b
      
         !!==== chgm_3 ====!! Use eq. (13.4.4) of Ref. [1]
         b = b - 1
         chgm_3 = (b * chgm_1_ - z * chgm_2_) / b
         chgm(l) = chgm_3
      
         !!==== Avoid unnecessary calculations ====!!
         if (l == m) exit
      
         !!==== chgm_2 ====!! Use eq. (13.4.2) of Ref. [1]
         b = b - 1
         chgm_2 = -(-(b+1)*(b+z)*chgm_3 + z*(b+1-a)*chgm_4) / b / (b+1)
      
         !!==== chgm_1 ====!! Use eq. (13.4.2) of Ref. [1]
         b = b - 1
         chgm_1 = -(-(b+1)*(b+z)*chgm_2 + z*(b+1-a)*chgm_3) / b / (b+1)
      
         !!==== For next iteration ====!!
         chgm_3_ = chgm_3
         chgm_2_ = chgm_2
         chgm_1_ = chgm_1
      enddo
            
      end
c     References:
c        [1] Handbook of Mathematical Functions by Abramowitz and Stegun.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chgm ( a, b, x, hg )

c*********************************************************************72
c
cc CHGM computes the confluent hypergeometric function M(a,b,x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    27 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HG, the value of M(a,b,x).
c
      implicit none

      double precision a
      double precision a0
      double precision a1
      double precision b
      double precision hg
      double precision hg1
      double precision hg2
      integer i
      integer j
      integer k
      integer la
      integer m
      integer n
      integer nl
      double precision pi
      double precision r
      double precision r1
      double precision r2
      double precision rg
      double precision sum1
      double precision sum2
      double precision ta
      double precision tb
      double precision tba
      double precision x
      double precision x0
      double precision xg
      double precision y0
      double precision y1

      pi = 3.141592653589793D+00
      a0 = a
      a1 = a
      x0 = x
      hg = 0.0D+00

      if ( b .eq. 0.0D+00 .or. b .eq. - abs ( int ( b ) ) ) then
        hg = 1.0D+300
      else if ( a .eq. 0.0D+00 .or. x .eq. 0.0D+00 ) then
        hg = 1.0D+00
      else if ( a .eq. -1.0D+00 ) then
        hg = 1.0D+00 - x / b
      else if ( a .eq. b ) then
        hg = exp ( x )
      else if ( a - b .eq. 1.0D+00 ) then
        hg = ( 1.0D+00 + x / b ) * exp ( x )
      else if ( a .eq. 1.0D+00 .and. b .eq. 2.0D+00 ) then
        hg = ( exp ( x ) - 1.0D+00 ) / x
      else if ( a .eq. int ( a ) .and. a .lt. 0.0D+00 ) then
        m = int ( - a )
        r = 1.0D+00
        hg = 1.0D+00
        do k = 1, m
          r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
          hg = hg + r
        end do
      end if

      if ( hg .ne. 0.0D+00 ) then
        return
      end if

      if ( x .lt. 0.0D+00 ) then
        a = b - a
        a0 = a
        x = abs ( x )
      end if

      if ( a .lt. 2.0D+00 ) then
        nl = 0
      end if

      if ( 2.0D+00 .le. a ) then
        nl = 1
        la = int ( a )
        a = a - la - 1.0D+00
      end if

      do n = 0, nl

        if ( 2.0D+00 .le. a0 ) then
          a = a + 1.0D+00
        end if

        if ( x .le. 30.0D+00 + abs ( b ) .or. a .lt. 0.0D+00 ) then

          hg = 1.0D+00
          rg = 1.0D+00
          do j = 1, 500
            rg = rg * ( a + j - 1.0D+00 ) 
     &        / ( j * ( b + j - 1.0D+00 ) ) * x
            hg = hg + rg
            if ( abs ( rg / hg ) .lt. 1.0D-15 ) then
              go to 10
            end if
          end do

10        continue

        else

          call wfat_gammaf ( a, ta )
          call wfat_gammaf ( b, tb )
          xg = b - a
          call wfat_gammaf ( xg, tba )
          sum1 = 1.0D+00
          sum2 = 1.0D+00
          r1 = 1.0D+00
          r2 = 1.0D+00
          do i = 1, 8
            r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
            r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
            sum1 = sum1 + r1
            sum2 = sum2 + r2
          end do
          hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
          hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
          hg = hg1 + hg2

        end if

        if ( n .eq. 0 ) then
          y0 = hg
        else if ( n .eq. 1 ) then
          y1 = hg
        end if

      end do

      if ( 2.0D+00 .le. a0 ) then
        do i = 1, la - 1
          hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
          y0 = y1
          y1 = hg
          a = a + 1.0D+00
        end do
      end if

      if ( x0 .lt. 0.0D+00 ) then
        hg = hg * exp ( x0 )
      end if

      a = a1
      x = x0

      return
      end      
c=======================================================

      
c=======================================================
c=======================================================  
      subroutine wfat_cchg ( a, b, z, chg )
      
c      !*****************************************************************************80
c      !
c      !! CCHG computes the confluent hypergeometric function.
c      !
c      ! Website : https://people.sc.fsu.edu/~jburkardt/f77_src/special_functions/special_functions.html
c      !
c      !  Discussion:
c      !
c      !    This function computes the confluent hypergeometric function
c      !    M(a,b,z) with real parameters a, b and complex argument z.
c      !
c      !  Licensing:
c      !
c      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
c      !    they give permission to incorporate this routine into a user program
c      !    provided that the copyright is acknowledged.
c      !
c      !  Modified:
c      !
c      !    26 July 2012
c      !
c      !  Author:
c      !
c      !    Shanjie Zhang, Jianming Jin
c      !   
c      !  Edited by:
c      !
c      !    Imam Wahyutama
c      !
c      !  Reference:
c      !
c      !    Shanjie Zhang, Jianming Jin,
c      !    Computation of Special Functions,
c      !    Wiley, 1996,
c      !    ISBN: 0-471-11963-6,
c      !    LC: QA351.C45.
c      !
c      !  Parameters:
c      !
c      !    Input, real ( kind = 8 ) A, B, parameter values.
c      !
c      !    Input, complex ( kind = 8 ) Z, the argument.
c      !
c      !    Output, complex ( kind = 8 ) CHG, the value of M(a,b,z).
c      !
      
      implicit none
      
      double precision a        ! [input]
      double precision a0
      double precision a1
      double precision b        ! [input]
      double precision ba
      double complex cfac
      double complex chg        ! [output]
      double complex chg1
      double complex chg2
      double complex chw
      double complex ci
      double complex cr
      double complex cr1
      double complex cr2
      double complex crg
      double complex cs1
      double complex cs2
      double complex cy0
      double complex cy1
      double precision g1
      double precision g2
      double precision g3
      integer i
      integer j
      integer k
      integer la
      integer m
      integer n
      integer nl
      integer ns
      double precision phi
      double precision pi
      double precision x
      double precision x0
      double precision y
      double complex z         ! [in/output]
      double complex z0
      
      
      pi = 3.141592653589793D+00
      ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
      a0 = a
      a1 = a
      z0 = z
      
      
      if ( b == 0.0D+00 .or. b == - int ( abs ( b ) ) ) then
         chg = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
      else if ( a == 0.0D+00 .or. z == 0.0D+00 ) then
         chg = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
      else if ( a == -1.0D+00 ) then
         chg = 1.0D+00 - z / b
      else if ( a == b ) then
         chg = exp ( z )
      else if ( a - b == 1.0D+00 ) then
         chg = ( 1.0D+00 + z / b ) * exp ( z )
      else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
         chg = ( exp ( z ) - 1.0D+00 ) / z
      else if ( a == int ( a ) .and. a < 0.0D+00 ) then
         m = int ( - a )
         cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
         chg = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
         do k = 1, m
            cr = cr * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * z
            chg = chg + cr
         end do
      else
      
         x0 = real ( z, kind = 8 )
         if ( x0 < 0.0D+00 ) then
            a = b - a
            a0 = a
            z = - z
         end if
      
         if ( a < 2.0D+00 ) then
            nl = 0
         else
            nl = 1
            la = int ( a )
            a = a - la - 1.0D+00
         end if
      
         do n = 0, nl
      
            if ( 2.0D+00 <= a0 ) then
               a = a + 1.0D+00
            end if
      
            if ( cdabs ( z ) < 20.0D+00 + abs ( b ) .or. a < 0.0D+00 )
     &           then
      
               chg = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
               crg = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
               do j = 1, 500
                  crg = crg * ( a + j - 1.0D+00 ) /
     &                  ( j * ( b + j - 1.0D+00 ) ) * z
                  chg = chg + crg
                  if ( abs ( ( chg - chw ) / chg ) < 1.0D-15 ) then
                     exit
                  end if
                  chw = chg
               end do
      
            else
      
               call wfat_gammaf ( a, g1 )
               call wfat_gammaf ( b, g2 )
               ba = b - a
               call wfat_gammaf ( ba, g3 )
               cs1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
               cs2 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
               cr1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
               cr2 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
      
               do i = 1, 8
                  cr1 = - cr1 * (     a + i - 1.0D+00 ) *
     &                  ( a - b + i ) / ( z * i )
                  cr2 =   cr2 * ( b - a + i - 1.0D+00 ) *
     &                  ( i - a ) / ( z * i )
                  cs1 = cs1 + cr1
                  cs2 = cs2 + cr2
               end do
      
               x = real ( z, kind = 8 )
               y = imag ( z )
      
               if ( x == 0.0D+00 .and. 0.0D+00 <= y ) then
                  phi = 0.5D+00 * pi
               else if ( x == 0.0D+00 .and. y <= 0.0D+00 ) then
                  phi = -0.5D+00 * pi
               else
                  phi = atan ( y / x )
               end if
      
               if ( -1.5D+00 * pi < phi .and. phi <= -0.5 * pi ) then
                  ns = -1
               else if ( -0.5D+00 * pi < phi .and. phi < 1.5D+00 * pi )
     &                 then
                  ns = 1
               end if
      
               if ( y == 0.0D+00 ) then
                  cfac = cos ( pi * a )
               else
                  cfac = exp ( ns * ci * pi * a )
               end if
      
               chg1 = g2 / g3 * z ** ( - a ) * cfac * cs1
               chg2 = g2 / g1 * exp ( z ) * z ** ( a - b ) * cs2
               chg = chg1 + chg2
      
            end if
      
            if ( n == 0 ) then
               cy0 = chg
            else if ( n == 1 ) then
               cy1 = chg
            end if
      
         end do
      
         if ( 2.0D+00 <= a0 ) then
            do i = 1, la - 1
               chg = ( ( 2.0D+00 * a - b + z ) * cy1 + ( b - a ) * cy0 )
     &               / a
               cy0 = cy1
               cy1 = chg
               a = a + 1.0D+00
            end do
         end if
      
         if ( x0 < 0.0D+00 ) then
            chg = chg * exp ( - z )
         end if
      
      end if
      
      a = a1
      z = z0
      
      return
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_gammaf ( x, ga )
      
c      !*****************************************************************************80
c      !
c      !! GAMMA evaluates the Gamma function.
c      !
c      ! Website : https://people.sc.fsu.edu/~jburkardt/f77_src/special_functions/special_functions.html
c      !
c      !  Licensing:
c      !
c      !    The original FORTRAN77 version of this routine is copyrighted by 
c      !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
c      !    incorporate this routine into a user program that the copyright 
c      !    is acknowledged.
c      !
c      !  Modified:
c      !
c      !    08 September 2007
c      !
c      !  Author:
c      !
c      !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
c      !    FORTRAN90 version by John Burkardt.
c      !
c      !  Reference:
c      !
c      !    Shanjie Zhang, Jianming Jin,
c      !    Computation of Special Functions,
c      !    Wiley, 1996,
c      !    ISBN: 0-471-11963-6,
c      !    LC: QA351.C45
c      !
c      !  Parameters:
c      !
c      !    Input, real ( kind = 8 ) X, the argument.
c      !    X must not be 0, or any negative integer.
c      !
c      !    Output, real ( kind = 8 ) GA, the value of the Gamma function.
c      !
      implicit none
      double precision g(26) 
      double precision ga
      double precision gr
      integer k
      integer m
      integer m1
      double precision pi
      parameter (pi = 3.141592653589793D+00)
      double precision r
      double precision x
      double precision z      
      


      g = (/ 1.0D+00, 
     &       0.5772156649015329D+00, 
     &      -0.6558780715202538D+00, 
     &      -0.420026350340952D-01, 
     &       0.1665386113822915D+00, 
     &      -0.421977345555443D-01, 
     &      -0.96219715278770D-02, 
     &       0.72189432466630D-02, 
     &      -0.11651675918591D-02, 
     &      -0.2152416741149D-03, 
     &       0.1280502823882D-03,  
     &      -0.201348547807D-04, 
     &      -0.12504934821D-05, 
     &       0.11330272320D-05, 
     &      -0.2056338417D-06,  
     &       0.61160950D-08, 
     &       0.50020075D-08, 
     &      -0.11812746D-08, 
     &       0.1043427D-09,  
     &       0.77823D-11, 
     &      -0.36968D-11, 
     &       0.51D-12, 
     &      -0.206D-13, 
     &      -0.54D-14, 
     &       0.14D-14, 
     &       0.1D-15 /)
      
      if ( x == aint ( x ) ) then
      
        if ( 0.0D+00 < x ) then
          ga = 1.0D+00
          m1 = int ( x ) - 1
          do k = 2, m1
            ga = ga * k
          enddo
        else
          ga = 1.0D+300
        endif
      
      else
      
        if ( 1.0D+00 < abs ( x ) ) then
          z = abs ( x )
          m = int ( z )
          r = 1.0D+00
          do k = 1, m
            r = r * ( z - real ( k, kind = 8 ) )
          enddo
          z = z - real ( m, kind = 8 )
        else
          z = x
        endif
      
        gr = g(26)
        do k = 25, 1, -1
          gr = gr * z + g(k)
        enddo
      
        ga = 1.0D+00 / ( gr * z )
      
        if ( 1.0D+00 < abs ( x ) ) then
          ga = ga * r
          if ( x < 0.0D+00 ) then
            ga = - pi / ( x* ga * sin ( pi * x ) )
          endif
        endif
      
      endif
      
      return
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_psif ( x, ps )
c*********************************************************************72
c
cc PSI computes the Psi function.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision PS, the value of the function.
c
      implicit none

      double precision a1
      double precision a2
      double precision a3
      double precision a4
      double precision a5
      double precision a6
      double precision a7
      double precision a8
      double precision el
      integer k
      integer n
      double precision pi
      double precision ps
      double precision s
      double precision x
      double precision x2
      double precision xa

      xa = abs ( x )
      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      s = 0.0D+00

      if ( x .eq. int ( x ) .and. x .le. 0.0D+00 ) then

        ps = 1.0D+300
        return

      else if ( xa .eq. int ( xa ) ) then

        n = int ( xa )
        do k = 1, n - 1
          s = s + 1.0D+00 / k
        end do
        ps = - el + s

      else if ( xa + 0.5D+00 .eq. int ( xa + 0.5D+00 ) ) then

        n = int ( xa - 0.5D+00 )
        do k = 1, n
          s = s + 1.0D+00 / ( 2.0D+00 * k - 1.0D+00 )
        end do
        ps = - el + 2.0D+00 * s - 1.386294361119891D+00

      else

        if ( xa .lt. 10.0D+00 ) then
          n = 10 - int ( xa )
          do k = 0, n - 1 
            s = s + 1.0D+00 / ( xa + k )
          end do
          xa = xa + n
        end if

        x2 = 1.0D+00 / ( xa * xa )
        a1 = -0.8333333333333D-01
        a2 = 0.83333333333333333D-02
        a3 = -0.39682539682539683D-02
        a4 = 0.41666666666666667D-02
        a5 = -0.75757575757575758D-02
        a6 = 0.21092796092796093D-01
        a7 = -0.83333333333333333D-01
        a8 = 0.4432598039215686D+00
        ps = log ( xa ) - 0.5D+00 / xa + x2 * (((((((
     &      a8   * x2
     &    + a7 ) * x2
     &    + a6 ) * x2
     &    + a5 ) * x2
     &    + a4 ) * x2
     &    + a3 ) * x2
     &    + a2 ) * x2
     &    + a1 )
        ps = ps-s

      end if

      if ( x .lt. 0.0D+00 ) then
        ps = ps - pi * cos ( pi * x ) / sin ( pi * x ) - 1.0D+00 / x
      end if

      return
      end      
c=======================================================
