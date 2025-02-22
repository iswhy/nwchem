c=======================================================
c=======================================================
      subroutine wfat_calc_o(n, m, erg, zq, o)

      implicit none
#include "stdio.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      double precision o(2)

      double precision kp, gmm, bt0, oi(2)

      
      kp = sqrt(2 * abs(erg))
      gmm = (1.0d0 - m**2) / 4.0d0
      bt0 = zq - kp*(n + (abs(m)+1)/2.0d0)
      
      call wfat_calc_o_intl(n, m, erg, zq, oi)
cdebug      call wfat_calc_o_intl_trapz(n, m, erg, zq, oi)
      
      o(1) = 1/(4*kp**3)*(1 + gmm**2/2 - gmm*(4*bt0/kp-(bt0/kp)**2)
     &       + 4*bt0/kp + (7.0d0/2)*(bt0/kp)**2 - 4*(bt0/kp)**3 
     &       + (bt0/kp)**4/2) + oi(1)
      o(2) = -1.0d0/(2*kp**2) + oi(2)

cdebug      write(luout, '("oo = ", 2es21.12)') oi(1), oi(2)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl(n, m, erg, zq, oi)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      double precision oi(2)

      logical l1, l2
      integer k, NNGR, l_wt, k_wt, l_pt, k_pt
      double precision eb1, eb2, kp, oi_p1(2), oi_p2(2), oi_p3(2), THR
      parameter (NNGR=3, THR=1.0d-4)
      integer ngr(NNGR)
      double precision oi0(2,NNGR)


      ngr = (/ 40, 50, 60 /)
      kp = sqrt(2 * abs(erg))
      eb1 = 1.0d0 / kp              ! so that zb = 1.0d0 inside wfat_calc_o_intl_part* below.
      eb1 = 0.3d0 / kp
      eb2 = 30.0d0 / kp

      
      call wfat_calc_o_intl_part1(n, m, erg, zq, eb1, oi_p1)
c      write(luout, '("oi_p1 : ", 2es21.12)') oi_p1          ! For debugging

      
      if (.not. MA_alloc_get(MT_DBL, maxval(ngr), 'lgl_pt', l_pt, k_pt))
     &     call errquit('WFAT: lgl_pt allocation failed', 32, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, maxval(ngr), 'lgl_wt', l_wt, k_wt))
     &     call errquit('WFAT: lgl_wt allocation failed', 33, MA_ERR)
      do k = 1, NNGR
         call wfat_lgl(-1.0d0, 1.0d0, ngr(k)-1, dbl_mb(k_pt),
     &        dbl_mb(k_wt))
         
         call wfat_calc_o_intl_part2(n, m, erg, zq, eb1, eb2, ngr(k),  
     &        dbl_mb(k_pt), dbl_mb(k_wt), oi0(:,k))
c         write(luout, '("oi_p2 : ", i3, 2es21.12)') k, oi0(:,k)          ! For debugging
         if (k > 1) then
            l1 = (abs(oi0(1,k) - oi0(1,k-1)) < THR)
            l2 = (abs(oi0(2,k) - oi0(2,k-1)) < THR)
            if (l1 .and. l2) then
               oi_p2 = oi0(:,k)
               exit
            endif
         endif

         if (k == NNGR)
     &        call errquit('WFAT: Calculations of oi_p2 never '//
     &        'converged', 19, UERR)
      enddo


c      write(luout, '("hered1")')
      call wfat_calc_o_intl_part3(n, m, erg, zq, eb2, oi_p3)
c      write(luout, '("oi_p3 : ", 2es21.12)') oi_p3          ! For debugging

      oi = oi_p1 + oi_p2 + oi_p3


      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate lgl_pt', 68, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate lgl_wt', 69, MA_ERR)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl_part1(n, m, erg, zq, eb, oi_p1)
c     DESCRIPTION:
c     Computes the small-$\eta$ part of $o_\nu$. The small-$\eta$ part of $o_\nu$ can
c     be written so as to separate the integral over the integrand that contains
c     $R_\nu^{(0)} O_\nu^{(0)}$ (hence the integral cannot be done analytically) and
c     over one that contains analytically integrable terms. The later part is referred
c     to as the direct integration part below because it can be directly integrated
c     analytically. To integrate the former part, we make use of the fact that the
c     integral only runs over small values of $\eta$ (in fact, we restricted it to
c     $\eta < 1/\kappa$), and thus the asymptotic forms of $R_\nu^{(0)}$ and
c     $O_\nu^{(0)}$ for small argument can be employed to obtain the non-analytically
c     integrable part.
      implicit none
#include "stdio.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      double precision eb
      double precision oi_p1(2)

      integer i, j, pmax
      double precision kp, zb, gmm, bt0, bt1_1, bt1_2
      double precision part_pow(2), part_dir(2)
      parameter (pmax=40)
      double precision o0i(pmax,2), o0il(pmax,2)


      kp = sqrt(2 * abs(erg))
      zb = kp * eb
      bt0 = zq - kp*(n + (abs(m)+1)/2.0d0)
      bt1_1 = -1.0d0/4/kp**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))
      gmm = (1.0d0 - m**2) / 4.0d0
      

c     ==== Do the power expansion part ====
      call wfat_o0_i_coeff_nomu(n, m, erg, zq, pmax, o0i, o0il)
      do j = 1, 2
         part_pow(j) = 0.0d0
         do i = 1, pmax
            part_pow(j) = part_pow(j) +
     &                    (o0i(i,j) + log(zb)*o0il(i,j)) * zb**i
         enddo
      enddo
cdebug      write(luout, '("part_pow = ", 2es18.8)') part_pow

      
c     ==== Do the direct integration part ====
      part_dir(1) = - 1/kp * ( eb**2/8 + bt0*eb/(2*kp**2) +
     &              (bt1_1 + gmm/(2*kp**2) + (3*bt0**2)/(2*kp**4)) *
     &              log(eb+1) )
      part_dir(2) = 1/kp * ( eb/2 - (bt1_2 - bt0/kp**2) * log(eb+1) )
cdebug      write(luout, '("part_dir = ", 2es18.8)') part_dir

      
      oi_p1 = part_pow + part_dir

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl_part2(n, m, erg, zq, eb1, eb2, ngr,
     &           pt, wt, oi)

      implicit none
      integer n
      integer m
      double precision erg
      double precision zq
      double precision eb1
      double precision eb2
      integer ngr
      double precision pt(ngr)
      double precision wt(ngr)
      double precision oi(2)

      integer i
      double precision et, wt0, ee(2)


      oi = 0.0d0
      do i = 1, ngr
         et = 0.5d0*((eb2 - eb1)*pt(i) + eb2 + eb1)
         wt0 = (eb2 - eb1)/2 * wt(i)

         call wfat_calc_o_intd(et, n, m, erg, zq, ee)
         oi = oi + wt0 * ee
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl_part2_OLD(n, m, erg, zq, eb, ngr,
     &           oi)

      implicit none
      integer n
      integer m
      double precision erg
      double precision zq
      double precision eb
      integer ngr
      double precision oi(2)

      integer i
      double precision kp, a, aa, eta, ee(2), pt(ngr), wt(ngr)

      
c     ==== Update the library of Gauss-Laguerre grid and ====
c     ====   nodes used in the call of wfat_chgu later   ====
      kp = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kp
      if (a-1 < 0.0d0) then
         aa = (a + abs(int(a-1)) + 2) - 1
         call wfat_update_lfqlib(aa)
         call wfat_update_lfqlib(aa-1)
      else
         aa = a - 1
         call wfat_update_lfqlib(aa)
      endif

      
c     ==== Calculate the integral using Gauss-Laguerre ====
c     ====  quadrature (whose purpose is unrelated to  ====
c     ====           the one mentioned above)          ====
      call wfat_lf_quadrature_rule(ngr, 0.0d0, pt, wt)
      oi = 0.0d0
      do i = 1, ngr
         eta = pt(i)/kp + eb
         call wfat_calc_o_intd_exp(eta, n, m, erg, zq, .true., ee)     ! 1)

         oi = oi + wt(i) * ee
      enddo
      oi = oi * exp(-kp*eb)/kp**2

      end
c     NOTES:
c     1) Inside wfat_calc_o_intd is where the grids and nodes updated in
c        the two sections above matters.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intd_exp(eta, n, m, erg, zq, exp1, ee)
c     With the exponential factor taken out.
      implicit none
#include "errquit.fh"
      double precision eta
      integer n
      integer m
      double precision erg
      double precision zq
      logical exp1
      double precision ee(2)

      integer j, pmax, md
      double precision kp, a, b, g, d, bt0, bt1_1, bt1_2, gmm, chgu, gm
      double precision v1, v2, ps1, ps2, ps3, dj, u0, chg_mu, expf
      double precision factorial, rise_fact
      double complex zz, chgm
      parameter (pmax=40)
      external factorial, rise_fact


      kp = sqrt(2 * abs(erg))
      if (kp*eta < 1.0d-14)
     &     call errquit('WFAT: When m=0, the argument eta of '//
     &     'wfat_calc_o_intd_exp must be such that kp*eta >= 1.0E-14',
     &     28, UERR)
      if (eta <= 0.0d0)
     &     call errquit('WFAT: The argument eta of '//
     &     'wfat_calc_o_intd_exp must be positive', 36, UERR)


      a = 1 + abs(m) + n - zq/kp
      b = dble(abs(m) + 1)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))
      bt0 = zq - kp*(n + (abs(m)+1)/2.0d0)
      bt1_1 = -1.0d0/4/kp**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))
      gmm = (1.0d0 - m**2) / 4.0d0


c     u = (kp*eta)^|m| * U(a,b,kp*eta)
      if (kp*eta > 0.0d0 .and. kp*eta < 1.0d0) then
         chgu = 0.0d0
         call wfat_gammaf(a-abs(m), gm)
         v1 = (-1)**b/factorial(abs(m))/gm
         do j = 0, pmax
            call wfat_psif(a+j, ps1)
            call wfat_psif(dble(1+j), ps2)
            call wfat_psif(dble(1+j+abs(m)), ps3)
            dj = ps1 - ps2 - ps3
            u0 = v1 * rise_fact(a,j) / rise_fact(b,j) / factorial(j)
     &           * (kp*eta)**(abs(m)+j) * (log(kp*eta) + dj)
            chgu = chgu + u0
         enddo

         if (abs(m) > 0) then
            v2 = 1.0d0 / g
            do j = 1, abs(m)
               u0 = v2 * factorial(j-1) * rise_fact(1-a+j,abs(m)-j)
     &              / factorial(abs(m)-j) * (kp*eta)**(abs(m)-j)
               chgu = chgu + u0
            enddo
         endif
      else
         call wfat_chgu(a, b, kp*eta, chgu, md)
         chgu = chgu * (kp*eta)**(abs(m))
      endif

      zz = cmplx(kp*eta, 0.0d0, kind=kind(1.0d0))
      call wfat_cchg(a, b, zz, chgm)
      chg_mu = real(chgm) * chgu

      
      if (exp1) then
         expf = 1.0d0
      else
         expf = exp(-kp*eta)
      endif
      
      ee(1) = d * (bt1_1*kp + kp*eta**2/4) * expf * chg_mu
     &        - (expf*exp(kp*eta)) * (eta/4 + bt0/2/kp**2 + (bt1_1 
     &        + gmm/2/kp**2 + (3*bt0**2)/(2*kp**4)) / (eta+1))

      ee(2) = d * (bt1_2*kp - kp*eta/2) * expf * chg_mu
     &        + (expf*exp(kp*eta)) * (1.0d0/2 - (bt1_2 - bt0/kp**2)
     &        / (eta+1))

      ee = ee / kp

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intd(eta, n, m, erg, zq, ee)

      implicit none
#include "errquit.fh"
      double precision eta
      integer n
      integer m
      double precision erg
      double precision zq
      double precision ee(2)

      double precision kp, bt0, bt1_1, bt1_2, gmm, r0o0
      

      kp = sqrt(2 * abs(erg))
      if (kp*eta < 1.0d-14)
     &     call errquit('WFAT: When m=0, the argument eta of '//
     &     'wfat_calc_o_intd must be such that kp*eta >= 1.0E-14',
     &     28, UERR)
      if (eta <= 0.0d0)
     &     call errquit('WFAT: The argument eta of wfat_calc_o_intd '//
     &     'must be positive', 36, UERR)


      bt0 = zq - kp*(n + (abs(m)+1)/2.0d0)
      bt1_1 = -1.0d0/4/kp**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))
      gmm = (1.0d0 - m**2) / 4.0d0
      

      call wfat_calc_r0o0(eta, n, m, erg, zq, r0o0)
      
      ee(1) = (bt1_1/eta + eta/4) * r0o0
     &        - eta/4 - bt0/2/kp**2 - (bt1_1 + gmm/2/kp**2 
     &        + (3*bt0**2)/(2*kp**4)) / (eta+1)

      ee(2) = (bt1_2/eta - 1.0d0/2) * r0o0
     &        + 1.0d0/2 - (bt1_2 - bt0/kp**2) / (eta+1)

      ee = ee / kp

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl_part3(n, m, erg, zq, eb, oi_p3)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      double precision eb
      double precision oi_p3(2)

      integer j, i, pmax
      double precision kp, zb
      parameter(pmax=30)
      double precision cf(pmax,2)


      kp = sqrt(2 * abs(erg))
      zb = kp * eb

      if (zb < 30.0d0 - 1.0d-13) then
         write(luout, '("zb = ", f21.14)') zb
         call errquit('WFAT: The argument eb of '//
     &     'wfat_calc_o_intl_part3 must be such that sqrt(2*|erg|)*eb'//
     &        ' is bigger than or equal to 30.0', 41, UERR)
      endif
      
      

c      write(luout, '("herew1")')
      call wfat_o_largee_coeff(n, m, erg, zq, pmax, cf)
c      write(luout, '("herew2")')

      do j = 1, 2
         oi_p3(j) = 0.0d0
         do i = 1, pmax
            oi_p3(j) = oi_p3(j) + cf(i,j) * (1/zb)**i
         enddo
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_o_largee_coeff(n, m, erg, zq, pmax, cf)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      integer pmax                   ! [input]
      double precision cf(pmax,2)    ! [output]

      logical l1, l2, l3
      integer i, j, k, k1, k2, k3
      double precision kp, a, b, bt0, bt1_1, bt1_2, gmm, w, w1, w2
      double precision rise_fact, factorial, binom_cf
      external rise_fact, factorial, binom_cf


      if (pmax < 1)
     &     call errquit('WFAT: Argument pmax of wfat_o_largee_coeff '//
     &     'cannot be smaller than 1.', 48, UERR)
      if (n < 0) call errquit('Argument n of '//
     &     'wfat_o_largee_coeff must be zero or positive.', 28, UERR)
      
      
      kp = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kp
      b = dble(abs(m) + 1)
      bt0 = zq - kp*(n + (abs(m)+1)/2.0d0)
      bt1_1 = -1.0d0/4/kp**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))
      gmm = (1.0d0 - m**2) / 4.0d0      


c      write(luout, '("herew3")')
      cf = 0.0d0
      do i = 0, pmax
         do j = 0, pmax
            k1 = i + j
            k2 = i + j - 2
            k3 = i + j - 1

c            write(luout, '(2i5)') i, j

            l1 = (-i-j-1 <= -2 .and. k1 <= pmax)
            l2 = (1-i-j <= -2 .and. k2 <= pmax)
            l3 = (-i-j <= -2 .and. k3 <= pmax)
            if (l1 .or. l2 .or. l3)
     &           w = rise_fact(b-a,i)*rise_fact(1-a,i)/factorial(i) *
     &               rise_fact(a,j)*rise_fact(1+a-b,j)/factorial(j) *
     &               (-1)**j / kp**2

            if (l1) then
               cf(k1,1) = cf(k1,1) + w*bt1_1*kp/(i+j)
               cf(k1,2) = cf(k1,2) + w*bt1_2*kp/(i+j)
            endif

            if (l2) cf(k2,1) = cf(k2,1) - w/(4*kp*(2-i-j))

            if (l3) cf(k3,2) = cf(k3,2) + w/(2*(-i-j+1))
         enddo
      enddo
      

      w1 = - (bt1_1 + gmm/(2*kp**2) + (3*bt0**2)/(2*kp**4)) / kp**2
      w2 = - (bt1_2 - bt0/kp**2) / kp**2
      do k = 1, pmax
         cf(k,1) = cf(k,1) + w1 * binom_cf(-1.0d0, k) * kp**(k+1) / k
         cf(k,2) = cf(k,2) + w2 * binom_cf(-1.0d0, k) * kp**(k+1) / k
      enddo
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl_trapz(n, m, erg, zq, oi)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      double precision oi(2)

      integer k, ngr
      double precision eb, eb_n, de, eta, wt, kp, oi_p1(2), oi_p2(2)
      double precision ee(2)
      parameter (ngr=100000)


      kp = sqrt(2 * abs(erg))
      eb = 1.0d0 / kp           ! so that zb = 1.0d0 inside wfat_calc_o_intl_part* below.
      eb = 0.5d0 / kp
      eb_n = 500.0d0 / kp
      
      call wfat_calc_o_intl_part1(n, m, erg, zq, eb, oi_p1)

      de = (eb_n - eb) / (ngr - 1)
      oi_p2 = 0.0d0
      do k = 1, ngr
         eta = eb + (k-1)*de
         if (k==1 .or. k==ngr) then
            wt = 0.5d0 * de
         else
            wt = 1.0d0 * de
         endif

         call wfat_calc_o_intd(eta, n, m, erg, zq, ee)
         oi_p2 = oi_p2 + wt * ee
      enddo

      oi = oi_p1 + oi_p2

      end
c=======================================================      
