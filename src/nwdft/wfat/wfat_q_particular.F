c=======================================================
c=======================================================
      subroutine wfat_q_particular(nr, r, n, m, lmax, erg, zq, nlgl,
     &                             q_p1, q_p2)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer nr
      double precision r(nr)    ! [input], increasingly ordered radial points
      integer n
      integer m
      integer lmax
      double precision erg
      double precision zq
      integer nlgl
      double precision q_p1(abs(m):lmax, nr)       ! [output], orientation-independent
      double precision q_p2(abs(m):lmax, nr)       ! [output], orientation-dependent

      logical first_r_1
      integer i, j, l, nn, nc, nci, i_gi1, i_gi2, i_hi1, i_hi2, i_ci
      integer i_cd, i_chi, i_chil, i_cgi, pmax, pmax_l
      integer l_pt, l_wt, l_gi, l_hi, l_gg, l_hh, l_ci, l_cd, l_chi
      integer l_chil, l_cgi
      integer k_pt, k_wt, k_gi, k_hi, k_gg, k_hh, k_ci, k_cd, k_chi
      integer k_chil, k_cgi
      double precision kappa, lmb, a, b, gl, f, g, h, r_, y, kp, km
      double precision om_m, om_p, om_e, ybound, rbound
      double precision chgm(abs(m):lmax), chgu(abs(m):lmax)
      double precision qq(abs(m):lmax), factorial
      external factorial
      parameter( pmax=40, ybound=1.0d0 )


cOLD      pmax = 30          !(2*lmax + 3) + 40

      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0
      rbound = ybound / lmb
      do l = abs(m), lmax
         call wfat_gammaf(l + 1 - zq/kappa, gl)
         qq(l) = kappa * 2**(2*l+2) * gl / factorial(2*l+1)
      enddo


c     ==== Obtain the unshifted LGL quadrature nodes and weights ====
      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_pt', l_pt, k_pt))
     &     call errquit('WFAT: lgl_pt allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_wt', l_wt, k_wt))
     &     call errquit('WFAT: lgl_wt allocation failed', 85, MA_ERR)
      call wfat_lgl(-1.0d0, 1.0d0, nlgl-1, dbl_mb(k_pt), dbl_mb(k_wt))

      
      nn = 2 * (lmax - abs(m) + 1)
      if (.not. MA_alloc_get(MT_DBL, nn, 'g_i', l_gi, k_gi))
     &     call errquit('WFAT: g_i allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nn, 'h_i', l_hi, k_hi))
     &     call errquit('WFAT: h_i allocation failed', 85, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nn, 'gg', l_gg, k_gg))
     &     call errquit('WFAT: gg allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nn, 'hh', l_hh, k_hh))
     &     call errquit('WFAT: hh allocation failed', 85, MA_ERR)

      
c     ==== Obtain various power series coefficients (for small r) ====
      nc = pmax * (lmax - abs(m) + 1)
      if (.not. MA_alloc_get(MT_DBL, nc, 'coeff_oi', l_ci, k_ci))
     &     call errquit('WFAT: coeff_oi allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nc, 'coeff_od', l_cd, k_cd))
     &     call errquit('WFAT: coeff_od allocation failed', 85, MA_ERR)
      nci = pmax * 2 * (lmax - abs(m) + 1)
      if (.not. MA_alloc_get(MT_DBL, nci, 'coeff_hi', l_chi, k_chi))
     &     call errquit('WFAT: coeff_hi allocation failed', 85, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nci, 'coeff_hil', l_chil, k_chil))
     &     call errquit('WFAT: coeff_hil allocation failed', 85, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nci, 'coeff_gi', l_cgi, k_cgi))
     &     call errquit('WFAT: coeff_gi allocation failed', 85, MA_ERR)
      do l = abs(m), lmax
         if (l == abs(m)) then
            km = 0.0d0
            om_m = 0.0d0
         else
            km = sqrt(dble(l**2-m**2) / (4*l**2-1))
            call wfat_calc_radnorm(l-1, n, m, int(zq), kappa, om_m)
         endif
         kp = sqrt(dble((l+1)**2-m**2) / (4*(l+1)**2-1))
         call wfat_calc_radnorm(l+1, n, m, int(zq), kappa, om_p)
         call wfat_calc_radnorm(l, n, m, int(zq), kappa, om_e)

c        == Coefficients of A-B ==
         i_ci = k_ci + (l-abs(m))*pmax
         i_cd = k_cd + (l-abs(m))*pmax
         call wfat_q_particular_coeff(l, m, km, kp, om_m, om_p, om_e,
     &        erg, zq, pmax, dbl_mb(i_ci), dbl_mb(i_cd))

c        == Coefficients of h_i ==
         i_chi = k_chi + (l-abs(m))*pmax*2
         i_chil = k_chil + (l-abs(m))*pmax*2
         call wfat_h_i_coeff(l, m, km, kp, om_m, om_p, om_e,
     &        erg, zq, pmax, dbl_mb(i_chi), dbl_mb(i_chil))

c        == Coefficients of g_i ==
         i_cgi = k_cgi + (l-abs(m))*pmax*2
         call wfat_g_i_coeff(l, m, km, kp, om_m, om_p, om_e,
     &        erg, zq, pmax, dbl_mb(i_cgi))
      enddo
      

      first_r_1 = .true.
c     ==== Loop over the increasingly-ordered radial points ====
      do i = 1, nr             
         if (lmb*r(i) < ybound) then
            y = lmb*r(i)

            do l = abs(m), lmax
               q_p1(l, i) = 0.0d0
               q_p2(l, i) = 0.0d0
               i_ci = k_ci + (l-abs(m))*pmax
               i_cd = k_cd + (l-abs(m))*pmax
               do j = 1, pmax
                  q_p1(l, i) = q_p1(l, i) + dbl_mb(i_ci+j-1) * y**j
                  q_p2(l, i) = q_p2(l, i) + dbl_mb(i_cd+j-1) * y**j
               enddo
               f = (y/2)**l * exp(-y/2)
               q_p1(l, i) = q_p1(l, i) * f
               q_p2(l, i) = q_p2(l, i) * f
            enddo
         else
            call wfat_chgm_1pt(r(i), abs(m), lmax, a, b, lmb, chgm)
            call wfat_chgu_1pt(r(i), abs(m), lmax, a, b, lmb, chgu)
            
            
c           ==== Calculate the integrals involving ====
c           ====   g and h within r(i-1) and r(i)  ====
            r_ = r(i-1)
            if (first_r_1) then
               first_r_1 = .false.
               r_ = rbound
               call wfat_calc_gh_i_smallr(rbound, erg, m, lmax, pmax, 
     &              dbl_mb(k_chi), dbl_mb(k_chil), dbl_mb(k_cgi), 
     &              dbl_mb(k_gi), dbl_mb(k_hi))
            endif
            
            
c           ==== Calculate the integrals involving ====
c           ====     g and h within 0 and r(i)     ====
            call wfat_calc_gh_i(r_, r(i), n, m, lmax, erg, zq, nlgl,
     &           dbl_mb(k_pt), dbl_mb(k_wt), dbl_mb(k_gg), dbl_mb(k_hh))
            dbl_mb(k_gi:k_gi+nn-1) = dbl_mb(k_gi:k_gi+nn-1) +
     &                               dbl_mb(k_gg:k_gg+nn-1)
            dbl_mb(k_hi:k_hi+nn-1) = dbl_mb(k_hi:k_hi+nn-1) +
     &                               dbl_mb(k_hh:k_hh+nn-1)
            
            
c           ==== Calculate q_p1 and q_p2 at r(i) ====
            do l = abs(m), lmax
               f = (kappa*r(i))**l * exp(-kappa*r(i))
               g = f * chgm(l)
               h = f * chgu(l)
            
               i_gi1 = k_gi + (l-abs(m))*2
               i_gi2 = i_gi1 + 1
               i_hi1 = k_hi + (l-abs(m))*2
               i_hi2 = i_hi1 + 1
               q_p1(l, i) = qq(l) * (g * dbl_mb(i_hi1) -
     &                      h * dbl_mb(i_gi1))
               q_p2(l, i) = qq(l) * (g * dbl_mb(i_hi2) -
     &                      h * dbl_mb(i_gi2))
            enddo

         endif
      enddo

      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate lgl_pt', 49, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate lgl_wt', 49, MA_ERR)
      if (.not. ma_free_heap(l_gi))
     &     call errquit('WFAT: Cannot deallocate g_i', 49, MA_ERR)
      if (.not. ma_free_heap(l_hi))
     &     call errquit('WFAT: Cannot deallocate h_i', 49, MA_ERR)
      if (.not. ma_free_heap(l_gg))
     &     call errquit('WFAT: Cannot deallocate gg', 49, MA_ERR)
      if (.not. ma_free_heap(l_hh))
     &     call errquit('WFAT: Cannot deallocate hh', 49, MA_ERR)
      if (.not. ma_free_heap(l_ci))
     &     call errquit('WFAT: Cannot deallocate coeff_oi', 49, MA_ERR)
      if (.not. ma_free_heap(l_cd))
     &     call errquit('WFAT: Cannot deallocate coeff_od', 49, MA_ERR)
      if (.not. ma_free_heap(l_chi))
     &     call errquit('WFAT: Cannot deallocate coeff_hi', 49, MA_ERR)
      if (.not. ma_free_heap(l_chil))
     &     call errquit('WFAT: Cannot deallocate coeff_hil', 49, MA_ERR)
      if (.not. ma_free_heap(l_cgi))
     &     call errquit('WFAT: Cannot deallocate coeff_gi', 49, MA_ERR)
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_gh_i(r1, r2, n, m, lmax, erg, zq, nlgl,
     &                          pt, wt, gg, hh)

      implicit none
#include "stdio.fh"
#include "global.fh"      
      double precision r1
      double precision r2
      integer n
      integer m
      integer lmax
      double precision erg
      double precision zq
      integer nlgl
      double precision pt(nlgl)
      double precision wt(nlgl)
      double precision gg(2,abs(m):lmax)
      double precision hh(2,abs(m):lmax)

      integer nfem, ngr, id, i, j, l
      double precision pt_4(4), wt_4(4), pt_6(6), wt_6(6)
      double precision dr, dfem, kappa, lmb, a, b, r1_, r2_, r, wt_, f
      double precision g, h, c, c_m, c_p, k_m, k_p, y, ind_oi, ind_od
      double precision chgm(abs(m):lmax+1), chgu(abs(m):lmax)
      double precision min_dfem
      parameter( min_dfem=2.0d0 )


      pt_4 = (/ -1.0d0,
     &          -0.44721359549996d0,
     &           0.44721359549996d0,
     &           1.0d0 /)
      wt_4 = (/ 0.16666666666667d0,
     &          0.83333333333333d0,
     &          0.83333333333333d0,
     &          0.16666666666667d0 /)
      pt_6 = (/ -1.0d0,
     &          -0.76505532392947d0,
     &          -0.28523151648065d0,
     &           0.28523151648065d0,
     &           0.76505532392947d0,
     &           1.0d0 /)
      wt_6 = (/ 0.06666666666667d0,
     &          0.37847495629785d0,
     &          0.55485837703549d0,
     &          0.55485837703549d0,
     &          0.37847495629785d0,
     &          0.06666666666667d0 /)

      
c     ==== Determine the no. of finite elements and grids ====
c     dr = the total interval of r to be integrated
c     dfem = the finite element interval
      dr = r2 - r1
      if (dr < 1.0d-12) then
         nfem = 1
         dfem = dr
         ngr = 2
         id = 0
      else if (dr >= 1.0d-12 .and. dr < 1.0d-3) then
         nfem = 1
         dfem = dr
         ngr = 4
         id = 1
      else if (dr >= 1.0d-3 .and. dr < 1.0d-1) then
         nfem = 1
         dfem = dr
         ngr = 6
         id = 2
      else
         nfem = int(max(dr,min_dfem+1.0d-13) / min_dfem)
         dfem = dr / nfem
         ngr = nlgl
         id = 3
      endif

cfordebug      if (ga_nodeid()==0)
cfordebug     &     write(luout, '("nfem, ngr, id, r1, r2, dr : ", '//
cfordebug     &     '3i4, 3es21.12)') nfem, ngr, id, r1, r2, dr
      
      
      
      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0

      gg = 0.0d0
      hh = 0.0d0
c     ==== Loop over finite elements ====
c     r1_ = lower bound of a FEM
c     r2_ = upper bound of a FEM
c     r = shifted quadrature grid
c     wt_ = shifted quadrature weight
      do i = 1, nfem
         r1_ = r1 + (i-1)*dfem
         r2_ = r1_ + dfem

c        ==== Loop over quadrature grids ====
         do j = 1, ngr
c           ==== Determine the shifted grids and weights ====
            if (id == 0) then
               r = r1_ + (j-1)*dfem
               wt_ = (r2_ - r1_)/2
            else if (id == 1) then
               r = 0.5d0*((r2_-r1_)*pt_4(j) + r2_ + r1_)
               wt_ = (r2_ - r1_)/2 * wt_4(j)
            else if (id == 2) then
               r = 0.5d0*((r2_-r1_)*pt_6(j) + r2_ + r1_)
               wt_ = (r2_ - r1_)/2 * wt_6(j)
            else if (id == 3) then
               r = 0.5d0*((r2_-r1_)*pt(j) + r2_ + r1_)
               wt_ = (r2_ - r1_)/2 * wt(j)
            endif
            if (r < 1.0d-13) cycle           ! When r is zero, the integrands are equal to zero.
            

c           ==== Calculate the necessary Hypergeometric functions ====
            call wfat_chgm_1pt(r, abs(m), lmax+1, a, b, lmb, chgm)     ! 1)
            call wfat_chgu_1pt(r, abs(m), lmax, a, b, lmb, chgu)     

            do l = abs(m), lmax
               f = (kappa*r)**l * exp(-kappa*r)
               g = f * chgm(l)
               h = f * chgu(l)

               call wfat_calc_radnorm(l, n, m, int(zq), kappa, c)
               call wfat_calc_radnorm(l+1, n, m, int(zq), kappa, c_p)
               k_p = sqrt(dble((l+1)**2-m**2) / (4*(l+1)**2-1))
               y = k_p * c_p*(kappa*r)**(l+1)*exp(-kappa*r) * chgm(l+1)
               if (l > abs(m)) then
                  call wfat_calc_radnorm(l-1, n, m, int(zq), kappa, c_m)
                  k_m = sqrt(dble(l**2-m**2) / (4*l**2-1))
                  y = y + k_m *
     &                c_m*(kappa*r)**(l-1)*exp(-kappa*r) * chgm(l-1)
               endif

c              ==== The g_l^\nu part ====
               ind_oi = g * r**3 * y
               ind_od = g * r**2 * c*g
               gg(1,l) = gg(1,l) + wt_*ind_oi
               gg(2,l) = gg(2,l) + wt_*ind_od

c              ==== The h_l^\nu part ====
               ind_oi = h * r**3 * y
               ind_od = h * r**2 * c*g
               hh(1,l) = hh(1,l) + wt_*ind_oi
               hh(2,l) = hh(2,l) + wt_*ind_od
            enddo

         enddo
      enddo
      

      end
c     NOTES:
c     1) lmax+1 instead of just lmax because the integrands involve the evaluation
c        of R_{l+1}^\nu, see below.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_q_particular_coeff(l, m, km, kp, om_m, om_p, om_e,
     &           erg, zq, pmax, c_oi, c_od)
c
c     Evaluate the expansion coefficients of eq. (306).
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer l                       ! [input],
      integer m                       ! [input],
      double precision km             ! [input],
      double precision kp             ! [input],
      double precision om_m           ! [input],
      double precision om_p           ! [input],
      double precision om_e           ! [input],
      double precision erg            ! [input],
      double precision zq             ! [input],
      integer pmax                    ! [input],
      double precision c_oi(1:pmax)   ! [output],
      double precision c_od(1:pmax)   ! [output],

      logical lg_m, lg_p, lg_e
      integer l_t_, l_e_, l_t, l_e
      integer k_t_, k_e_, k_t, k_e
      integer ubound, i, j, k, k_, s2, ii, jj, kk, sg, k3, s1
      double precision kappa, a, b, ga, w_, mk, tm_, tp_, u_, w, dk, tm
      double precision tp, u, gl, ff, ps1, ps2, ps3
      double precision factorial, rise_fact
      external factorial, rise_fact


      if (l < abs(m))
     &     call errquit('WFAT: Argument l of wfat_q_particular_coeff '//
     &     'cannot be smaller than |m|.', 48, UERR)
      if (pmax < 2)
     &     call errquit('WFAT: Argument pmax of '//
     &     'wfat_q_particular_coeff cannot be smaller than 2.',
     &     48, UERR)

      
      kappa = sqrt(2*abs(erg))
      a = l + 1 - zq/kappa
      b = dble(2*l + 2)
      call wfat_gammaf(a, ga)


      if (.not. MA_alloc_get(MT_DBL, pmax, 'fin_k_oi', l_t_, k_t_))
     &     call errquit('WFAT: fin_k_oi allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, pmax, 'fin_k_od', l_e_, k_e_))
     &     call errquit('WFAT: fin_k_od allocation failed', 85, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, pmax, 'inf_k_oi', l_t, k_t))
     &     call errquit('WFAT: inf_k_oi allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, pmax, 'inf_k_od', l_e, k_e))
     &     call errquit('WFAT: inf_k_od allocation failed', 85, MA_ERR)
      
      
c     ==== Finite k part ====
      dbl_mb(k_t_:k_t_+pmax-1) = 0.0d0
      dbl_mb(k_e_:k_e_+pmax-1) = 0.0d0
      ubound = pmax - 2
      do i = 0, ubound
         do j = 0, ubound
            do k = 2*l+1, 1, -1
               do k_ = 0, ubound
                  s2 = i + j - k + k_
                  ii = 2*l + 3 + s2
                  jj = 2*l + 5 + s2
                  kk = 2*l + 3 + s2
                  lg_m = (ii >= 2 .and. ii <= pmax .and. l > abs(m))
                  lg_p = (jj >= 4 .and. jj <= pmax)
                  lg_e = (kk >= 2 .and. kk <= pmax)

                  if (lg_m .or. lg_p .or. lg_e) then
                     w_ = ( (-1)**i/factorial(i) )                ! eq. (307d)
     &                    *( 1.0d0/factorial(j) )
     &                    *( factorial(k-1)*rise_fact(1-a+k,2*l+1-k)
     &                    /factorial(2*l+1-k) )
     &                    *( rise_fact(a,k_)/rise_fact(b,k_)
     &                    /factorial(k_) ) * (k_+k)/(2*kappa)**4
                  endif

                  if (lg_m) then
                     sg = 2*l + 3 + i + j
                     mk = rise_fact(a-1,j) / rise_fact(b-2,j) /
     &                    (sg-k) / (sg+k_)
                     tm_ = km * om_m / 2**(2*l-1) * mk            ! eq. (307e)
                     k3 = k_t_ + ii - 1
                     dbl_mb(k3) = dbl_mb(k3) + w_ * tm_
                  endif

                  if (lg_p) then
                     sg = 2*l + 5 + i + j
                     mk = rise_fact(a+1,j) / rise_fact(b+2,j) /
     &                    (sg-k) / (sg+k_)
                     tp_ = kp * om_p / 2**(2*l+1) * mk            ! eq. (307e)
                     k3 = k_t_ + jj - 1
                     dbl_mb(k3) = dbl_mb(k3) + w_ * tp_
                  endif

                  if (lg_e) then
                     sg = 2*l + 3 + i + j
                     mk = rise_fact(a,j) / rise_fact(b,j) /
     &                    (sg-k) / (sg+k_)
                     u_ = 2*kappa * om_e / 2**(2*l) * mk          ! eq. (307f)
                     k3 = k_e_ + kk - 1
                     dbl_mb(k3) = dbl_mb(k3) + w_ * u_
                  endif

               enddo
            enddo
         enddo
      enddo
      dbl_mb(k_t_:k_t_+pmax-1) = dbl_mb(k_t_:k_t_+pmax-1) / ga
      dbl_mb(k_e_:k_e_+pmax-1) = dbl_mb(k_e_:k_e_+pmax-1) / ga


c     ==== Infinite k part ====
      ubound = pmax - 2*l - 3
      dbl_mb(k_t:k_t+pmax-1) = 0.0d0
      dbl_mb(k_e:k_e+pmax-1) = 0.0d0
      do i = 0, ubound
         do j = 0, ubound
            do k = 0, ubound
               do k_ = 0, ubound
                  s1 = i + j + k + k_
                  ii = 2*l + 3 + s1
                  jj = 2*l + 5 + s1
                  kk = 2*l + 3 + s1
                  lg_m = (ii >= 2*l+3 .and. ii <= pmax .and. l > abs(m))
                  lg_p = (jj >= 2*l+5 .and. jj <= pmax)
                  lg_e = (kk >= 2*l+3 .and. kk <= pmax)

                  if (lg_m .or. lg_p .or. lg_e) then
                     w = ( (-1)**i/factorial(i) )                ! eq. (307a)
     &                   *( 1.0d0/factorial(j) )
     &                   *( rise_fact(a,k)/rise_fact(b,k)/factorial(k) )
     &                   *( rise_fact(a,k_)/rise_fact(b,k_)
     &                   /factorial(k_) ) / (2*kappa)**4
                     call wfat_psif(k+a, ps1)
                     call wfat_psif(dble(k+1), ps2)
                     call wfat_psif(dble(k+b), ps3)
                     dk = ps1 - ps2 - ps3
                  endif
                  
                  if (lg_m) then
                     sg = 2*l + 3 + i + j
                     mk = rise_fact(a-1,j)/rise_fact(b-2,j)
     &                    *(dk*dble(k_-k)/(sg+k)/(sg+k_) -
     &                    1.0d0/(sg+k)**2)
                     tm = km * om_m / 2**(2*l-1) * mk            ! eq. (307b)
                     k3 = k_t + ii - 1
                     dbl_mb(k3) = dbl_mb(k3) + w * tm
                  endif
                  
                  if (lg_p) then
                     sg = 2*l + 5 + i + j
                     mk = rise_fact(a+1,j)/rise_fact(b+2,j)
     &                    *(dk*dble(k_-k)/(sg+k)/(sg+k_) -
     &                    1.0d0/(sg+k)**2)
                     tp = kp * om_p / 2**(2*l+1) * mk            ! eq. (307b)
                     k3 = k_t + jj - 1
                     dbl_mb(k3) = dbl_mb(k3) + w * tp
                  endif

                  if (lg_e) then
                     sg = 2*l + 3 + i + j
                     mk = rise_fact(a,j)/rise_fact(b,j)
     &                    *(dk*dble(k_-k)/(sg+k)/(sg+k_) -
     &                    1.0d0/(sg+k)**2)
                     u = 2*kappa * om_e / 2**(2*l) * mk          ! eq. (307c)
                     k3 = k_e + kk - 1
                     dbl_mb(k3) = dbl_mb(k3) + w * u
                  endif
                  
               enddo
            enddo
         enddo
      enddo
      call wfat_gammaf(-l-zq/kappa, gl)
      dbl_mb(k_t:k_t+pmax-1) = dbl_mb(k_t:k_t+pmax-1)
     &                         / factorial(2*l+1) / gl
      dbl_mb(k_e:k_e+pmax-1) = dbl_mb(k_e:k_e+pmax-1)
     &                         / factorial(2*l+1) / gl


c     ==== The final value of the coefficients ====
      ff = kappa * 2**(2*l+2) * ga / factorial(2*l+1)
      c_oi = ff * (dbl_mb(k_t:k_t+pmax-1) + dbl_mb(k_t_:k_t_+pmax-1))
      c_od = ff * (dbl_mb(k_e:k_e+pmax-1) + dbl_mb(k_e_:k_e_+pmax-1))


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_t_))
     &     call errquit('WFAT: Cannot deallocate fin_k_oi', 49, MA_ERR)
      if (.not. ma_free_heap(l_e_))
     &     call errquit('WFAT: Cannot deallocate fin_k_od', 49, MA_ERR)
      if (.not. ma_free_heap(l_t))
     &     call errquit('WFAT: Cannot deallocate inf_k_oi', 49, MA_ERR)
      if (.not. ma_free_heap(l_e))
     &     call errquit('WFAT: Cannot deallocate inf_k_od', 49, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_h_i_coeff(l, m, km, kp, om_m, om_p, om_e, erg, zq,
     &           pmax, hi, hil)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      integer l
      integer m                       ! [input],
      double precision km             ! [input],
      double precision kp             ! [input],
      double precision om_m           ! [input],
      double precision om_p           ! [input],
      double precision om_e           ! [input],
      double precision erg            ! [input],
      double precision zq             ! [input],
      integer pmax                    ! [input],
      double precision hi(pmax, 2)    ! 1st row: orientation-independent
                                      ! 2nd row: orientation-dependent
      double precision hil(pmax, 2)   ! 1st row: ln-multiplied orientation-independent
                                      ! 2nd row: ln-multiplied orientation-dependent

      logical lg_m, lg_p, lg_e
      integer i, j, k, ii, jj, kk, ubound
      double precision kappa, a, b, gl, ga, v1, v2, w_, mk, t_m_, t_p_
      double precision u_, w, ps1, ps2, ps3, dk, tl_m, t_m, tl_p, t_p
      double precision ul, u, factorial, rise_fact
      external factorial, rise_fact
      

      if (l < abs(m))
     &     call errquit('WFAT: Argument l of wfat_h_i_coeff cannot be'//
     &     'smaller than |m|.', 48, UERR)
      if (pmax < 2)
     &     call errquit('WFAT: Argument pmax of wfat_h_i_coeff '//
     &     'cannot be smaller than 2.', 48, UERR)
      
      
      kappa = sqrt(2*abs(erg))
      a = l + 1 - zq/kappa
      b = dble(2*l + 2)
      call wfat_gammaf(-l-zq/kappa, gl)
      v1 = 1.0d0 / factorial(2*l+1) / gl
      call wfat_gammaf(a, ga)
      v2 = 1.0d0 / ga

      
      hi = 0.0d0
      hil = 0.0d0
      ubound = pmax - 2
      do i = 0, ubound
         do j = 0, ubound
            do k = 2*l+1, 1, -1
               ii = 2*l + 3 + i + j - k
               jj = 2*l + 5 + i + j - k
               kk = 2*l + 3 + i + j - k
               lg_m = (ii >= 2 .and. ii <= pmax .and. l > abs(m))
               lg_p = (jj >= 4 .and. jj <= pmax)
               lg_e = (kk >= 2 .and. kk <= pmax)

               if (lg_m .or. lg_p .or. lg_e) then
                  w_ = v2 * ( (-1)**i/factorial(i) )
     &                 *( 1.0d0/factorial(j) )
     &                 *( factorial(k-1)*rise_fact(1-a+k,2*l+1-k)
     &                 /factorial(2*l+1-k) )
     &                 /(2*kappa)**4
               endif

               if (lg_m) then
                  mk = rise_fact(a-1,j) / rise_fact(b-2,j) / ii
                  t_m_ = km * om_m / 2**(2*l-1) * mk
                  hi(ii,1) = hi(ii,1) + w_ * t_m_
               endif

               if (lg_p) then
                  mk = rise_fact(a+1,j) / rise_fact(b+2,j) / jj
                  t_p_ = kp * om_p / 2**(2*l+1) * mk
                  hi(jj,1) = hi(jj,1) + w_ * t_p_
               endif

               if (lg_e) then
                  mk = rise_fact(a,j) / rise_fact(b,j) / kk
                  u_ = 2*kappa * om_e / 2**(2*l) * mk
                  hi(kk,2) = hi(kk,2) + w_ * u_
               endif
               
            enddo
         enddo
      enddo
      

      ubound = pmax - 2*l - 3
      do i = 0, ubound
         do j = 0, ubound
            do k = 0, ubound
               ii = 2*l + 3 + i + j + k
               jj = 2*l + 5 + i + j + k
               kk = 2*l + 3 + i + j + k
               lg_m = (ii >= 2*l+3 .and. ii <= pmax .and. l > abs(m))
               lg_p = (jj >= 2*l+5 .and. jj <= pmax)
               lg_e = (kk >= 2*l+3 .and. kk <= pmax)
                  
               if (lg_m .or. lg_p .or. lg_e) then
                  w = v1 * ( (-1)**i/factorial(i) )
     &                *( 1.0d0/factorial(j) )
     &                *( rise_fact(a,k)/rise_fact(b,k)/factorial(k) )
     &                /(2*kappa)**4
                  call wfat_psif(k+a, ps1)
                  call wfat_psif(dble(k+1), ps2)
                  call wfat_psif(dble(k+b), ps3)
                  dk = ps1 - ps2 - ps3
               endif

               if (lg_m) then
c                 == Terms with ln ==
                  mk = rise_fact(a-1,j) / rise_fact(b-2,j) / ii
                  tl_m = km * om_m / 2**(2*l-1) * mk
                  hil(ii, 1) = hil(ii, 1) + w * tl_m

c                 == Terms without ln ==
                  mk = rise_fact(a-1,j) / rise_fact(b-2,j) *
     &                 (dk/ii - 1.0d0/ii**2)
                  t_m = km * om_m / 2**(2*l-1) * mk
                  hi(ii, 1) = hi(ii, 1) + w * t_m
               endif

               if (lg_p) then
c                 == Terms with ln ==
                  mk = rise_fact(a+1,j) / rise_fact(b+2,j) / jj
                  tl_p = kp * om_p / 2**(2*l+1) * mk
                  hil(jj, 1) = hil(jj, 1) + w * tl_p

c                 == Terms without ln ==
                  mk = rise_fact(a+1,j) / rise_fact(b+2,j) *
     &                 (dk/jj - 1.0d0/jj**2)
                  t_p = kp * om_p / 2**(2*l+1) * mk
                  hi(jj, 1) = hi(jj, 1) + w * t_p
               endif

               if (lg_e) then
c                 == Terms with ln ==
                  mk = rise_fact(a,j) / rise_fact(b,j) / kk
                  ul = 2*kappa * om_e / 2**(2*l) * mk
                  hil(kk, 2) = hil(kk, 2) + w * ul

c                 == Terms without ln ==
                  mk = rise_fact(a,j) / rise_fact(b,j) *
     &                 (dk/kk - 1.0d0/kk**2)
                  u = 2*kappa * om_e / 2**(2*l) * mk
                  hi(kk, 2) = hi(kk, 2) + w * u
               endif
               
            enddo
         enddo
      enddo
      

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_g_i_coeff(l, m, km, kp, om_m, om_p, om_e, erg, zq,
     &           pmax, gi)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      integer l
      integer m                       ! [input],
      double precision km             ! [input],
      double precision kp             ! [input],
      double precision om_m           ! [input],
      double precision om_p           ! [input],
      double precision om_e           ! [input],
      double precision erg            ! [input],
      double precision zq             ! [input],
      integer pmax                    ! [input],
      double precision gi(pmax, 2)      ! 1st row: orientation-independent
                                        ! 2nd row: orientation-dependent

      logical lg_m, lg_p, lg_e
      integer i, j, k_, ii, jj, kk, ubound
      double precision kappa, a, b, mk, w, t_m, t_p, u
      double precision factorial, rise_fact
      external factorial, rise_fact
      

      if (l < abs(m))
     &     call errquit('WFAT: Argument l of wfat_h_i_coeff cannot be'//
     &     'smaller than |m|.', 48, UERR)
cc      if (pmax < 2*l + 3)                  !! THIS CONDITIONAL SEEMS TO BE SAFE TO IGNORE. IF THIS CONDITIONAL IS FALSE THEN THE OUTPUT gi CONTAINS JUST ZEROS.
cc     &     call errquit('WFAT: Argument pmax of wfat_g_i_coeff '//
cc     &     'cannot be smaller than 2l+3.', 48, UERR)
      
      
      kappa = sqrt(2*abs(erg))
      a = l + 1 - zq/kappa
      b = dble(2*l + 2)
      

      gi = 0.0d0
      ubound = pmax - 2*l - 3
      do i = 0, ubound
         do j = 0, ubound
            do k_ = 0, ubound
               ii = 2*l + 3 + i + j + k_
               jj = 2*l + 5 + i + j + k_
               kk = 2*l + 3 + i + j + k_
               lg_m = (ii >= 2*l+3 .and. ii <= pmax .and. l > abs(m))
               lg_p = (jj >= 2*l+5 .and. jj <= pmax)
               lg_e = (kk >= 2*l+3 .and. kk <= pmax)

               if (lg_m .or. lg_p .or. lg_e) then
                  w = ( (-1)**i/factorial(i) )
     &                *( 1.0d0/factorial(j) )
     &                *( rise_fact(a,k_)/rise_fact(b,k_)/factorial(k_) )
     &                /(2*kappa)**4
               endif

               if (lg_m) then
                  mk = rise_fact(a-1,j) / rise_fact(b-2,j) / ii
                  t_m = km * om_m / 2**(2*l-1) * mk
                  gi(ii, 1) = gi(ii, 1) + w * t_m
               endif

               if (lg_p) then
                  mk = rise_fact(a+1,j) / rise_fact(b+2,j) / jj
                  t_p = kp * om_p / 2**(2*l+1) * mk
                  gi(jj, 1) = gi(jj, 1) + w * t_p
               endif

               if (lg_e) then
                  mk = rise_fact(a,j) / rise_fact(b,j) / kk
                  u = 2*kappa * om_e / 2**(2*l) * mk
                  gi(kk, 2) = gi(kk, 2) + w * u
               endif

            enddo
         enddo
      enddo
      

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_gh_i_smallr(r, erg, m, lmax, pmax, chi, chil,
     &           cgi, gi, hi)

      implicit none
      double precision r                           ! [input],
      double precision erg                         ! [input],
      integer m                                    ! [input],
      integer lmax                                 ! [input],
      integer pmax                                 ! [input],
      double precision chi(pmax, 2, abs(m):lmax)   ! [input],
      double precision chil(pmax, 2, abs(m):lmax)  ! [input],
      double precision cgi(pmax, 2, abs(m):lmax)   ! [input],
      double precision gi(2, abs(m):lmax)          ! [output],
      double precision hi(2, abs(m):lmax)          ! [output],

      integer l, i, j
      double precision kappa, y, ln_y


      kappa = sqrt(2 * abs(erg))
      y = 2*kappa*r
      ln_y = log(y)

      gi = 0.0d0
      hi = 0.0d0
      do l = abs(m), lmax
         do i = 1, 2
            do j = 1, pmax
               gi(i,l) = gi(i,l) + cgi(j,i,l) * y**j
               hi(i,l) = hi(i,l) +
     &                   (ln_y*chil(j,i,l) + chi(j,i,l)) * y**j
            enddo
         enddo
      enddo

      end
c=======================================================
      
