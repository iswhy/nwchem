c=======================================================
c=======================================================
      subroutine wfat_r1(neta, eta, n, m, erg, zq, bt, gm, chmmt, nlgl,
     &           r1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer neta
      double precision eta(neta)    ! [input], increasingly ordered eta.
      integer n
      integer m
      double precision erg
      double precision zq
      double precision bt           ! [input], beta orientation angle in degree.
      double precision gm           ! [input], gamma orientation angle in degree.
      double precision chmmt(3)
      integer nlgl
      double precision r1(neta)

      logical first_eta_1
      integer i, j, pmax, md
      integer l_pt, l_wt, l_cr, l_coi, l_coil, l_cri
      integer k_pt, k_wt, k_cr, k_coi, k_coil, k_cri
      double precision protmat(3,3), mu(3), kappa, a, b, g, d, z, q
      double precision eta_, ebound, ri, oi, rr, oo, r0, o0, zbound
      double precision chgu, factorial
      double complex zz, chgm_
      parameter (zbound=1.0d0, pmax=40)
      external factorial
      

c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      mu = matmul(transpose(protmat), -chmmt)
      
      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      ebound = zbound / kappa
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))


      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_pt', l_pt, k_pt))
     &     call errquit('WFAT: lgl_pt allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_wt', l_wt, k_wt))
     &     call errquit('WFAT: lgl_wt allocation failed', 85, MA_ERR)
      call wfat_lgl(-1.0d0, 1.0d0, nlgl-1, dbl_mb(k_pt), dbl_mb(k_wt))

      
c     ==== Coefficients of A1-A2 ====
      if (.not. MA_alloc_get(MT_DBL, pmax, 'r1_coef', l_cr, k_cr))
     &     call errquit('WFAT: r1_coef allocation failed', 62, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, pmax, 'o0i_coef', l_coi, k_coi))
     &     call errquit('WFAT: o0i_coef allocation failed', 62, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, pmax, 'o0il_coef', l_coil, k_coil))
     &     call errquit('WFAT: o0il_coef allocation failed', 62, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, pmax, 'r0i_coef', l_cri, k_cri))
     &     call errquit('WFAT: r0i_coef allocation failed', 62, MA_ERR)
      call wfat_r1_coeff(n, m, erg, zq, mu(3), pmax, dbl_mb(k_cr))
      call wfat_o0_i_coeff(n, m, erg, zq, mu(3), pmax, dbl_mb(k_coi),
     &     dbl_mb(k_coil))
      call wfat_r0_i_coeff(n, m, erg, zq, mu(3), pmax, dbl_mb(k_cri))
      
      
      first_eta_1 = .true.
c     ==== Loop over the increasingly-ordered radial points ====
      do i = 1, neta
         z = kappa * eta(i)
         if (z < zbound) then
            q = 0.0d0
            do j = 1, pmax
               q = q + dbl_mb(k_cr+j-1) * z**j
            enddo
            r1(i) = z**(b/2) * exp(-z/2) * q
         else
c           ==== Calculate the integrals involving R^{(0)} ====
c           ====   and O^{(0)} within eta(i-1) and eta(i)  ====
            eta_ = eta(i-1)
            if (first_eta_1) then
               first_eta_1 = .false.
               eta_ = ebound
               call wfat_calc_r0o0_i_smalle(ebound, erg, pmax,
     &              dbl_mb(k_cri), dbl_mb(k_coi), dbl_mb(k_coil),
     &              ri, oi)
            endif

            call wfat_calc_r0o0_i(eta_, eta(i), n, m, erg, zq, mu(3),
     &           nlgl, dbl_mb(k_pt), dbl_mb(k_wt), rr, oo)
            ri = ri + rr
            oi = oi + oo

            zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
            call wfat_cchg(a, b, zz, chgm_)
            r0 = d * z**(b/2) * exp(-z/2) * real(chgm_)
            call wfat_chgu(a, b, z, chgu, md)
            o0 = z**(b/2) * exp(-z/2) * chgu

            r1(i) = o0*ri - r0*oi
         endif
      enddo

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate lgl_pt', 49, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate lgl_wt', 49, MA_ERR)
      if (.not. ma_free_heap(l_cr))
     &     call errquit('WFAT: Cannot deallocate r1_coef', 49, MA_ERR)
      if (.not. ma_free_heap(l_coi))
     &     call errquit('WFAT: Cannot deallocate o0i_coef', 49, MA_ERR)
      if (.not. ma_free_heap(l_coil))
     &     call errquit('WFAT: Cannot deallocate o0il_coef', 49, MA_ERR)
      if (.not. ma_free_heap(l_cri))
     &     call errquit('WFAT: Cannot deallocate cri_coef', 49, MA_ERR)
      
      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_r1_coeff(n, m, erg, zq, mu, pmax, cr1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      double precision mu            ! [input]
      integer pmax                   ! [input]
      double precision cr1(1:pmax)   ! [output]

      logical l_i, l_j, l_k
      integer l_a, k_a, l_b, k_b
      integer i, j, k, k_, ii, jj, kk, k3, s1, s2, ubound, sg
      double precision kappa, a, b, z, beta1, ga, gm, mk, w_, t_, u_, v_
      double precision w, t, u, v, dk, ps1, ps2, ps3
      double precision factorial, rise_fact
      external factorial, rise_fact


      if (pmax < 1) call errquit('Argument pmax of '//
     &     'wfat_r1_smalleta must be larger than 1.', 28, UERR)
      if (n < 0) call errquit('Argument n of '//
     &     'wfat_r1_smalleta must be zero or positive.', 28, UERR)

      cr1 = 0.0d0
      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, ga)


      if (.not. MA_alloc_get(MT_DBL, pmax, 'part_a', l_a, k_a))
     &     call errquit('WFAT: part_a allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, pmax, 'part_b', l_b, k_b))
     &     call errquit('WFAT: part_b allocation failed', 85, MA_ERR)
      

c     ==== Finite k part ====
      dbl_mb(k_b:k_b+pmax-1) = 0.0d0
      if (abs(m) > 0) then
         ubound = pmax - 1
         do i = 0, ubound
            do j = 0, ubound
               do k = abs(m), 1, -1
                  do k_ = 0, ubound
                     s2 = i + j - k + k_
                     ii = abs(m) + s2 + 1
                     jj = abs(m) + s2 + 2
                     kk = abs(m) + s2 + 3
                     l_i = (ii >= 1 .and. ii <= pmax)
                     l_j = (jj >= 2 .and. jj <= pmax)
                     l_k = (kk >= 3 .and. kk <= pmax)

                     if (l_i .or. l_j .or. l_k) then
                        w_ = (-1)**i/factorial(i)
     &                       *( rise_fact(a,j)/factorial(j)
     &                       /rise_fact(b,j) )
     &                       *( rise_fact(a,k_)/factorial(k_)
     &                       /rise_fact(b,k_) )
     &                       *( factorial(k-1)*rise_fact(1-a+k,abs(m)-k)
     &                       /factorial(abs(m)-k) )
     &                       *(k+k_) 
                     endif

                     if (l_i) then
                        sg = abs(m) + 1 + i + j
                        mk = 1.0d0 / (sg+k_) / (sg-k)
                        t_ = w_ * beta1 * kappa * mk
                        k3 = k_b + ii - 1
                        dbl_mb(k3) = dbl_mb(k3) + t_
                     endif

                     if (l_j) then
                        sg = abs(m) + 2 + i + j
                        mk = 1.0d0 / (sg+k_) / (sg-k)
                        u_ = -w_ * mu/2 * mk
                        k3 = k_b + jj - 1
                        dbl_mb(k3) = dbl_mb(k3) + u_
                     endif

                     if (l_k) then
                        sg = abs(m) + 3 + i + j
                        mk = 1.0d0 / (sg+k_) / (sg-k)
                        v_ = w_ / (4*kappa) * mk
                        k3 = k_b + kk - 1
                        dbl_mb(k3) = dbl_mb(k3) + v_
                     endif

                  enddo
               enddo
            enddo
         enddo
         dbl_mb(k_b:k_b+pmax-1) = dbl_mb(k_b:k_b+pmax-1) / ga
      endif
      

c     ==== Infinite k part ====
      ubound = pmax - abs(m) - 1
      dbl_mb(k_a:k_a+pmax-1) = 0.0d0
      do i = 0, ubound
         do j = 0, ubound
            do k = 0, ubound
               do k_ = 0, ubound
                  s1 = i + j + k + k_
                  ii = abs(m) + s1 + 1
                  jj = abs(m) + s1 + 2
                  kk = abs(m) + s1 + 3
                  l_i = (ii >= abs(m)+1 .and. ii <= pmax)
                  l_j = (jj >= abs(m)+2 .and. jj <= pmax)
                  l_k = (kk >= abs(m)+3 .and. kk <= pmax)

                  if (l_i .or. l_j .or. l_k) then
                     w = (-1)**i/factorial(i)
     &                    *( rise_fact(a,j)/rise_fact(b,j)
     &                    /factorial(j) )
     &                    *( rise_fact(a,k)/rise_fact(b,k)
     &                    /factorial(k) )
     &                    *( rise_fact(a,k_)/rise_fact(b,k_)
     &                    /factorial(k_) )
                     call wfat_psif(k+a, ps1)
                     call wfat_psif(dble(k+1), ps2)
                     call wfat_psif(dble(k+1+abs(m)), ps3)
                     dk = ps1 - ps2 - ps3
                  endif
                  
                  if (l_i) then
                     sg = abs(m) + 1 + i + j
                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
                     t = w * beta1 * kappa * mk
                     k3 = k_a + ii - 1
                     dbl_mb(k3) = dbl_mb(k3) + t
                  endif
                  
                  if (l_j) then
                     sg = abs(m) + 2 + i + j
                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
                     u = -w * mu/2 * mk
                     k3 = k_a + jj - 1
                     dbl_mb(k3) = dbl_mb(k3) + u
                  endif

                  if (l_k) then
                     sg = abs(m) + 3 + i + j
                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
                     v = w / (4*kappa) * mk
                     k3 = k_a + kk - 1
                     dbl_mb(k3) = dbl_mb(k3) + v
                  endif
                  
               enddo
            enddo
         enddo
      enddo
      call wfat_gammaf(a-abs(m), gm)
      dbl_mb(k_a:k_a+pmax-1) = dbl_mb(k_a:k_a+pmax-1)
     &     * (-1)**(abs(m)+1) / factorial(abs(m)) / gm


c     ==== The final value of R_\nu^{(1)}(eta) ====
      cr1 = (ga/kappa/factorial(abs(m)))**2 * 
     &     (dbl_mb(k_a:k_a+pmax-1) - dbl_mb(k_b:k_b+pmax-1))


cfordebug      write(luout, '("R1 coefficients in wfat_r1")')
cfordebug      do j = 1, pmax
cfordebug         write(luout, '(i3, 3es15.6)') j, cr1(j),
cfordebug     &        dbl_mb(k_a+j-1), dbl_mb(k_b+j-1)
cfordebug      enddo

      
c     z**(b/2) * exp(-z/2) *

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_a))
     &     call errquit('WFAT: Cannot deallocate part_a', 49, MA_ERR)
      if (.not. ma_free_heap(l_b))
     &     call errquit('WFAT: Cannot deallocate part_b', 49, MA_ERR)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_o0_i_coeff(n, m, erg, zq, mu, pmax, o0i, o0il)
c     Calculate the power series coefficients of 
c        A2(\eta) = \frac{1}{\kappa} \int_0^\eta d\eta' O_\nu^{(0)}(\eta')
c                   (\beta_\nu^{(1)}/\eta' - \mu_z/2 + \eta'/4) R_\nu^{(0)}(\eta')
      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      double precision mu            ! [input]
      integer pmax                   ! [input]
      double precision o0i(1:pmax)   ! [output]
      double precision o0il(1:pmax)  ! [output]

      logical l_i, l_j, l_k
      integer ubound, i, j, k, ii, jj, kk, s1, s2
      double precision kappa, a, b, beta1, v1, v2, ga, gm, w_, t_, u_
      double precision v_, w, ps1, ps2, ps3, dk, mk, t, u, v
      double precision factorial, rise_fact
      external factorial, rise_fact


      if (pmax < 1)
     &     call errquit('WFAT: Argument pmax of wfat_o0_i_coeff '//
     &     'cannot be smaller than 1.', 48, UERR)
      if (n < 0) call errquit('Argument n of '//
     &     'wfat_o0_i_coeff must be zero or positive.', 28, UERR)


      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)

      v1 = 1.0d0 / kappa**2 / factorial(abs(m))
      call wfat_gammaf(a, ga)
      call wfat_gammaf(a-abs(m), gm)
      v2 = (-1)**(abs(m)+1) * ga / (kappa*factorial(abs(m)))**2 / gm

      
c     ==== Finite k part ====
      o0i = 0.0d0
      o0il = 0.0d0
      if (abs(m) > 0) then
         ubound = pmax - 1
         do i = 0, ubound
            do j = 0, ubound
               do k = abs(m), 1, -1
                  s2 = i + j - k
                  ii = abs(m) + s2 + 1
                  jj = abs(m) + s2 + 2
                  kk = abs(m) + s2 + 3
                  l_i = (ii >= 1 .and. ii <= pmax)
                  l_j = (jj >= 2 .and. jj <= pmax)
                  l_k = (kk >= 3 .and. kk <= pmax)

                  if (l_i .or. l_j .or. l_k) then
                     w_ = v1 * (-1)**i/factorial(i)
     &                    *( rise_fact(a,j)/factorial(j)
     &                    /rise_fact(b,j) )
     &                    *( factorial(k-1)*rise_fact(1-a+k,abs(m)-k)
     &                    /factorial(abs(m)-k) )
                  endif

                  if (l_i) then
                     t_ = w_ * beta1 * kappa / ii
                     o0i(ii) = o0i(ii) + t_
                  endif

                  if (l_j) then
                     u_ = -w_ * (mu/2) / jj
                     o0i(jj) = o0i(jj) + u_
                  endif

                  if (l_k) then
                     v_ = w_ / (4*kappa) / kk
                     o0i(kk) = o0i(kk) + v_
                  endif
               enddo
            enddo
         enddo
      endif      


c     ==== Infinite k part ====
      ubound = pmax - abs(m) - 1
      do i = 0, ubound
         do j = 0, ubound
            do k = 0, ubound
               s1 = i + j + k
               ii = abs(m) + s1 + 1
               jj = abs(m) + s1 + 2
               kk = abs(m) + s1 + 3
               l_i = (ii >= abs(m)+1 .and. ii <= pmax)
               l_j = (jj >= abs(m)+2 .and. jj <= pmax)
               l_k = (kk >= abs(m)+3 .and. kk <= pmax)

               if (l_i .or. l_j .or. l_k) then
                  w = v2 * (-1)**i/factorial(i)
     &                *( rise_fact(a,j)/rise_fact(b,j)
     &                /factorial(j) )
     &                *( rise_fact(a,k)/rise_fact(b,k)
     &                /factorial(k) )
                  call wfat_psif(k+a, ps1)
                  call wfat_psif(dble(k+1), ps2)
                  call wfat_psif(dble(k+1+abs(m)), ps3)
                  dk = ps1 - ps2 - ps3
               endif
               
               if (l_i) then
                  mk = dk/ii - 1.0d0/ii**2
                  t = w * beta1 * kappa * mk
                  o0i(ii) = o0i(ii) + t

                  t = w * beta1 * kappa / ii
                  o0il(ii) = o0il(ii) + t
               endif
               
               if (l_j) then
                  mk = dk/jj - 1.0d0/jj**2
                  u = -w * mu/2 * mk
                  o0i(jj) = o0i(jj) + u

                  u = -w * (mu/2) / jj
                  o0il(jj) = o0il(jj) + u
               endif

               if (l_k) then
                  mk = dk/kk - 1.0d0/kk**2
                  v = w / (4*kappa) * mk
                  o0i(kk) = o0i(kk) + v

                  v = w / (4*kappa) / kk
                  o0il(kk) = o0il(kk) + v
               endif
            enddo
         enddo
      enddo
            

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_r0_i_coeff(n, m, erg, zq, mu, pmax, r0i)
c     Calculate the power series coefficients of
c        A1(\eta) = \frac{1}{\kappa} \int_0^\eta d\eta' (R_\nu^{(0)}(\eta'))^2
c                   (\beta_\nu^{(1)}/\eta' - \mu_z/2 + \eta'/4)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      double precision mu            ! [input]
      integer pmax                   ! [input]
      double precision r0i(1:pmax)

      logical l_i, l_j, l_k
      integer i, j, k_, s1, ii, jj, kk, ubound
      double precision kappa, a, b, beta1, ga, v1, w, t, u, v
      double precision factorial, rise_fact
      external factorial, rise_fact


      if (n < 0) call errquit('Argument n of '//
     &     'wfat_r0_i_coeff must be zero or positive.', 28, UERR)
      

      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, ga)
      v1 = ( ga / kappa / factorial(abs(m)) )**2


c     ==== Infinite k part ====
      r0i = 0.0d0
      ubound = pmax - abs(m) - 1
      do i = 0, ubound
         do j = 0, ubound
            do k_ = 0, ubound
               s1 = i + j + k_
               ii = abs(m) + s1 + 1
               jj = abs(m) + s1 + 2
               kk = abs(m) + s1 + 3
               l_i = (ii >= abs(m)+1 .and. ii <= pmax)
               l_j = (jj >= abs(m)+2 .and. jj <= pmax)
               l_k = (kk >= abs(m)+3 .and. kk <= pmax)

               if (l_i .or. l_j .or. l_k) then
                  w = v1 * (-1)**i/factorial(i)
     &                *( rise_fact(a,j)/rise_fact(b,j)
     &                /factorial(j) )
     &                *( rise_fact(a,k_)/rise_fact(b,k_)
     &                /factorial(k_) )
               endif
               
               if (l_i) then
                  t = w * beta1 * kappa / ii
                  r0i(ii) = r0i(ii) + t
               endif
               
               if (l_j) then
                  u = -w * (mu/2) / jj
                  r0i(jj) = r0i(jj) + u
               endif

               if (l_k) then
                  v = w / (4*kappa) / kk
                  r0i(kk) = r0i(kk) + v
               endif
            enddo
         enddo
      enddo      

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_r0o0_i(eta1, eta2, n, m, erg, zq, mu, nlgl,
     &           pt, wt, rr, oo)

c     RR = \frac{1}{\kappa} \int_{\eta_1}^{\eta2} d\eta' (R_\nu^{(0)}(\eta'))^2
c          (\beta_\nu^{(1)}/\eta' - \mu_z/2 + \eta'/4)
c     OO = \frac{1}{\kappa} \int_{\eta_1}^{\eta2} d\eta' O_\nu^{(0)}(\eta')
c          (\beta_\nu^{(1)}/\eta' - \mu_z/2 + \eta'/4) R_\nu^{(0)}(\eta')
      
      implicit none
#include "stdio.fh"
      double precision eta1
      double precision eta2
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      double precision mu            ! [input]
      integer nlgl
      double precision pt(nlgl)
      double precision wt(nlgl)
      double precision rr
      double precision oo

      integer nfem, ngr, id, i, j, md
      double precision pt_4(4), wt_4(4), pt_6(6), wt_6(6)
      double precision kappa, a, b, beta1, g, d, de, dfem, eta1_, eta2_
      double precision eta, wt_, ee, z, r0, o0, ind, min_dfem, factorial
      double precision chgu
      double complex zz, chgm_
      parameter( min_dfem=2.0d0 )
      external factorial
      

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


      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))
      

c     ==== Determine the no. of finite elements and grids ====
c     de = the total interval of r to be integrated
c     dfem = the finite element interval
      de = eta2 - eta1
      if (de < 1.0d-12) then
         nfem = 1
         dfem = de
         ngr = 2
         id = 0
      else if (de >= 1.0d-12 .and. de < 1.0d-3) then
         nfem = 1
         dfem = de
         ngr = 4
         id = 1
      else if (de >= 1.0d-3 .and. de < 1.0d-1) then
         nfem = 1
         dfem = de
         ngr = 6
         id = 2
      else
         nfem = int(max(de,min_dfem+1.0d-13) / min_dfem)
         dfem = de / nfem
         ngr = nlgl
         id = 3
      endif


      rr = 0.0d0
      oo = 0.0d0
c     ==== Loop over finite elements ====
c     eta1_ = lower bound of a FEM
c     eta2_ = upper bound of a FEM
c     eta = shifted quadrature grid
c     wt_ = shifted quadrature weight
      do i = 1, nfem
         eta1_ = eta1 + (i-1)*dfem
         eta2_ = eta1_ + dfem

c        ==== Loop over quadrature grids ====
         do j = 1, ngr
c           ==== Determine the shifted grids and weights ====
            if (id == 0) then
               eta = eta1_ + (j-1)*dfem
               wt_ = (eta2_ - eta1_)/2
            else if (id == 1) then
               eta = 0.5d0*((eta2_-eta1_)*pt_4(j) + eta2_ + eta1_)
               wt_ = (eta2_ - eta1_)/2 * wt_4(j)
            else if (id == 2) then
               eta = 0.5d0*((eta2_-eta1_)*pt_6(j) + eta2_ + eta1_)
               wt_ = (eta2_ - eta1_)/2 * wt_6(j)
            else if (id == 3) then
               eta = 0.5d0*((eta2_-eta1_)*pt(j) + eta2_ + eta1_)
               wt_ = (eta2_ - eta1_)/2 * wt(j)
            endif
            if (eta < 1.0d-13) cycle           ! When eta is zero, the integrands are equal to zero.
            ee = (beta1/eta - mu/2 + eta/4)
            

c           ==== Calculate the necessary Hypergeometric functions ====
            z = kappa * eta
            zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
            call wfat_cchg(a, b, zz, chgm_)
            r0 = d * z**(b/2) * exp(-z/2) * real(chgm_)
            call wfat_chgu(a, b, z, chgu, md)
            o0 = z**(b/2) * exp(-z/2) * chgu
            

c           ==== The R^{(0)} part ====
            ind = r0 * r0 * ee
            rr = rr + wt_*ind

c           ==== The O^{(0)} part ====
            ind = o0 * r0 * ee
            oo = oo + wt_*ind
         enddo
      enddo

      rr = rr / kappa
      oo = oo / kappa
      

      end
c=======================================================
      
c=======================================================
c=======================================================
      subroutine wfat_calc_r0o0_i_smalle(eta, erg, pmax, cri, coi, coil,
     &           ri, oi)

      implicit none
      double precision eta            ! [input],
      double precision erg            ! [input],
      integer pmax                    ! [input],
      double precision cri(pmax)      ! [input],
      double precision coi(pmax)      ! [input],
      double precision coil(pmax)     ! [input],
      double precision ri
      double precision oi
      
      integer j
      double precision kappa, z, ln_z


      kappa = sqrt(2 * abs(erg))
      z = kappa*eta
      ln_z = log(z)

      ri = 0.0d0
      oi = 0.0d0
      do j = 1, pmax
         ri = ri + cri(j) * z**j
         oi = oi + (ln_z*coil(j) + coi(j)) * z**j
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_o0_i_coeff_nomu(n, m, erg, zq, pmax, o0i, o0il)
c     A big portion of this subroutine was copied from wfat_o0_i_coeff because
c     the underlying algebraic structure is identical
      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      integer pmax                   ! [input]
      double precision o0i(pmax,2)   ! [output]
      double precision o0il(pmax,2)  ! [output]

      logical l_i, l_j, l_k
      integer ubound, i, j, k, ii, jj, kk, s1, s2
      double precision kappa, a, b, bt1_1, bt1_2, v1, v2, ga, gm, w_, t_
      double precision u_, v_, w, ps1, ps2, ps3, dk, mk, t, u, v
      double precision factorial, rise_fact
      external factorial, rise_fact
      
      
      if (pmax < 1)
     &     call errquit('WFAT: Argument pmax of wfat_o0_i_coeff_nomu '//
     &     'cannot be smaller than 1.', 48, UERR)
      if (n < 0) call errquit('Argument n of '//
     &     'wfat_o0_i_coeff_nomu must be zero or positive.', 28, UERR)


      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      bt1_1 = -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      bt1_2 = -1.0d0/2/kappa * (2*n+1+abs(m))
      v1 = 1.0d0 / kappa**2 / factorial(abs(m))
      call wfat_gammaf(a, ga)
      call wfat_gammaf(a-abs(m), gm)
      v2 = (-1)**(abs(m)+1) * ga / (kappa*factorial(abs(m)))**2 / gm

            
c     ==== Finite k part ====
      o0i = 0.0d0
      o0il = 0.0d0
      if (abs(m) > 0) then
         ubound = pmax - 1
         do i = 0, ubound
            do j = 0, ubound
               do k = abs(m), 1, -1
                  s2 = i + j - k
                  ii = abs(m) + s2 + 1
                  jj = abs(m) + s2 + 2
                  kk = abs(m) + s2 + 3
                  l_i = (ii >= 1 .and. ii <= pmax)
                  l_j = (jj >= 2 .and. jj <= pmax)
                  l_k = (kk >= 3 .and. kk <= pmax)

                  if (l_i .or. l_j .or. l_k) then
                     w_ = v1 * (-1)**i/factorial(i)
     &                    *( rise_fact(a,j)/factorial(j)
     &                    /rise_fact(b,j) )
     &                    *( factorial(k-1)*rise_fact(1-a+k,abs(m)-k)
     &                    /factorial(abs(m)-k) )
                  endif

                  if (l_i) then
                     t_ = w_ * bt1_1 * kappa / ii
                     o0i(ii,1) = o0i(ii,1) + t_
                     t_ = w_ * bt1_2 * kappa / ii
                     o0i(ii,2) = o0i(ii,2) + t_
                  endif

                  if (l_j) then
                     u_ = -w_ * (1.0d0/2) / jj
                     o0i(jj,2) = o0i(jj,2) + u_
                  endif

                  if (l_k) then
                     v_ = w_ / (4*kappa) / kk
                     o0i(kk,1) = o0i(kk,1) + v_
                  endif
               enddo
            enddo
         enddo
      endif      


c     ==== Infinite k part ====
      ubound = pmax - abs(m) - 1
      do i = 0, ubound
         do j = 0, ubound
            do k = 0, ubound
               s1 = i + j + k
               ii = abs(m) + s1 + 1
               jj = abs(m) + s1 + 2
               kk = abs(m) + s1 + 3
               l_i = (ii >= abs(m)+1 .and. ii <= pmax)
               l_j = (jj >= abs(m)+2 .and. jj <= pmax)
               l_k = (kk >= abs(m)+3 .and. kk <= pmax)

               if (l_i .or. l_j .or. l_k) then
                  w = v2 * (-1)**i/factorial(i)
     &                *( rise_fact(a,j)/rise_fact(b,j)
     &                /factorial(j) )
     &                *( rise_fact(a,k)/rise_fact(b,k)
     &                /factorial(k) )
                  call wfat_psif(k+a, ps1)
                  call wfat_psif(dble(k+1), ps2)
                  call wfat_psif(dble(k+1+abs(m)), ps3)
                  dk = ps1 - ps2 - ps3
               endif
               
               if (l_i) then
                  mk = dk/ii - 1.0d0/ii**2
                  t = w * bt1_1 * kappa * mk
                  o0i(ii,1) = o0i(ii,1) + t
                  t = w * bt1_2 * kappa * mk
                  o0i(ii,2) = o0i(ii,2) + t

                  t = w * bt1_1 * kappa / ii
                  o0il(ii,1) = o0il(ii,1) + t
                  t = w * bt1_2 * kappa / ii
                  o0il(ii,2) = o0il(ii,2) + t
               endif
               
               if (l_j) then
                  mk = dk/jj - 1.0d0/jj**2
                  u = -w * 1.0d0/2 * mk
                  o0i(jj,2) = o0i(jj,2) + u

                  u = -w * (1.0d0/2) / jj
                  o0il(jj,2) = o0il(jj,2) + u
               endif

               if (l_k) then
                  mk = dk/kk - 1.0d0/kk**2
                  v = w / (4*kappa) * mk
                  o0i(kk,1) = o0i(kk,1) + v

                  v = w / (4*kappa) / kk
                  o0il(kk,1) = o0il(kk,1) + v
               endif
            enddo
         enddo
      enddo      

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_r0o0(eta, n, m, erg, zq, ro)

      implicit none
#include "errquit.fh"
      double precision eta
      integer n
      integer m
      double precision erg
      double precision zq
      double precision ro

      integer i, j, imax, jmax, pmax, md
      double precision kp, a, b, g, d, z, chgu, rnu, onu, x
      double precision factorial, rise_fact
      double complex zz, chgm
      parameter (imax=20, jmax=20, pmax=40)
      external factorial, rise_fact


      kp = sqrt(2*abs(erg))
      z = kp * eta
      if (z <= 0.0d0)
     &     call errquit('WFAT: The argument eta of wfat_calc_r0o0 '//
     &     'must be such that sqrt(2*|erg|)*eta is positive', 36, UERR)
      
      
      a = 1 + abs(m) + n - zq/kp
      b = dble(abs(m) + 1)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))

      
      if (z < 30.0d0) then
         zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
         call wfat_cchg(a, b, zz, chgm)
         rnu = d * z**(b/2) * exp(-z/2) * real(chgm)
         call wfat_chgu(a, b, z, chgu, md)
         onu = z**(b/2) * exp(-z/2) * chgu

         ro = rnu * onu
      else

         ro = 0.0d0
         do i = 0, imax
            do j = 0, jmax
               if (i+j <= pmax) then
                  x = (rise_fact(b-a,i) * rise_fact(1-a,i) /
     &                factorial(i)) *
     &                (rise_fact(a,j) * rise_fact(1+a-b,j) /
     &                factorial(j)) * (-1)**j * 
     &                z**(-i-j)
                  ro = ro + x
               endif
            enddo
         enddo

      endif

      end
c=======================================================
      
      
