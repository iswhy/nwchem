c==========================================================
c==========================================================
      subroutine wfat_2a1e_solve(m, nr, nt, geom, ke, mult, adapt_ke, 
     &     maxiter, id_target, pottype, nel, potparam, g_cf, erg)

      implicit none
#include "global.fh"
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      integer nt
      integer geom
      double precision ke                    ! [in/output]
      double precision mult
      logical adapt_ke
      integer maxiter                        ! [input], only used when adapt_ke=.true.
      integer id_target                      ! [input], only used when adapt_ke=.true.
      character*3 pottype
      integer nel                    ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)
      integer g_cf                           ! [input], GA of size nr*nt x nr*nt
      double precision erg(nr*nt)            ! [output],
      
      integer l_lr, l_ur, l_hr, l_gr, l_lt, l_ut, l_ht, l_gt
      integer l_x, l_t, l_rw, l_tw
      integer k_lr, k_ur, k_hr, k_gr, k_lt, k_ut, k_ht, k_gt
      integer k_x, k_t, k_rw, k_tw
      integer ndid, i, j, iter, i_hr, i_x, maxiter_, g_h
      double precision e0, deg2rad
      parameter(deg2rad=0.01745329251994d0)

      
      ndid = ga_nodeid()

c     =====================
c     ==== Radial part ====
c     =====================
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'lapr', l_lr, k_lr))
     &     call errquit('WFAT: Cannot allocate lapr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'u_r', l_ur, k_ur))
     &     call errquit('WFAT: Cannot allocate u_r', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'g_r', l_gr, k_gr))
     &     call errquit('WFAT: Cannot allocate g_r', 1, MA_ERR)
      call wfat_laplr_fbrr(nr, abs(m), dbl_mb(k_lr))
      call wfat_dvr_fbr_r(abs(m), nr, dbl_mb(k_ur))

c     transpose(U_r) . T_r . U_r
      call wfat_represent_sym(nr, nr, dbl_mb(k_ur), dbl_mb(k_lr),
     &     dbl_mb(k_gr))
      
c     ==== Get the radial quadrature nodes ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'h_r', l_hr, k_hr))
     &     call errquit('WFAT: Cannot allocate h_r', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'xweight', l_rw, k_rw))
     &     call errquit('WFAT: Cannot allocate xweight', 1, MA_ERR)
      call wfat_lf_quadrature_rule(nr, abs(dble(m))+2, dbl_mb(k_x),
     &     dbl_mb(k_rw))
      do j = 1, nr
         do i = 1, nr
            i_hr = k_hr + (j-1)*nr + (i-1)
            if (i == j) then
               i_x = k_x + i - 1
               dbl_mb(i_hr) = 1/(dbl_mb(i_x)**2)
            else
               dbl_mb(i_hr) = 0.0d0
            endif
         enddo
      enddo
c     =====================

      
c     ====================
c     ==== Theta part ====
c     ====================
      if (.not. ma_alloc_get(MT_DBL, nt, 'lapt', l_lt, k_lt))
     &     call errquit('WFAT: Cannot allocate lapt', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'u_t', l_ut, k_ut))
     &     call errquit('WFAT: Cannot allocate u_t', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'i_t', l_ht, k_ht))
     &     call errquit('WFAT: Cannot allocate i_t', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'g_t', l_gt, k_gt))
     &     call errquit('WFAT: Cannot allocate g_t', 1, MA_ERR)

      do j = 1, nt
         dbl_mb(k_lt+j-1) = -dble((abs(m)+j-1) * (abs(m)+j))
      enddo
      call wfat_dvr_fbr_t(abs(m), nt, dbl_mb(k_ut))
c     transpose(U_t) . T_t . U_t
      call wfat_represent_diag(nt, nt, dbl_mb(k_ut), dbl_mb(k_lt),
     &     dbl_mb(k_gt))
      call wfat_get_identity(nt, dbl_mb(k_ht))

c     ==== Get the polar quadrature nodes ====
      if (.not. ma_alloc_get(MT_DBL, nt, 'tquad', l_t, k_t))
     &     call errquit('WFAT: Cannot allocate tquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt, 'tweight', l_tw, k_tw))
     &     call errquit('WFAT: Cannot allocate tweight', 1, MA_ERR)
      call wfat_assocleg_quad(-1.0d0, 1.0d0, abs(m), nt, dbl_mb(k_t),
     &     dbl_mb(k_tw))
      dbl_mb(k_t:k_t+nt-1) = acos(dbl_mb(k_t:k_t+nt-1)) / deg2rad
c     ====================

      
c     ==== Hamiltonian construction ====
      if (adapt_ke) then
         maxiter_ = maxiter
      else
         maxiter_ = 1
      endif
      
      e0 = 0.0d0
      do iter = 1, maxiter_
         if (.not. ga_create(mt_dbl, nr*nt, nr*nt, 'hamil', nr*nt, 0,
     &        g_h)) call errquit('WFAT: Cannot create hamil GA.',
     &        1, GA_ERR)
         call wfat_2a1e_hamiltonian(nr, nt, dbl_mb(k_x), dbl_mb(k_t),
     &        dbl_mb(k_gr), dbl_mb(k_ht), dbl_mb(k_hr), dbl_mb(k_gt), 
     &        ke, geom, pottype, nel, potparam, g_h)
      
c        ==== Hamiltonian diagonalization ====
         call ga_diag_std(g_h, g_cf, erg)
         call ga_sync()

         if (adapt_ke) then
            if (ndid == 0)
     &           write(luout, '(i5, ") energy, delta energy = ", '//
     &           '2f15.8)') iter, erg(id_target), erg(id_target) - e0
            e0 = erg(id_target)
c            ke = sqrt( 2 * abs(erg(id_target)) ) * 2
            ke = sqrt( 2 * abs(erg(id_target)) ) * mult
c            ke = sqrt( 2 * abs(erg(id_target)) ) * 9
         endif
      enddo
      
      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_lr))
     &     call errquit('WFAT: Cannot deallocate lapr', 1, MA_ERR)
      if (.not. ma_free_heap(l_ur))
     &     call errquit('WFAT: Cannot deallocate u_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_hr))
     &     call errquit('WFAT: Cannot deallocate i_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_gr))
     &     call errquit('WFAT: Cannot deallocate g_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_x))
     &     call errquit('WFAT: Cannot deallocate xquad', 1, MA_ERR)
      if (.not. ma_free_heap(l_rw))
     &     call errquit('WFAT: Cannot deallocate xweight', 1, MA_ERR)

      if (.not. ma_free_heap(l_lt))
     &     call errquit('WFAT: Cannot deallocate lapt', 1, MA_ERR)
      if (.not. ma_free_heap(l_ut))
     &     call errquit('WFAT: Cannot deallocate u_t', 1, MA_ERR)
      if (.not. ma_free_heap(l_ht))
     &     call errquit('WFAT: Cannot deallocate i_t', 1, MA_ERR)
      if (.not. ma_free_heap(l_gt))
     &     call errquit('WFAT: Cannot deallocate g_t', 1, MA_ERR)
      if (.not. ma_free_heap(l_t))
     &     call errquit('WFAT: Cannot deallocate tquad', 1, MA_ERR)
      if (.not. ma_free_heap(l_tw))
     &     call errquit('WFAT: Cannot deallocate tweight', 1, MA_ERR)

      if (.not. ga_destroy(g_h))
     &     call errquit('WFAT: Cannot destroy hamil.', 40, GA_ERR)
      
      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_2a1e_hamiltonian(nr, nt, xquad, tquad, gr, ht, hr,
     &     gt, ke, geom, pottype, nel, potparam, g_h)

      implicit none
#include "mafdecls.fh"      
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer nr
      integer nt
      double precision xquad(nr)     ! [input], scaled radial quadrature nodes.
      double precision tquad(nt)     ! [input], theta quadrature nodes in degrees.
      double precision gr(nr, nr)
      double precision ht(nt, nt)
      double precision hr(nr, nr)
      double precision gt(nt, nt)
      double precision ke
      integer geom
      character*3 pottype
      integer nel                    ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)      ! 1)
      integer g_h
      
c     1) When pottype = 'cou', potparam(1) = soft Coulomb parameter, while potparam(2) is unused.
c        When pottype = 'eff', potparam(1:2) = the two screening parameters.

      integer i, g_tm1, g_tm2, g_pot, i_r, i_t, lo, hi, ndid
      double precision x, t, vpt(3), pot
      double precision deg2rad
      parameter(deg2rad=0.01745329251994d0)

      
      call ga_sync()
      ndid = ga_nodeid()

      if (.not. ga_create(mt_dbl, nr*nt, nr*nt, 'tmp1', nr*nt, 0,
     &     g_tm1)) call errquit('WFAT: Cannot create tmp1 global array.'
     &     , 1, GA_ERR)
      if (.not. ga_create(mt_dbl, nr*nt, nr*nt, 'tmp2', nr*nt, 0,
     &     g_tm2)) call errquit('WFAT: Cannot create tmp2 global array.'
     &     , 1, GA_ERR)
      call wfat_ga_kron_prod(nr, nr, nt, nt, gr, ht, g_tm1)
      call wfat_ga_kron_prod(nr, nr, nt, nt, hr, gt, g_tm2)

c     ==== Kinetic part ====
      call ga_add(1.0d0, g_tm1, 1.0d0, g_tm2, g_h)
      call ga_scale(g_h, -0.5d0*ke**2)
c     g_h now contains kinetic energy matrix.
      
c     ==== Nuclear attraction part ====
      if (.not. nga_create(mt_dbl, 1, nr*nt, 'pot', 0, g_pot)) 
     &     call errquit('WFAT: Cannot create tmp2 global array.', 
     &     1, GA_ERR)           ! 1)
c     1) For some reason, ga_add_diagonal appearing later requires strictly
c        one dimensional array (vector), and the 1-D array created through
c        ga_create (instead of nga_create) turns out to be unconforming to
c        this requirement.
      call ga_zero(g_pot)
      call nga_distribution(g_pot, ga_nodeid(), lo, hi)

      do i = lo, hi
         call wfat_invert_kron_id(i, nt, i_r, i_t)
         x = xquad(i_r)
         t = tquad(i_t) * deg2rad

         vpt = (/ x/ke*sin(t), 0.0d0, x/ke*cos(t) /)
         call wfat_2a1e_potential(geom, vpt, pottype, nel, potparam, 
     &        pot)
         call nga_put(g_pot, i, i, pot, 1)
      enddo
      call ga_sync()

c     ==== Total Hamiltonian ====
      call ga_add_diagonal(g_h, g_pot)
c     Only add potential to the diagonal elements since in DVR basis,
c     potential matrix is diagonal.
                  
c     ==== Deallocations ====
      if (.not. ga_destroy(g_tm1))
     &     call errquit('WFAT: Cannot destroy tmp1.', 1, GA_ERR)
      if (.not. ga_destroy(g_tm2))
     &     call errquit('WFAT: Cannot destroy tmp2.', 1, GA_ERR)
      if (.not. ga_destroy(g_pot))
     &     call errquit('WFAT: Cannot destroy pot.', 1, GA_ERR)

      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_laplr_fbrr(n, m, lapr)
c     DESCRIPTION:
c     Calculates the matrix element of the radial derivative term of the Laplacian
c     operator in the radial FBR basis based the generalized Laugerre polynomials.
c     $$
c        \Delta^r_{ij} = \int_0^\infty dr r^2 R_i^{|m|}(r) T(R_j^{|m|}(r))
c     $$
c     where
c     \begin{align}
c        T(f(r)) &= \frac{1}{r} \frac{d^2}{dr^2} (rf(r)) ,\\
c        R_i^{|m|}(r) &= \frac{k^{3/2}}{ \sqrt{h_{i-1,|m|+2}} } (kr)^{|m|/2} \exp(-kr/2)  \\
c                        L_{i-1}^{|m|+2}(kr) ,\\
c        h_{i,m} &= \Gamma(|m|+1) \binom{|m|+i}{i}.
c     \end{align}
c     for $k=1$. The matrix for $k \neq 1$ (but $k>0$) may be obtained by multiplying the 
c     output of this subroutine (lapr) with $k^2$.
      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n                    ! [input], the size of the FBR bases, n >= 1.
      integer m                    ! [input], the generalization parameter of the Laguerre polynomials, m >= 0.
      double precision lapr(n,n)   ! [output], $\Delta^r$, the matrix representation of radial derivative part of the
                                   !           Laplacian operator in the FBR basis.

      integer i, j
      double precision mm, v1, v2, v3, v4, v5, v6, normi, normj
      double precision laguerre_norm_i
      external laguerre_norm_i


      if (n < 1) call errquit('WFAT: The input argument n of '//
     &     'wfat_laplr_fbrr must be such that n >= 1', 58, UERR)
      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_laplr_fbrr must be such that m >= 0', 58, UERR)

      mm = dble(m)
      do j = 1, n
         normj = sqrt(laguerre_norm_i(j-1,m+2))
         do i = 1, n
c         do j = 1, i

            normi = sqrt(laguerre_norm_i(i-1,m+2))
            
            call wfat_laguerre_intg1(mm+2, mm+2, mm+0, i-1, j-1, v1)
            call wfat_laguerre_intg1(mm+2, mm+2, mm+1, i-1, j-1, v2)
            call wfat_laguerre_intg1(mm+2, mm+2, mm+2, i-1, j-1, v3)
            if (j-2 < 0) then
               v4 = 0.0d0
               v5 = 0.0d0
            else
               call wfat_laguerre_intg1(mm+2, mm+3, mm+1, i-1, j-2, v4)
               call wfat_laguerre_intg1(mm+2, mm+3, mm+2, i-1, j-2, v5)
            endif
            if (j-3 < 0) then
               v6 = 0.0d0
            else
               call wfat_laguerre_intg1(mm+2, mm+4, mm+2, i-1, j-3, v6)
            endif

            lapr(i,j) = 1/normi/normj *
     &           ((mm/2+1)*mm/2*v1 - (mm/2+1)*v2 + 0.25d0*v3 -
     &           2*(mm/2+1)*v4 + v5 + v6)

c            lapr(j,i) = lapr(i,j)
         enddo
      enddo

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_eval_fbrr(m, nr, ke, r, rfun)
c     DESCRIPTION:
c     Evaluates the radial FBR basis based on the associated Laguerre functions.
      implicit none
#include "errquit.fh"
      integer m
      integer nr
      double precision ke
      double precision r
      double precision rfun(nr)

      integer i
      double precision mm, x, laguerre_norm_i
      external laguerre_norm_i

      x = ke * r
      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_eval_fbrr must be such that m >= 0', 58, UERR)
      if (x < 0.0d0)
     &     call errquit('WFAT: The input argument r of '//
     &     'wfat_eval_fbrr must be such that r >= 0.', 80, UERR)
      
      mm = dble(m)
      call wfat_lagfun(1, nr-1, mm+2, x, rfun)
c     At this point, rfun(1:nr) contains the generalized Laguerre polynomials
c     L_n^{|m|+2} from n=0 to n=nr-1 evaluated at x.

      do i = 1, nr
         rfun(i) = sqrt(ke**3) * 1/sqrt(laguerre_norm_i(i-1, m+2)) * 
     &        x**(mm/2) * exp(-x/2) * rfun(i)
      enddo




cc      mm = dble(m)
cc      call wfat_lagfun_e(1, nr-1, mm+2, x, exp(-x/2), rfun)
ccc     At this point, rfun(1:nr) contains the generalized Laguerre polynomials
ccc     L_n^{|m|+2} from n=0 to n=nr-1 evaluated at x.
cc
cc      do i = 1, nr
cc         rfun(i) = sqrt(ke**3) * 1/sqrt(laguerre_norm_i(i-1, m+2)) * 
cc     &        x**(mm/2) * rfun(i)
cc      enddo

      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_eval_fbrt(m, nt, theta, tfun)
c     DESCRIPTION:
c     Evaluates the theta FBR basis based on the associated Legendre functions.
      implicit none
#include "errquit.fh"
      integer m                     ![input], the magnetic quantum number, m >= 0.
      integer nt                    ![input], the number of theta FBR basis states.
      double precision theta        ![input], polar angle in degrees to evaluate tfun.
      double precision tfun(1:nt)   ![output], the theta FBR basis states.

      integer i, l
      double precision x, th, assocleg_norm, deg2rad
      parameter(deg2rad=0.01745329251994d0)
      external assocleg_norm


      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_eval_fbrt must be such that m >= 0', 58, UERR)
      if (theta < 0.0d0 .or. theta > 180.0d0)
     &     call errquit('WFAT: The input argument theta of '//
     &     'wfat_eval_fbrt must be such that 0 <= theta <= 180.',
     &     80, UERR)

      th = theta * deg2rad
      x = cos(th)
      call wfat_assocleg_fixm(m, m+nt-1, x, tfun)

      do i = 1, nt
         l = m + i - 1
         tfun(i) = tfun(i) / sqrt(assocleg_norm(l, m))
      enddo

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_eval_dvrr(m, nr, ke, r, dvrr)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      double precision ke
      double precision r
      double precision dvrr(nr)

      integer i, l_fb, k_fb, l_u, k_u, i_u


      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_eval_dvrr must be such that m >= 0', 58, UERR)

      if (.not. ma_alloc_get(MT_DBL, nr, 'fbrr', l_fb, k_fb))
     &     call errquit('WFAT: Cannot allocate fbrr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'u_r', l_u, k_u))
     &     call errquit('WFAT: Cannot allocate u_r', 1, MA_ERR)


      call wfat_dvr_fbr_r(m, nr, dbl_mb(k_u))
      call wfat_eval_fbrr(m, nr, ke, r, dbl_mb(k_fb))

      do i = 1, nr
         i_u = k_u + (i-1)*nr
         dvrr(i) = sum(dbl_mb(i_u:i_u+nr-1) * dbl_mb(k_fb:k_fb+nr-1))
      enddo


      if (.not. ma_free_heap(l_fb))
     &     call errquit('WFAT: Cannot deallocate fbrr', 1, MA_ERR)
      if (.not. ma_free_heap(l_u))
     &     call errquit('WFAT: Cannot deallocate u_r', 1, MA_ERR)

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_eval_dvrt(m, nt, theta, dvrt)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m                  ![input], the magnetic quantum number, m >= 0.
      integer nt                 ![input], the number of theta DVR basis states.
      double precision theta     ![input], polar angle in degrees to evaluate dvrt.
      double precision dvrt(nt)  ![output], the theta DVR basis states.

      integer i, l_fb, k_fb, l_u, k_u, i_u


      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_eval_dvrt must be such that m >= 0', 58, UERR)

      if (.not. ma_alloc_get(MT_DBL, nt, 'fbrt', l_fb, k_fb))
     &     call errquit('WFAT: Cannot allocate fbrt', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'u_t', l_u, k_u))
     &     call errquit('WFAT: Cannot allocate u_t', 1, MA_ERR)
      

      call wfat_dvr_fbr_t(m, nt, dbl_mb(k_u))
      call wfat_eval_fbrt(m, nt, theta, dbl_mb(k_fb))

      do i = 1, nt
         i_u = k_u + (i-1)*nt
         dvrt(i) = sum(dbl_mb(i_u:i_u+nt-1) * dbl_mb(k_fb:k_fb+nt-1))
      enddo

      
      if (.not. ma_free_heap(l_fb))
     &     call errquit('WFAT: Cannot deallocate fbrt', 1, MA_ERR)
      if (.not. ma_free_heap(l_u))
     &     call errquit('WFAT: Cannot deallocate u_t', 1, MA_ERR)

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_2a1e_eval_phi(s, m, phi, f)

      implicit none
      character*1 s
      integer m
      double precision phi
      double precision f

      double precision pi_const
      parameter(pi_const=3.14159265359d0)

      if (s == 'c') then
         if (m == 0) then
            f = 1.0d0 / sqrt(2*pi_const)
         else
            f = cos(abs(m)*phi) / sqrt(pi_const)
         endif
      else if (s == 's') then
         f = sin(abs(m)*phi) / sqrt(pi_const)
      endif

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_dvr_fbr_r(m, nr, u)
c     DESCRIPTION:
c     Calculates the transformation matrix from radial DVR basis to radial
c     FBR basis. The (i,j)-th element of this transformation matrix is given by
c     $$
c     u_{ij} = \sqrt{w_j} L_{i-1}^{|m|+2}(x_j) / \sqrt{h_{i-1,|m|+2}}
c     $$
c     where $x_i$ and $w_i$ is the i-th zero of $L_{nr}^{|m|+2}$ and the associated
c     quadrature weight in the Gauss-Laguerre quadrature rule. $h_{i-1,|m|+2}$ is the
c     normalization factor for $L_{i-1}^{|m|+2}$.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      double precision u(nr,nr)

      integer l_pt, k_pt, l_wt, k_wt, l_fb, k_fb, j, i
      double precision mm, w, x, lag, laguerre_norm_i
      external laguerre_norm_i

      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_dvr_fbr_r must be such that m >= 0', 58, UERR)

      mm = dble(m)
      if (.not. MA_alloc_get(MT_DBL, nr, 'pt', l_pt, k_pt))
     &     call errquit('WFAT: pt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nr, 'wt', l_wt, k_wt))
     &     call errquit('WFAT: wt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nr, 'fbr', l_fb, k_fb))
     &     call errquit('WFAT: fbr allocation failed', 59, MA_ERR)


c     ==== Get the quadrature nodes and weights ====
      call wfat_lf_quadrature_rule(nr, mm+2, dbl_mb(k_pt),
     &     dbl_mb(k_wt))

      
c     ==== Compute the transformation matrix elements ====
      do j = 1, nr      ! Loop over DVR basis indices (columns)
         w = dbl_mb(k_wt+j-1)
         x = dbl_mb(k_pt+j-1)
         call wfat_lagfun(1, nr-1, mm+2, x, dbl_mb(k_fb))
         do i = 1, nr      ! Loop over FBR basis indices (rows)
            lag = dbl_mb(k_fb+i-1)
            u(i,j) = sqrt(w) * lag / sqrt(laguerre_norm_i(i-1, m+2))
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate pt', 32, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate wt', 32, MA_ERR)
      if (.not. ma_free_heap(l_fb))
     &     call errquit('WFAT: Cannot deallocate fb', 32, MA_ERR)
      
      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_dvr_fbr_t(m, nt, u)
c     DESCRIPTION:
c     Calculates the transformation matrix from theta DVR basis to theta
c     FBR basis.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nt
      double precision u(nt,nt)

      integer j, i, l_pt, k_pt, l_wt, k_wt, l_fb, k_fb
      double precision mm, w, x, fbrt, deg2rad
      parameter (deg2rad=0.01745329251994d0)

      
      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_dvr_fbr_t must be such that m >= 0', 58, UERR)
      mm = dble(m)

      if (.not. MA_alloc_get(MT_DBL, nt, 'pt', l_pt, k_pt))
     &     call errquit('WFAT: pt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nt, 'wt', l_wt, k_wt))
     &     call errquit('WFAT: wt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nt, 'fbr', l_fb, k_fb))
     &     call errquit('WFAT: fbr allocation failed', 59, MA_ERR)

      
c     ==== Get the quadrature nodes and weights ====
      call wfat_assocleg_quad(-1.0d0, 1.0d0, m, nt, dbl_mb(k_pt),
     &     dbl_mb(k_wt))


c     ==== Compute the transformation matrix elements ====
      do j = 1, nt              ! Loop over FBR bases.
         w = dbl_mb(k_wt+j-1)
         x = dbl_mb(k_pt+j-1)
c         write(luout, '("X = ", es16.8)') x
         call wfat_eval_fbrt(m, nt, acos(x)/deg2rad, dbl_mb(k_fb))
         do i = 1, nt           ! Loop over DVR bases.
            fbrt = dbl_mb(k_fb+i-1)
            u(i,j) = sqrt(w) * (1-x**2)**(-mm/2) * fbrt

c            write(luout, '(2i5, 3es16.8)')
c     &           i, j, fbrt, w, (1-x**2)**(-mm/2)
         enddo
      enddo

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate pt', 32, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate wt', 32, MA_ERR)
      if (.not. ma_free_heap(l_fb))
     &     call errquit('WFAT: Cannot deallocate fb', 32, MA_ERR)
      
      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_eval_fbrr_mask(m, nr, ke, r, a, b, zi, kp, rfun)
c     DESCRIPTION:
c     Evaluates the radial FBR basis based on the associated Laguerre functions.
      implicit none
#include "errquit.fh"
      integer m
      integer nr
      double precision ke
      double precision r
      double precision a
      double precision b
      double precision zi
      double precision kp
      double precision rfun(nr)

      integer i
      double precision mm, x, laguerre_norm_i, er, p, d
      external laguerre_norm_i

      x = ke * r
      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_eval_fbrr must be such that m >= 0', 58, UERR)
      if (x < 0.0d0)
     &     call errquit('WFAT: The input argument r of '//
     &     'wfat_eval_fbrr must be such that r >= 0.', 80, UERR)
      
      mm = dble(m)
      call wfat_lagfun_e(1, nr-1, mm+2, x, exp(-x/2), rfun)
c     At this point, rfun(1:nr) contains the generalized Laguerre polynomials
c     L_n^{|m|+2} from n=0 to n=nr-1 evaluated at x.

      er = erfc(a*(x-b))
      do i = 1, nr
         ! WARNING: Below, one might need to divide by x**i instead of x**(i-1)
         ! in order for the radial basis to behave like hydrogen eigenfunctions
         ! at large distances, which is r**(n-1) exp(-r/n).
         p = 1 / x**(i-1) / x**(mm/2) * x**(zi/kp)
         d = p - (p-1) / 2 * er
         rfun(i) = sqrt(ke**3) * 1/sqrt(laguerre_norm_i(i-1, m+2)) * 
     &        x**(mm/2) * rfun(i) * d
      enddo

      end
c     NOTES:
c     1) The division by x**(i-1) is to neutralize the highest polynomial term
c        in the generalized Laguerre polynomial. While the division by x**(mm/2)
c        is to cancel the same prefactor in the numerator. As a result, the
c        basis goes to x**(zi/kp) exp(-x/2) as x goes to infinity.
c==========================================================
      

c==========================================================
c==========================================================
      subroutine wfat_eval_dvrr_mask(m, nr, ke, r, a, b, zi, kp, dvrr)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      double precision ke
      double precision r
      double precision a
      double precision b
      double precision zi
      double precision kp
      double precision dvrr(nr)

      integer i, l_fb, k_fb, l_u, k_u, i_u


      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_eval_dvrr must be such that m >= 0', 58, UERR)

      if (.not. ma_alloc_get(MT_DBL, nr, 'fbrr', l_fb, k_fb))
     &     call errquit('WFAT: Cannot allocate fbrr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'u_r', l_u, k_u))
     &     call errquit('WFAT: Cannot allocate u_r', 1, MA_ERR)


      call wfat_dvr_fbr_r(m, nr, dbl_mb(k_u))
      call wfat_eval_fbrr_mask(m, nr, ke, r, a, b, zi, kp, dbl_mb(k_fb))

      do i = 1, nr
         i_u = k_u + (i-1)*nr
         dvrr(i) = sum(dbl_mb(i_u:i_u+nr-1) * dbl_mb(k_fb:k_fb+nr-1))
      enddo


      if (.not. ma_free_heap(l_fb))
     &     call errquit('WFAT: Cannot deallocate fbrr', 1, MA_ERR)
      if (.not. ma_free_heap(l_u))
     &     call errquit('WFAT: Cannot deallocate u_r', 1, MA_ERR)

      end
c==========================================================
      
