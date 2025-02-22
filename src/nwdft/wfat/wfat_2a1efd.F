c==========================================================
c==========================================================
      subroutine wfat_2a1efd_solve(m, ni, n, maxn, xq, wq, nt, geom,
     &     pottype, nel, potparam, g_cf, erg)

      implicit none
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer ni
      integer n(ni)
      integer maxn
      double precision xq(maxn,ni)
      double precision wq(maxn,ni)
      integer nt
      integer geom
      character*3 pottype
      integer nel                    ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)
      integer g_cf               ! [output], GA of size nr*nt x nr*nt
      double precision erg(*)    ! [output], size nr*nt
      
      integer l_hr, l_gr, l_lt, l_ut, l_ht, l_gt
      integer l_x, l_t, l_tw
      integer k_hr, k_gr, k_lt, k_ut, k_ht, k_gt
      integer k_x, k_t, k_tw
      integer nr, i, j, i_hr, i_x, ndid, g_h
      double precision tt1, tt2, deg2rad
      parameter(deg2rad=0.01745329251994d0)


      ndid = ga_nodeid()

c     =====================
c     ==== Radial part ====
c     =====================
      call wfat_get_fedvrnr(ni, n, nr)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'g_r', l_gr, k_gr))
     &     call errquit('WFAT: Cannot allocate g_r', 3, MA_ERR)
      call wfat_2a1efd_d2r(ni, n, maxn, xq, wq, nr, dbl_mb(k_gr))
      
c     ==== Get the radial quadrature nodes ====
      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad_', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad_', 3, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'h_r', l_hr, k_hr))
     &     call errquit('WFAT: Cannot allocate h_r', 3, MA_ERR)
      call wfat_fedvr_basnodes(ni, n, maxn, xq, dbl_mb(k_x))
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
      if (.not. ga_duplicate(g_cf, g_h, 'hamil'))
     &     call errquit('WFAT: Cannot create hamil.', 36, GA_ERR)
c      if (.not. ga_create(mt_dbl, nr*nt, nr*nt, 'hamil', nr*nt, 0, g_h))
c     &     call errquit('WFAT: Cannot create hamil global array.',
c     &     1, GA_ERR)
      call wfat_2a1e_hamiltonian(nr, nt, dbl_mb(k_x), dbl_mb(k_t),
     &     dbl_mb(k_gr), dbl_mb(k_ht), dbl_mb(k_hr), dbl_mb(k_gt), 
     &     1.0d0, geom, pottype, nel, potparam, g_h)
      
c     ==== Hamiltonian diagonalization ====
      call ga_diag_std(g_h, g_cf, erg)
      call ga_sync()
      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_gr))
     &     call errquit('WFAT: Cannot deallocate g_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_x))
     &     call errquit('WFAT: Cannot deallocate xquad_', 1, MA_ERR)
      if (.not. ma_free_heap(l_hr))
     &     call errquit('WFAT: Cannot deallocate i_r', 1, MA_ERR)
      
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
     &     call errquit('WFAT: Cannot destroy hamil.', 45, GA_ERR)
      
      end
c==========================================================
      

c=======================================================
c=======================================================
      subroutine wfat_2a1efd_d2r(ni, n, maxn, xq, wq, nb, t)
c     Computes the matrix representation of the $d^2/dr^2$ operator in FEDVR basis.
      implicit none
#include "errquit.fh"
      integer ni
      integer n(ni)
      integer maxn             ! maxn = maxval(n)
      double precision xq(maxn,ni)
      double precision wq(maxn,ni)
      integer nb               ! nb=sum(n)-ni*2+ni-1
      double precision t(nb,nb)

      logical brg1, brg2, include1, include2
      integer nb_, j1, i1, m1, j2, i2, m2


      nb_ = (sum(n) - ni*2) + (ni - 1)   ! 1)
      if (nb /= nb_) call errquit('WFAT: Inconsistency found among '//
     &     'the value of the arguments nb, n, and ni of '//
     &     'wfat_2a1efd_kin.', 95, UERR)
      

c     ==== Loop over ket ====
      j2 = 0
      do i2 = 1, ni
         do m2 = 1, n(i2)
c           == Bridge or regular? ==
            brg2 = .false.
            if (m2 == n(i2)) brg2 = .true.

c           == Exclude (1,1) (the 1st basis), (ni,n(ni)) (the last ==
c           == basis), and (:,1) (have been absorbed into bridges) ==
            include2 = .true.
            if ((i2==1 .and. m2==1) .or. (i2==ni .and. m2==n(ni))
     &           .or. m2==1) include2 = .false.
            if (include2) j2 = j2 + 1

c           ==== Loop over bra ====
            j1 = 0
            do i1 = 1, ni
               do m1 = 1, n(i1)
c                 == Bridge or regular? ==
                  brg1 = .false.
                  if (m1 == n(i1)) brg1 = .true.

c                 == Exclude (1,1) (the 1st basis), (ni,n(ni)) (the last ==
c                 == basis), and (:,1) (have been absorbed into bridges) ==
                  include1 = .true.
                  if ((i1==1 .and. m1==1) .or. (i1==ni .and. m1==n(ni))
     &                 .or. m1==1) include1 = .false.
                  if (include1) j1 = j1 + 1

c                 == Calculate the kinetic energy matrix element ==
                  if (include1 .and. include2) then
                     call wfat_2a1efd_d2r_elm(brg1, i1, m1, brg2, i2, 
     &                    m2, ni, n, maxn, xq, wq, t(j1,j2))
                  endif

               enddo
            enddo

         enddo
      enddo

      if (j1 /= nb)
     &     call errquit('WFAT: j1 /= nb, this is a bug.', 40, UERR)
      if (j2 /= nb)
     &     call errquit('WFAT: j2 /= nb, this is a bug.', 40, UERR)
      
      end
c     NOTES:
c     1) sum(n)-ni*2 is the number of regular FEDVR bases while ni-1 is the number of bridge
c        FEDVR bases.
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_2a1efd_d2r_elm(brg1, i1, m1, brg2, i2, m2, ni, n,
     &     maxn, xq, wq, t)
c     Computes the matrix element of the $d^2/dr^2$ operator in FEDVR basis, that is,
c     it calculates
c     $$
c       \int_0^\infty \chi_{i,m}(r) \frac{d^2}{dr^2} \chi_{i',m'}(r) .
c     $$
c     where $\chi_{i,m}$ is the FEDVR basis that corresponds to the i-th FE and m-th DVR
c     within that FE.
      implicit none
#include "errquit.fh"
      logical brg1
      integer i1
      integer m1        ! 2 <= m1 <= n(i1)-1. Ignored if brg1=.true..
      logical brg2
      integer i2
      integer m2        ! 2 <= m2 <= n(i2)-1, Ignored if brg2=.true..
      integer ni
      integer n(ni)
      integer maxn             ! maxn = maxval(n)
      double precision xq(maxn,ni)
      double precision wq(maxn,ni)
      double precision t

      integer i, m, n1, n1_, n2, n2_
      double precision df1, df2
      

      n1 = n(i1)
      n1_ = n(i1+1)
      n2 = n(i2)
      n2_ = n(i2+1)

      if ((m1 == 1 .or. m1 == n1) .and. .not. brg1)
     &     call errquit('WFAT: The argument m1 of wfat_2a1efd_kin_elm'//
     &     ' must be such that 2 <= m1 <= n(i1)-1.', 17, UERR)
      if ((m2 == 1 .or. m2 == n2) .and. .not. brg2)
     &     call errquit('WFAT: The argument m2 of wfat_2a1efd_kin_elm'//
     &     ' must be such that 2 <= m2 <= n(i2)-1.', 17, UERR)
      

c     ==== Perform quadrature summation ====
      t = 0.0d0
      do i = 1, ni     ! Loop over FEs.
         do m = 1, n(i)       ! Loop over quadrature points in the i-th FE.
            
c           ==== Evaluate Bra at the (i,m)-th point ====
            if (brg1) then
               if (i == i1 .or. i == i1+1) then
                  call wfat_eval1_fedvr_bridge(i1, n1, xq(:,i1), 
     &                 wq(n1,i1), n1_, xq(:,i1+1), wq(1,i1+1),
     &                 i, m, df1)
               else
                  df1 = 0.0d0
               endif
            else
               if (i == i1) then
                  call wfat_eval1_fedvr(m1, n(i1), xq(:,i1), 
     &                 wq(m1,i1), m, df1)
               else
                  df1 = 0.0d0
               endif
            endif

c           ==== Evaluate Ket at the (i,m)-th point ====
            if (brg2) then
               if (i == i2 .or. i == i2+1) then
                  call wfat_eval1_fedvr_bridge(i2, n2, xq(:,i2), 
     &                 wq(n2,i2), n2_, xq(:,i2+1), wq(1,i2+1),
     &                 i, m, df2)
               else
                  df2 = 0.0d0
               endif
            else
               if (i == i2) then
                  call wfat_eval1_fedvr(m2, n(i2), xq(:,i2), 
     &                 wq(m2,i2), m, df2)
               else
                  df2 = 0.0d0
               endif
            endif
            
c           ==== Quadrature sum ====
            t = t + wq(m,i) * df1 * df2
         enddo
      enddo

      t = -t

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_eval1_fedvr(m, n, x, w, mm, f)
c     DESCRIPTION:
c     Evaluates a regular (non-bridge) FEDVR basis at one of the FEDVR quadrature points
c     WITHIN THE SAME FE. The value of the regular FEDVR basis in the other FEs is zero.
c     The FEDVR basis to be evaluated is characterized by m, n, x, and w, while the
c     evaluation point is specified by mm.
      implicit none
#include "errquit.fh"
      integer m             ! [input]. The DVR index of the FEDVR basis in the i-th FE.
      integer n             ! [input]. The number of DVR points in the i-th FE.
      double precision x(*) ! [input]. The DVR quadrature points in the i-th FE. Needed when m /= mm. When needed, its minimum size should size should be n.
      double precision w    ! [input]. The (i,m)-th quadrature weight.
      integer mm            ! [input]. The DVR index of the evaluation point in the i-th FE.
      double precision f    ! [input]. The value of the (i,m)-th FEDVR basis at the (i,mm)-th quadrature point.

      double precision f0

      if (m < 1 .or. m > n)
     &     call errquit('WFAT: The argument m of '//
     &     'wfat_eval1_lobatto_dm must be in between 1 and n.', 1, UERR)
      if (mm < 1 .or. mm > n)
     &     call errquit('WFAT: The argument mm of '//
     &     'wfat_eval1_lobatto_dm must be in between 1 and n.', 1, UERR)

      if (m /= mm) then
         call wfat_eval1_lobatto_dm(m, n, mm, x, f0)
      else
         call wfat_eval1_lobatto_em(m, n, w, f0)
      endif
      f = f0 / sqrt(w)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval1_fedvr_bridge(i, n1, x1, w1, n2, x2, w2, ii,
     &     mm, f)
c     DESCRIPTION:
c     Evaluates a bridge FEDVR basis at one of the FEDVR quadrature points WITHIN THE SAME
C     OR THE NEXT FE. The value of the bridge FEDVR basis outside of these two FEs is zero.
c     The given value of i represents the bridge function formed by the (i,n1)- and (i+1,1)-th
c     Lobatto shape functions, i.e. the bridge function is the one that bridges the i-th and
c     the (i+1)-th FEs. The bridge FEDVR basis to be evaluated is characterized by i,
c     n1, x1, w1, n2, x2, and w2, while the evaluation point is specified by ii and mm.
      implicit none
#include "errquit.fh"
      integer i
      integer n1             ! [input]. The number of DVR points in the i-th FE.
      double precision x1(*) ! [input]. The DVR quadrature points in the i-th FE. Needed when m /= mm. When needed, its minimum size should be n1.
      double precision w1    ! [input]. The (i,n1)-th quadrature weight.
      integer n2             ! [input]. The number of DVR points in the (i+1)-th FE.
      double precision x2(*) ! [input]. The DVR quadrature points in the (i+1)-th FE. Needed when m /= mm. When needed, its minimum size should be n2.
      double precision w2    ! [input]. The (i+1,1)-th quadrature weight.
      integer ii             ! [input]. The FE index of the evaluation point.
      integer mm             ! [input]. The DVR index of the evaluation point.
      double precision f     ! [input]. The value of the bridge FEDVR basis at the (i,mm)-th quadrature point.

      double precision f1, f2


      if (ii == i) then
         if (mm < 1 .or. mm > n1)
     &        call errquit('WFAT: Since ii=i, the argument mm of  '//
     &        'wfat_eval1_fedvr_bridge must be in between '//
     &        '1 and n1.', 1, UERR)
      else if (ii == i+1) then
         if (mm < 1 .or. mm > n2)
     &        call errquit('WFAT: Since ii=i+1, the argument mm of  '//
     &        'wfat_eval1_fedvr_bridge must be in between '//
     &        '1 and n2.', 1, UERR)
      endif
         
      
c     ==== i side ====
      if (ii == i) then
         if (mm /= n1) then
            call wfat_eval1_lobatto_dm(n1, n1, mm, x1, f1)
         else
            call wfat_eval1_lobatto_em(n1, n1, w1, f1)
         endif
         f2 = 0.0d0
c     ==== i+1 side ====
      else if (ii == i+1) then
         f1 = 0.0d0
         if (mm /= 1) then
            call wfat_eval1_lobatto_dm(1, n2, mm, x2, f2)
         else
            call wfat_eval1_lobatto_em(1, n2, w2, f2)
         endif
      else
         f1 = 0.0d0
         f2 = 0.0d0
      endif
      

      f = (f1 + f2) / sqrt(w1 + w2)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval1_lobatto_dm(m, n, mm, x, f)
c     Different m.
      implicit none
#include "errquit.fh"
      integer m
      integer n
      integer mm
      double precision x(n)
      double precision f

      integer k
      double precision a, b

      if (m < 1 .or. m > n)
     &     call errquit('WFAT: The argument m of '//
     &     'wfat_eval1_lobatto_dm must be in between 1 and n.', 1, UERR)
      if (mm < 1 .or. mm > n)
     &     call errquit('WFAT: The argument mm of '//
     &     'wfat_eval1_lobatto_dm must be in between 1 and n.', 1, UERR)
      
      a = 1.0d0 / (x(m) - x(mm))

      b = 1.0d0
      do k = 1, n
         if (k /= m .and. k /= mm) b = b * (x(mm)-x(k)) / (x(m)-x(k))
      enddo
      f = a * b

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_eval1_lobatto_em(m, n, w, f)
c     Equal m.
      implicit none
#include "errquit.fh"
      integer m
      integer n
      double precision w
      double precision f

      integer delkron
      external delkron

      if (m < 1 .or. m > n)
     &     call errquit('WFAT: The argument m of '//
     &     'wfat_eval1_lobatto_em must be in between 1 and n.', 1, UERR)
      
      f = 1.0d0 / (2*w) * (delkron(m,n) - delkron(m,1))
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval00_fedvr(m, n, x, w, y, f)

      implicit none
#include "errquit.fh"
      integer m
      integer n
      double precision x(n)
      double precision w
      double precision y
      double precision f

      if (m < 1 .or. m > n)
     &     call errquit('WFAT: The argument m of '//
     &     'wfat_eval1_fedvr must be in between 1 and n.', 1, UERR)
      
      call wfat_eval00_lobatto(m, n, x, y, f)
      f = f / sqrt(w)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_eval00_fedvr_bridge(n1, x1, w1, n2, x2, w2, y, f)

      implicit none
#include "errquit.fh"
      integer n1             ! [input]. The number of DVR points in the i-th FE.
      double precision x1(*) ! [input]. The DVR quadrature points in the i-th FE. Needed when m /= mm. When needed, its minimum size should be n1.
      double precision w1    ! [input]. The (i,n1)-th quadrature weight.
      integer n2             ! [input]. The number of DVR points in the (i+1)-th FE.
      double precision x2(*) ! [input]. The DVR quadrature points in the (i+1)-th FE. Needed when m /= mm. When needed, its minimum size should be n2.
      double precision w2    ! [input]. The (i+1,1)-th quadrature weight.
      double precision y     ! [input]. The evaluation point.
      double precision f     ! [input]. The value of the bridge FEDVR basis at the (i,mm)-th quadrature point.

      double precision f1, f2

      if (abs(x1(n1)-x2(1)) < 1.0d-12)
     &     call errquit('WFAT: In wfat_eval0_fedvr_bridge, '//
     &     'the last point of the i-th element is too far from the '//
     &     'first point of the (i+1)-th element.', 1, UERR)
      
      call wfat_eval00_lobatto(n1, n1, x1, y, f1)
      call wfat_eval00_lobatto(1, n2, x2, y, f2)

      f = (f1 + f2) / sqrt(w1 + w2)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval00_lobatto(m, n, x, y, f)

      implicit none
#include "errquit.fh"
      integer m
      integer n
      double precision x(n)
      double precision y
      double precision f

      integer j

      if (m < 1 .or. m > n)
     &     call errquit('WFAT: The argument m of '//
     &     'wfat_eval0_lobatto must be in between 1 and n.', 1, UERR)

      if (y < x(1) .or. y > x(n)) then
         f = 0.0d0
         return
      endif

      f = 1.0d0
      do j = 1, n
         if (j /= m) f = f * (y-x(j)) / (x(m)-x(j))
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval0_fedvr(ni, n, maxn, xq, wq, nb, y, scr, f)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer ni
      integer n(ni)
      integer maxn
      double precision xq(maxn,ni)
      double precision wq(maxn,ni)
      integer nb
      double precision y
      double precision scr(*)   ! Its minimum size is maxn*maxn + maxn.
      double precision f(nb)

      logical incl
      integer j, i, nn, pw, pf, m
      double precision f_


c      f(1:nb) = 0.0d0
      j = 0
      do i = 1, ni
         nn = n(i)
         pw = 1
         pf = pw + nn**2
         call wfat_eval0_lobatto(nn, xq(:,i), y, scr(pw), scr(pf))
c         write(luout, '("scr_pf: ", i5)') i
c         call wfat_dprintmat(1, nn, 'e', scr(pf))


         do m = 1, nn
            if (i >= 2 .and. m == 1)
     &           f(j) = (f_+scr(pf)) / sqrt(wq(n(i-1),i-1) + wq(1,i))
               
            incl = .true.
            if ((i==1 .and. m==1) .or. (i==ni .and. m==n(ni)) .or. m==1)
     &           incl = .false.
            if (incl) then
               j = j + 1
               if (m >= 2 .and. m <= nn-1) then
                  f(j) = scr(pf+m-1) / sqrt(wq(m,i))
               else if (i <= ni-1 .and. m == nn) then
                  f_ = scr(pf+nn-1)
               endif
            endif
         enddo
         
cc         do m = 1, nn
cc            if (i > 1 .and. m == 1) then
cc               f(j) = (f_ + scr(pf)) / sqrt(wq(nn,i-1) + wq(1,i))
cc               write(luout, '("f_ : ", 3i5, 4es14.6)')
cc     &              i, m, j, f_, scr(pf), wq(nn,i-1), wq(1,i)
cc            else if (m > 1 .and. m < nn) then
cc               f(j) = scr(pf+m-1) / sqrt(wq(m,i))
cc            else if (i < ni .and. m == nn) then
cc               f_ = scr(pf+nn-1)
cc            endif
cc
cc            incl = .true.
cc            if ((i==1 .and. m==1) .or. (i==ni .and. m==n(ni)) .or. m==1)
cc     &           incl = .false.
cc            if (incl) then
cc               write(luout, '("i,m,j = ", 3i5, es14.6)')
cc     &              i, m, j, f(j)
cc               j = j + 1
cc            endif
cc         enddo
         
      enddo

      if (j /= nb) then
         write(luout, '("j = ", i4, 4x, "nb = ", i4)') j, nb+1
         call errquit('WFAT: j /= nb, this is a bug.', 74, UERR)
      endif
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_eval0_lobatto(n, x, y, scr, f)

      implicit none
#include "stdio.fh"
      integer n
      double precision x(n)
      double precision y
      double precision scr(*)   ! minimum size is n*n.
      double precision f(n)

      integer j, m, i_w
      double precision lx

      do m = 1, n
         if (abs(y - x(m)) < 1.0d-12) then
            f = 0.0d0
            f(m) = 1.0d0
c            write(luout, '("match m")')
c           == Return if the evaluation point is ==
c           ==   equal to one of the DVR grids   ==
            return
         endif
      enddo
      
c     ==== If the evaluation point is not in the same FE ====
      if (y < x(1) .or. y > x(n)) then
         f = 0.0d0
         return
      endif

c      write(luout, '("herev1")')
      lx = 1.0d0
      do m = 1, n
         lx = lx * (y - x(m))
      enddo

c      write(luout, '("herev2")')
      call wfat_barycentric_w(n-1, x, scr)
      do m = 1, n
         i_w = 1 + (n-1)*n + (m-1)
         f(m) = lx * scr(i_w) / (y-x(m))
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_barycentric_w(n, x, w)
c     Implements the algorithm in Section 3 of Ref. [1] to compute the barycentric weights.
      implicit none
#include "stdio.fh"
      integer n
      double precision x(0:n)
      double precision w(0:n,0:n)

      integer j, k

c      write(luout, '("heres0")')
      w(0,0) = 1.0d0
      do j = 1, n
         do k = 0, j-1
            w(k,j) = (x(k)-x(j)) * w(k,j-1)
         enddo
         w(j,j) = 1.0d0
         do k = 0, j-1
            w(j,j) = w(j,j) * (x(j)-x(k))
         enddo
      enddo

      do j = 0, n
         w(j,n) = 1 / w(j,n)
      enddo
      
      end
c     REFERENCES:
c     [1] SIAM REVIEW Vol. 46, No. 3, pp. 501–517.
c=======================================================

      
