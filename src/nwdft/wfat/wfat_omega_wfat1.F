c=======================================================
c=======================================================
      subroutine wfat_omega_wfat1(geom, ozn, nch, nqn, erg, chpos, nid, 
     &     mid, lmax, nat, nptp, centid, ptlim, gq, nlgl, zq, wfat0, 
     &     omg1, omg2)

      implicit none
#include "util.fh"
#include "global.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                                         ! [input],
      integer ozn                                          ! [input],
      integer nch                                          ! [input],
      integer nqn                                          ! [input],
      double precision erg(nch)                            ! [input],
      double precision chpos(3,nch)                        ! [input],
      integer nid(nqn)                                     ! [input],
      integer mid(nqn)                                     ! [input],
      integer lmax                                         ! [input],
      integer nat                                          ! [input],
      integer nptp                                         ! [input],
      integer centid(nat)                                  ! [input],
      integer ptlim(2,nat)                                 ! [input],
      integer gq                                           ! [input],
      integer nlgl                                         ! [input],
      double precision zq                                  ! [input],
      logical wfat0(nqn,nch)                               ! [input],
      double complex omg1((1+lmax)**2, nqn, nch, nptp)     ! [output],
      double complex omg2((1+lmax)**2, nqn, nch, nptp)     ! [output],

      character*12 dt_s
      integer l_r, l_t, l_p, l_sp, l_id, l_q1, l_q2
      integer k_r, k_t, k_p, k_sp, k_id, k_q1, k_q2
      integer ndid, nsh, i, ii, j, n, m, nl, i_sp, k
      double precision t_now, t_last, t_start, dt_all, dt_prog
      double complex czero

      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


      ndid = ga_nodeid()
      if (.not. MA_alloc_get(MT_DBL, nptp, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 't', l_t, k_t))
     &     call errquit('WFAT: t allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 'p', l_p, k_p))
     &     call errquit('WFAT: p allocation failed', 59, MA_ERR)
      nsh = (lmax + 1)**2 * nptp
      if (.not. MA_alloc_get(MT_DCPL, nsh, 'sphm', l_sp, k_sp))
     &     call errquit('WFAT: sphm allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_INT, nptp, 'indx', l_id, k_id))
     &     call errquit('WFAT: indx allocation failed', 59, MA_ERR)
      

      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Computation of Q '//
     &        'functions', 50, .true.)
         write(luout, '("")')
         write(luout, 901)
         write(luout, '("")')
      endif
 901  format(
     &     5x, 'The computation of Q functions needed by WFAT(1) ',
     &     'has started. This part can take up to several minutes, ', /,
     &     5x, 'depending on the number of channels, parabolic ',
     &     'quantum numbers, maximum angular momentum, and integration',
     &      /, 5x,'grids, so please wait. The progress of ',
     &     'this calculation will be shown below.')

      
      t_last = util_cpusec()
      t_start = t_last
      do i = 1, nch
c        ==== Get the grid points with the current channel center ====
         call wfat_get_gridpt_cents(geom, ozn, nat, centid, ptlim, 
     &        chpos(:,i), gq, dbl_mb(k_r), dbl_mb(k_t), dbl_mb(k_p))
         do k = 1, nptp
            int_mb(k_id+k-1) = k
         enddo
         call wfat_dsortx(dbl_mb(k_r), 1, nptp, int_mb(k_id))

         
c        The code below will calculate the WFAT(1) parabolic channel functions
c        regardless of whether the evaluation point belongs to a nuclei, parabolic
c        channel center with grid, or gridless parabolic channel center. The
c        reason the last one is not skipped is that it may help in the integration
c        involved to calculate the WFAT(1) parabolic channel functions at points
c        with a greater radial distance than that of the current evaluation point
c        that belongs to a parabolic channel center.
         
c        ==== Calculate the spherical harmonics for all (l,m) ====
c        ====   where l=[0,lmax] and m=[-l,l] at each point   ====
         do ii = 1, nptp
            i_sp = k_sp + (ii-1)*(1+lmax)**2
            call wfat_sphm_1pt(lmax, dbl_mb(k_t+ii-1), dbl_mb(k_p+ii-1),
     &           dcpl_mb(i_sp))
         enddo

         do j = 1, nqn
            n = nid(j)
            m = mid(j)
            
            if (wfat0(j,i)) then
               omg1(:,j,i,:) = czero
               omg2(:,j,i,:) = czero
            else
               nl = lmax - abs(m) + 1
               if (.not. MA_alloc_get(MT_DBL, nl*nptp, 'q1',
     &              l_q1, k_q1)) call errquit('WFAT: q1 allocation '//
     &              'failed', 94, MA_ERR)
               if (.not. MA_alloc_get(MT_DBL, nl*nptp, 'q2',
     &              l_q2, k_q2)) call errquit('WFAT: q2 allocation '//
     &              'failed', 94, MA_ERR)
               
               
c              ==== Calculate Q_{(l,m)}^\nu at all points above ====
               call wfat_calc_q(nptp, dbl_mb(k_r), n, m, lmax, nlgl, 
     &              erg(i), zq, dbl_mb(k_q1), dbl_mb(k_q2))
               call wfat_reorder_q(lmax, m, nptp, int_mb(k_id),
     &              dbl_mb(k_q1), dbl_mb(k_q2))
               
               
c              ==== Calculate WFAT(1) Omega for all  ====
c              ====    (l,m) and all points above    ====
               call wfat_omega_wfat1_nm(nptp, lmax, m, dbl_mb(k_q1), 
     &              dbl_mb(k_q2), dcpl_mb(k_sp), omg1(:,j,i,:),
     &              omg2(:,j,i,:))
               
               
               if (.not. MA_free_heap(l_q1)) call errquit('WFAT: q1 '//
     &              'deallocation failed', 32, MA_ERR)
               if (.not. MA_free_heap(l_q2)) call errquit('WFAT: q2 '//
     &              'deallocation failed', 32, MA_ERR)
            endif

c           ==== Print timing ====
            t_now = util_cpusec()
            dt_prog = t_now - t_last
            t_last = t_now
            call wfat_omega_wfat1_progress(i, n, m, dt_prog)
            call ga_sync() 
         enddo
      enddo
      dt_all = t_now - t_start
      write(dt_s, '(f12.2)') dt_all
      if (ndid == 0) write(luout, '(5x, "The computation of Q '//
     &     'functions for WFAT(1) takes ", a, " s.")')
     &     trim(adjustl(dt_s))

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_t))
     &     call errquit('WFAT: t deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_p))
     &     call errquit('WFAT: p deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_sp))
     &     call errquit('WFAT: sphm deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_id))
     &     call errquit('WFAT: indx deallocation failed', 61, MA_ERR)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_omega_wfat1_nm(nptp, lmax, m, q_1, q_2, sphm,
     &     omg1, omg2)

c     Calculate omg1 and omg2 for all (l,m) and all points at a particular channel 
c     and parabolic quantum number.
      implicit none
#include "stdio.fh"
      integer nptp
      integer lmax
      integer m
      double precision q_1(abs(m):lmax, nptp)       ! [input], orientation-independent
      double precision q_2(abs(m):lmax, nptp)       ! [input], orientation-dependent
      double complex sphm((1+lmax)**2, nptp)        ! [input]
      double complex omg1((1+lmax)**2, nptp)        ! [output]
      double complex omg2((1+lmax)**2, nptp)        ! [output]
      
      integer ii, l, m_, i_sh
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )

      omg1 = czero
      omg2 = czero
      do ii = 1, nptp
         do l = abs(m), lmax
            do m_ = -l, l
               i_sh = l**2 + l + m_ + 1
               omg1(i_sh, ii) = q_1(l, ii) * sphm(i_sh, ii)
               omg2(i_sh, ii) = q_2(l, ii) * sphm(i_sh, ii)

               if ( isnan(abs( omg1(i_sh,ii) )) ) omg1(i_sh,ii) = czero
               if ( isnan(abs( omg2(i_sh,ii) )) ) omg2(i_sh,ii) = czero
            enddo
         enddo
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_omega_wfat1_progress(ich, n, m, dt)

      implicit none
#include "global.fh"
#include "stdio.fh"
      integer ich
      integer n
      integer m
      double precision dt

      character*12 i_s, n_s, m_s, dt_s

      write(i_s, '(i12)') ich
      write(n_s, '(i12)') n
      write(m_s, '(i12)') m
      write(dt_s, '(f12.2)') dt

      if (ga_nodeid() == 0) then
         write(luout, '(7x, "Channel: ", a, ",  n = ", a, ",  m = ",'//
     &        'a, ",  time needed = ", a, " s")')
     &        trim(adjustl(i_s)), trim(adjustl(n_s)),
     &        trim(adjustl(m_s)), trim(adjustl(dt_s))
      endif
      
      end
c=======================================================
      
      
