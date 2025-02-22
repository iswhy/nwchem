c============================================================
c============================================================
      subroutine wfat_print_wfat0_ints(nch, nqn, nat, lmax, nid, mid,
     &     intl)

      implicit none
#include "stdio.fh"
      integer nch
      integer nqn
      integer nat
      integer lmax
      integer nid(nqn)
      integer mid(nqn)
      double complex intl((1+lmax)**2, nat, nqn, nch)       ! [input],

      character*12 i_s, n_s, m_s
      integer i, iqn, l, l_in, l_fin, n_lgroup, ng, g, im, m, i_sh
      parameter(n_lgroup=4)

      
      write(luout, '(2/)')
      call util_print_centered(luout, 'Spherical harmonics '//
     &     'components of WFAT integrals', 53, .true.)
      write(luout, '("")')

      write(luout, 480)
 480  format(
     &     5x, 'Use the tables below to determine the convergence of ',
     &     'your WFAT calculation with repsect to the maximum ', /,
     &     5x, 'angular momentum (lmax). In particular, for a given ',
     &     'channel and parabolic quantum numbers, the integral ', /,
     &     5x, 'values in the last table (where l=lmax) must be ',
     &     'sufficiently small.')
      write(luout, '("")')

         
      do i = 1, nch
         do iqn = 1, nqn

            write(i_s, '(i12)') i
            write(n_s, '(i12)') nid(iqn)
            write(m_s, '(i12)') mid(iqn)
            write(luout, '(5x, "Channel: ", a, ",  n = ", a, '//
     &           '",  m = ", a)') trim(adjustl(i_s)), 
     &           trim(adjustl(n_s)), trim(adjustl(m_s))
            write(luout, '("")')

            l_in = abs(mid(iqn))
            l_fin = min(l_in + n_lgroup - 1, lmax)
            ng = int((lmax+1)/n_lgroup) + 1
            do g = 1, ng

c              ==== Column title ====
               write(luout, '(5x, "")', advance='no')    ! Indentation
               do l = l_in, l_fin
                  write(luout, '(2a4)', advance='no') 'l', 'm'
                  write(luout, '(1x, 2a12)', advance='no')
     &                 'I (Re.)', 'I (Im.)'
                  if (l < l_fin) then
                     write(luout, '("   ")', advance='no')    ! hori. space with the next l in the same group.
                  else if (l == l_fin) then
                     write(luout, '("")')
                  endif
               enddo

c              ==== Horizontal lines ====
               write(luout, '(5x, "")', advance='no')    ! Indentation
               do l = l_in, l_fin
                  call wfat_print_repeat(luout, 33, 0, '-')
                  if (l < l_fin) then
                     write(luout, '("   ")', advance='no')    ! hori. space with the next l in the same group.
                  else if (l == l_fin) then
                     write(luout, '("")')
                  endif
               enddo

c              ==== Integral values ====
               do im = 1, (2*l_fin+1)
                  write(luout, '(5x, "")', advance='no')    ! Indentation
                  do l = l_in, l_fin
                     m = im - 1 - l
                     if (m <= l) then
                        i_sh = l**2 + l + m + 1
                        write(luout, '(2i4)', advance='no') l, m
                        write(luout, '(1x, 2f12.8)', advance='no')
     &                       sum(intl(i_sh, :, iqn, i))
                     else
                        write(luout, '(33x, "")', advance='no')    ! blank space when m exceeds +l.
                     endif

                     if (l < l_fin) then
                        write(luout, '("   ")', advance='no')    ! hori. space with the next l in the same group.
                     else if (l == l_fin) then
                        write(luout, '("")')      ! move to the next line
                     endif
                  enddo
               enddo

               l_in = l_in + n_lgroup
               l_fin = min(l_fin + n_lgroup, lmax)
               write(luout, '("")')      ! vert. space with the next group
            enddo
            write(luout, '("")')    ! vert. space with the next (i,iqn)

         enddo
      enddo

      end
c============================================================
      

c============================================================
c============================================================
      subroutine wfat_print_wfat1_ints(nch, nqn, nat, lmax, nid, mid,
     &     wfat1, intl, intl_j, intl_k1, intl_k2)

      implicit none
#include "stdio.fh"
      integer nch
      integer nqn
      integer nat
      integer lmax
      integer nid(nqn)
      integer mid(nqn)
      logical wfat1
      double complex intl((1+lmax)**2, nat, nqn, nch)       ! [input],
      double complex intl_j(1:*)     ! [input], only used when wfat1=.true.
      double complex intl_k1(1:*)    ! [input], only used when wfat1=.true.
      double complex intl_k2(1:*)    ! [input], only used when wfat1=.true.

      character*12 i_s, n_s, m_s
      integer i, iqn, l, m, i_sh, j


      write(luout, '(2/)')
      call util_print_centered(luout, 'The spherical harmonics '//
     &     'components of WFAT integrals', 53, .true.)
      write(luout, '("")')

      write(luout, 480)
 480  format(
     &     5x, 'Use the tables below to determine the convergence of ',
     &     'your WFAT calculation with repsect to the maximum ', /,
     &     5x, 'angular momentum (lmax). In particular, for a given ',
     &     'channel and parabolic quantum numbers, the integral ', /,
     &     5x, 'values in the last table (where l=lmax) must be ',
     &     'sufficiently small.')
      write(luout, '("")')

         
      do i = 1, nch
         do iqn = 1, nqn

            write(i_s, '(i12)') i
            write(n_s, '(i12)') nid(iqn)
            write(m_s, '(i12)') mid(iqn)
            write(luout, '(5x, "Channel: ", a, ",  n = ", a, '//
     &           '",  m = ", a)') trim(adjustl(i_s)), 
     &           trim(adjustl(n_s)), trim(adjustl(m_s))

            call wfat_print_repeat(luout, 2*4+3+2*12+5*(3+2*12), 5, '-')
            write(luout, '("")')
            write(luout, '(5x, 2a4)', advance='no') 'l', 'm'
            write(luout, '(3x, 2a12)', advance='no')
     &           'I (Re.)', 'I (Im.)'
            
            if (wfat1) then
               write(luout, ('(3x, 2a12)'), advance='no')
     &              'J_x (Re.)', 'J_x (Im.)'
               write(luout, ('(3x, 2a12)'), advance='no') 
     &              'J_y (Re.)', 'J_y (Im.)'
               write(luout, ('(3x, 2a12)'), advance='no') 
     &              'J_z (Re.)', 'J_z (Im.)'
               write(luout, ('(3x, 2a12)'), advance='no') 
     &              'K_1 (Re.)', 'K_1 (Im.)'
               write(luout, ('(3x, 2a12)'), advance='no') 
     &              'K_2 (Re.)', 'K_2 (Im.)'
            endif
            write(luout, '("")')
            call wfat_print_repeat(luout, 2*4+3+2*12+5*(3+2*12), 5, '-')
            write(luout, '("")')


            do l = abs(mid(iqn)), lmax
               do m = -l, l
                  i_sh = l**2 + l + m + 1
                  write(luout, '(5x, 2i4)', advance='no') l, m
                  write(luout, '(3x, 2f12.6)', advance='no')
     &                 sum(intl(i_sh, :, iqn, i))

                  if (wfat1)
     &                 call wfat_print_wfat1_ints_(nch, nqn, nat, lmax,
     &                 i_sh, iqn, i, intl_j, intl_k1, intl_k2)

                  write(luout, '(" ")')        ! new line
               enddo
               call wfat_print_repeat(luout, 2*4+3+2*12+5*(3+2*12), 5,
     &              '-')
               write(luout, '("")')
            enddo
            write(luout, '(" ")')       ! gap between different (i,iqn)
            write(luout, '(" ")')       ! gap between different (i,iqn)

         enddo
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_print_wfat1_ints_(nch, nqn, nat, lmax, i_sh, i_qn,
     &     i_ch, intl_j, intl_k1, intl_k2)

      implicit none
#include "stdio.fh"
      integer nch
      integer nqn
      integer nat
      integer lmax
      integer i_sh
      integer i_qn
      integer i_ch
      double complex intl_j((1+lmax)**2, nat, nqn, 3, nch)  ! [input],   only used when wfat1=.true.
      double complex intl_k1((1+lmax)**2, nat, nqn, nch)    ! [input],   only used when wfat1=.true.
      double complex intl_k2((1+lmax)**2, nat, nqn, nch)    ! [input],   only used when wfat1=.true.

      integer j

      do j = 1, 3
         write(luout, '(3x, 2f12.6)', advance='no')
     &        sum(intl_j(i_sh, :, i_qn, j, i_ch))
      enddo
      write(luout, '(3x, 2f12.6)', advance='no')
     &     sum(intl_k1(i_sh, :, i_qn, i_ch))
      write(luout, '(3x, 2f12.6)', advance='no')
     &     sum(intl_k2(i_sh, :, i_qn, i_ch))

      end
c============================================================
      
