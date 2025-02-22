c==========================================================
c==========================================================
      subroutine wfat_2a1e_stateinfo(m, sphi, spin, orbid)

      implicit none
#include "stdio.fh"
      integer m, sphi, spin, orbid

      character*8, sm, so
      character*255 ss, xzr


      write(luout, '(2/)')
      call util_print_centered(luout, '2a1e State Information',
     &     23, .true.)
      write(luout, '("")')

      write(sm, '(i8)') m
      if (sphi == 1) then
         xzr = 'symmetric'
      else if (sphi == 2) then
         xzr = 'antisymmetric'
      endif
      if (spin == 1) then
         ss = 'alpha'
      else if (spin == 2) then
         ss = 'beta'
      endif
      write(so, '(i8)') orbid
      
      write(luout, '(4x, "Magnetic quantum number = ", a)')
     &     trim(adjustl(sm))
      write(luout, '(4x, "xz-plane reflection symmetry = ", a)')
     &     trim(adjustl(xzr))
      write(luout, '(4x, "Spin = ", a)') trim(adjustl(ss))
      write(luout, '(4x, "Energy eigenstate ID = ", a)')
     &     trim(adjustl(so))

      write(luout, '("")')
      
      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_2a1efd_print_rgrid(ni, n, maxn, xq, wq)

      implicit none
#include "stdio.fh"
      integer ni
      integer n(ni)
      integer maxn
      double precision xq(maxn,ni)
      double precision wq(maxn,ni)

      character*256 t
      character*14 sf
      character*9 si
      logical incl
      integer k, ife, ig, nline, maxline, nr
      parameter(maxline=110)


      write(luout, '(2/)')
      call util_print_centered(luout, 'Radial DVR grids and weights',
     &     58, .true.)
      write(luout, '("")')
      write(luout, '(4x, "DVR type = FEDVR")')
      call wfat_get_fedvrnr(ni, n, nr)
      write(luout, '(4x, "No. of radial DVR bases = ", i4)') nr
      write(luout, '("")')


c     ==== Print DVR grids ====
      write(luout, '(4x, a)')
     &     'Radial DVR grids (values inside parentheses correspond '//
     &     'to a bridge point) = '
      k = 0
      nline = 0
      write(luout, '(6x)', advance='no')
      do ife = 1, ni
         do ig = 2, n(ife)

            incl = .true.
            if ((ife==1 .and. ig==1) .or. (ife==ni .and. ig==n(ni))
     &           .or. ig==1) incl = .false.
            if (incl) then
               k = k + 1

               write(si, '(i9)') k
               if (xq(ig, ife) < 100.0d0) then
                  write(sf, '(f9.5)') xq(ig, ife)
               else
                  write(sf, '(es14.5)') xq(ig, ife)
               endif

               if (ig == n(ife)) then
                  write(t, '(2x, a)') "("//trim(adjustl(si))//'. '//
     &                 trim(adjustl(sf))//")"
               else
                  write(t, '(2x, a)') " "//trim(adjustl(si))//'. '//
     &                 trim(adjustl(sf))//" "
               endif
               nline = nline + len(trim(t))
               if (nline > maxline) then
                  write(luout, '("")')
                  write(luout, '(6x)', advance='no')
                  nline = len(trim(t))
               endif

               write(luout, '(a)', advance='no') trim(t)
            endif

         enddo
      enddo

      write(luout, '("")')
      write(luout, '("")')

      
c     ==== Print DVR weights ====
      write(luout, '(4x, a)')
     &     'Radial DVR weights (values inside parentheses correspond '//
     &     'to the first, last, or a bridge point) = '
      k = 0
      nline = 0
      write(luout, '(6x)', advance='no')
      do ife = 1, ni
         do ig = 1, n(ife)

            incl = .true.
            if ((ife==1 .and. ig==1) .or. (ife==ni .and. ig==n(ni))
     &           .or. ig==1) incl = .false.
            if (incl) then
               k = k + 1
            endif
            
            write(si, '(i9)') k
            write(sf, '(es14.5)') wq(ig, ife)

            if (ig == 1) then
               if (ife == 1) then
                  write(t, '(2x, a)') "("//trim(adjustl(sf))//")"
               else                  
                  write(t, '(2x, a)') " "//trim(adjustl(si))//'. '//
     &                 trim(adjustl(sf))//")"
               endif
            else if (ig == n(ife)) then
               if (ife == ni) then
                  write(t, '(2x, a)') "("//trim(adjustl(sf))//")"
               else 
                  write(t, '(2x, a)') "("//trim(adjustl(si))//'. '//
     &                 trim(adjustl(sf))//" "
               endif
            else
               write(t, '(2x, a)') " "//trim(adjustl(si))//'. '//
     &              trim(adjustl(sf))//" "
            endif
            nline = nline + len(trim(t))
            if (nline > maxline) then
               write(luout, '("")')
               write(luout, '(6x)', advance='no')
               nline = len(trim(t))
            endif

            write(luout, '(a)', advance='no') trim(t)

         enddo
      enddo

      write(luout, '("")')
      
      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_2a1e_print_rgrid(m, nr)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer m
      integer nr

      character*256 t
      character*14 sf
      character*9 si
      integer l_x, k_x, l_rw, k_rw, i, nline, maxline
      double precision xq
      parameter(maxline = 110)

      
      write(luout, '(2/)')
      call util_print_centered(luout, 'Radial DVR grids and weights',
     &     58, .true.)
      write(luout, '("")')
      write(luout, '(4x, "DVR type = Generalized Laguerre")')
      write(luout, '(4x, "No. of radial DVR bases (subscript) = ", i4)')
     &     nr
      write(luout, '(4x, "Superscript index = ", i4)') m + 2
      write(luout, '("")')


      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'xweight', l_rw, k_rw))
     &     call errquit('WFAT: Cannot allocate xweight', 1, MA_ERR)
      call wfat_lf_quadrature_rule(nr, abs(dble(m))+2, dbl_mb(k_x),
     &     dbl_mb(k_rw))


c     ==== Print DVR grids ====
      write(luout, '(4x, a)') 'Radial DVR grids = '
      nline = 0
      write(luout, '(6x)', advance='no')
      do i = 1, nr

         write(si, '(i9)') i
         xq = dbl_mb(k_x+i-1)
         if (xq < 100.0d0) then
            write(sf, '(f9.5)') xq
         else
            write(sf, '(es14.5)') xq
         endif

         write(t, '(2x, a)') trim(adjustl(si))//'. '//trim(adjustl(sf))
         nline = nline + len(trim(t))
         if (nline > maxline) then
            write(luout, '("")')
            write(luout, '(6x)', advance='no')
            nline = len(trim(t))
         endif

         write(luout, '(a)', advance='no') trim(t)
      enddo

      write(luout, '("")')
      write(luout, '("")')

      
c     ==== Print DVR weights ====
      write(luout, '(4x, a)') 'Radial DVR weights = '
      nline = 0
      write(luout, '(6x)', advance='no')
      do i = 1, nr
         write(si, '(i9)') i
         write(sf, '(es14.5)') dbl_mb(k_rw+i-1)
         write(t, '(2x, a)') trim(adjustl(si))//'. '//trim(adjustl(sf))

         nline = nline + len(trim(t))
         if (nline > maxline) then
            write(luout, '("")')
            write(luout, '(6x)', advance='no')
            nline = len(trim(t))
         endif

         write(luout, '(a)', advance='no') trim(t)
      enddo

      write(luout, '("")')      

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_x))
     &     call errquit('WFAT: xquad deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_rw))
     &     call errquit('WFAT: xweight deallocation failed', 3, MA_ERR)
      
      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_2a1e_print_tgrid(m, nt)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer m
      integer nt

      character*256 t
      character*14 sf
      character*9 si
      integer l_t, k_t, l_tw, k_tw, i, nline, maxline
      double precision deg2rad
      parameter(maxline = 110, deg2rad=0.01745329251994d0)


      write(luout, '(2/)')
      call util_print_centered(luout, 'Angular DVR grids and weights',
     &     58, .true.)
      write(luout, '("")')
      write(luout, '(4x, "DVR type = Associated Legendre")')
      write(luout, '(4x, "No. of angular DVR bases = ", i4)') nt
      write(luout, '(4x, "Magnetic quantum number = ", i4)') m
      write(luout, '("")')


      if (.not. ma_alloc_get(MT_DBL, nt, 'tquad', l_t, k_t))
     &     call errquit('WFAT: Cannot allocate tquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt, 'tweight', l_tw, k_tw))
     &     call errquit('WFAT: Cannot allocate tweight', 1, MA_ERR)
      call wfat_assocleg_quad(-1.0d0, 1.0d0, abs(m), nt, dbl_mb(k_t),
     &     dbl_mb(k_tw))
      dbl_mb(k_t:k_t+nt-1) = acos(dbl_mb(k_t:k_t+nt-1)) / deg2rad


c     ==== Print DVR grids ====
      write(luout, '(4x, a)') 'Angular DVR grids = '
      nline = 0
      write(luout, '(6x)', advance='no')
      do i = 1, nt
         write(si, '(i9)') i
         write(sf, '(f9.5)') dbl_mb(k_t+i-1)
         write(t, '(2x, a)') trim(adjustl(si))//'. '//trim(adjustl(sf))

         nline = nline + len(trim(t))
         if (nline > maxline) then
            write(luout, '("")')
            write(luout, '(6x)', advance='no')
            nline = len(trim(t))
         endif

         write(luout, '(a)', advance='no') trim(t)
      enddo

      write(luout, '("")')
      write(luout, '("")')


c     ==== Print DVR weights ====
      write(luout, '(4x, a)') 'Angular DVR weights = '
      nline = 0
      write(luout, '(6x)', advance='no')
      do i = 1, nt
         write(si, '(i9)') i
         write(sf, '(es14.5)') dbl_mb(k_tw+i-1)
         write(t, '(2x, a)') trim(adjustl(si))//'. '//trim(adjustl(sf))

         nline = nline + len(trim(t))
         if (nline > maxline) then
            write(luout, '("")')
            write(luout, '(6x)', advance='no')
            nline = len(trim(t))
         endif

         write(luout, '(a)', advance='no') trim(t)
      enddo

      write(luout, '("")')      
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_t))
     &     call errquit('WFAT: tquad deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_tw))
     &     call errquit('WFAT: tweight deallocation failed', 3, MA_ERR)
      
      end
c==========================================================
      
      
