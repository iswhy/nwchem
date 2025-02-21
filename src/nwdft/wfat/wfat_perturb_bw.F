c===============================================================
c===============================================================
      subroutine wfat_perturb_bw(geom, bas, nb, ns, g_mo, erg0, g_hi, 
     &     nch, orbid, sid, dgnthr, offthr, e1thr, itrthr, maxitr, 
     &     ordthr, maxord, ord, numitr, numord, max_c, max_c_i, 
     &     erg, g_cmo)
      
      implicit none
#include "stdio.fh"      
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"      
      integer geom
      integer bas
      integer nb
      integer ns
      integer g_mo(ns)                       ! [input], all unperturbed MO's.
      double precision erg0(nb,ns)           ! [input], unperturbed energies of all MO's.
      integer g_hi                           ! [input], interaction Hamiltonian in standard GTO basis.
      integer nch
      integer orbid(nch)
      integer sid(nch)
      double precision dgnthr                ! [input], degeneracy threshold.
      double precision offthr                ! [input], off-diagonal threshold.
      double precision e1thr                 ! [input], threshold for 1st order to be used in getting guess energies.
      double precision itrthr                ! [input], threshold used for solving the implicit equation for corrected energies.
      integer maxitr                         ! [input], maximum no. iterations for solving the implicit equation for corrected energies.
      double precision ordthr                ! [input], threshold used for determining no. of significant orders.
      integer maxord                         ! [input], maximum order of the corrections used for numerically solving the corrected energies.
      integer ord                            ! [input], order of corrections of MO coefficients returned in g_cmo. Naturally, 'ord' must be smaller than'maxord'.
      integer numitr(nch)                    ! [output], no. of iterations (per itrthr or maxitr) to yield erg.
      integer numord(nch)                    ! [output], no. of significant orders (per ordthr or maxord) in the last iteration.
      double precision max_c(nch)            ! [output], 
      integer max_c_i(2,nch)                 ! [output], 
      double precision erg(0:maxord, nch)    ! [output], corrections of MO energies in Hartree.
      integer g_cmo                          ! [output], coefficients of corrected orbitals. The size should be [nb]x[nch*(ord+1)]

      character*9 movname(2)
      character*7 hiname(2)
      integer is, ndid, g_mo2(ns), g_mohi(ns)
      

      call ga_sync()
      ndid = ga_nodeid()


      movname = (/ 'mov_alpha', 'mov_beta' /)
      hiname = (/ 'h_int_a', 'h_int_b' /)
      
      if (ord > maxord)
     &     call errquit('WFAT: Argument ''ord'' of '//
     &     '''wfat_perturb_bw'' cannot be bigger than ''maxord''.',
     &     1, UERR)
      

c     ==== Generate new MO's based on degeneracy ====
      do is = 1, ns
         if (.not. ga_duplicate(g_mo(is), g_mo2(is), movname(is)))
     &        call errquit('WFAT: Cannot duplicate '//movname(is)//'.',
     &        1, GA_ERR)
         call ga_copy(g_mo(is), g_mo2(is))
      enddo      
      call wfat_diag_perturb(nb, ns, g_mo2, erg0, g_hi, dgnthr, offthr)
      

c     ==== Construct H_int in degeneracy-diagonialized ====
c     ====            MO's representation              ====
      do is = 1, ns
         if (.not. ga_duplicate(g_hi, g_mohi(is), hiname(is)))
     &        call errquit('WFAT: Cannot create '//hiname(is)//'.',
     &        1, GA_ERR)
         call ga_zero(g_mohi(is))
         call wfat_ga_vrepresent(g_mo2(is), g_hi, g_mohi(is))

cforcheck         if (ndid==0)
cforcheck     &        write(luout, '("Matrix of dipole interaction '//
cforcheck     &        'Hamiltonian in the adapted ", a, " representation :")')
cforcheck     &        trim(adjustl(movname(is)))
cforcheck         call ga_print(g_mohi(is))
      enddo


c     ==== Numerically solve for corrected energies ====
      call wfat_bw_solveerg(geom, bas, nb, ns, erg0, nch, orbid,
     &           sid, g_mohi, e1thr, itrthr, maxitr, ordthr, maxord, 
     &           numitr, numord, max_c, max_c_i, erg)


c     ==== Calculate corrected MO's ====
      call wfat_bw_calcwfn(nb, ns, g_mo2, erg0, nch, orbid,
     &     sid, sum(erg, 1), g_mohi, ord, g_cmo)


c     ==== Deallocations ====
      do is = 1, ns
         if (.not. ga_destroy(g_mo2(is)))
     &        call errquit('WFAT: Cannot destroy '//movname(is)//'.',
     &        1, GA_ERR)
         if (.not. ga_destroy(g_mohi(is)))
     &        call errquit('WFAT: Cannot destroy '//hiname(is)//'.',
     &        1, GA_ERR)
      enddo

      call ga_sync()
      
      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_bw_solveerg(geom, bas, nb, ns, erg0, nch, orbid,
     &           sid, g_mohi, e1thr, itrthr, maxitr, ordthr, maxord, 
     &           numitr, numord, max_c, max_c_i, erg)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"      
      integer geom                           ! [input], geometry handle.         
      integer bas                            ! [input], basis handle.          
      integer nb                             ! [input], number of bases in bas.
      integer ns
      double precision erg0(nb,ns)           ! [input], unperturbed energies in Hartree.
      integer nch                            ! [input], number of MO's of interest.
      integer orbid(nch)                     ! [input], index of MO whose energy corrections are computed.
      integer sid(nch)
      integer g_mohi(ns)                     ! [input], interaction Hamiltonian matrix represented in unperturbed (orthonormal) MO basis.
      double precision e1thr                 ! [input], threshold for 1st order to be used in getting guess energies.
      double precision itrthr                ! [input], threshold used for solving the implicit equation for corrected energies.
      integer maxitr                         ! [input], maximum no. iterations for solving the implicit equation for corrected energies.
      double precision ordthr                ! [input], threshold used for determining no. of significant orders.
      integer maxord                         ! [input], maximum order of the corrections.
      integer numitr(nch)                    ! [output], no. of iterations (per itrthr or maxitr) to yield erg.
      integer numord(nch)                    ! [output], no. of significant orders (per ordthr or maxord) in the last iteration.
      double precision max_c(nch)            ! [output], 
      integer max_c_i(2,nch)                 ! [output], 
      double precision erg(0:maxord, 1:nch)  ! [output], corrections of MO energies in Hartree.

      integer ich, iv, is, i, nord
      integer l_e1, k_e1, g_pwhi, g_ms, g_it
      double precision ein, ein_prev


      call ga_sync()
      
      if (maxord < 2)
     &     call errquit('WFAT: Argument ''maxord'' of '//
     &     '''wfat_bw_ergsum'' must be larger than or equal to 2',
     &     1, UERR)
      

c     ==== Obtain 1st order corrections ====
      if (.not. ma_alloc_get(mt_dbl, nch, 'e1', l_e1, k_e1))
     &     call errquit('WFAT: Cannot allocate e1', 1, MA_ERR)
      do ich = 1, nch
         iv = orbid(ich)
         is = sid(ich)
         call ga_get(g_mohi(is), iv, iv, iv, iv, dbl_mb(k_e1+ich-1), 1)
      enddo


c     ==== Initialization of scratch arrays ====
      if (.not. ga_create(mt_dbl, nb, 1, 'mask', 0, 0, g_ms))
     &     call errquit('WFAT: Cannot create mask', 1, GA_ERR)
      call ga_zero(g_ms)

      if (.not. ga_duplicate(g_mohi(1), g_pwhi, 'power_h_int'))
     &     call errquit('WFAT: Cannot create power_h_int', 1, GA_ERR)
      
      if (.not. ga_duplicate(g_mohi(1), g_it, 'itermat'))
     &     call errquit('WFAT: Cannot create itermat', 1, GA_ERR)
      

c     ==== Loop over channels ====
      do ich = 1, nch
         iv = orbid(ich)
         is = sid(ich)
         erg(0:1, ich) = (/ erg0(iv,is), dbl_mb(k_e1+ich-1) /)

c        ==== Initiate the guess energies ====         
         if (abs(dbl_mb(k_e1+ich-1)) > e1thr) then
            ein = erg0(iv,is) + dbl_mb(k_e1+ich-1)
         else
            ein = erg0(iv,is) + e1thr  !1.0d-2
         endif

c        ==== Iterate the sum ====
         do i = 1, maxitr
            ein_prev = ein
            call ga_copy(g_mohi(is), g_pwhi)
            call ga_copy(g_mohi(is), g_it)
            call wfat_bw_ergsum(geom, bas, nb, iv, ein, erg0(:,is), 
     &           g_pwhi, g_ms, g_it, maxord, ordthr, nord,
     &           max_c(ich), max_c_i(:,ich), erg(2:maxord, ich))

            ein = sum(erg(0:maxord,ich))
            numord(ich) = nord
            numitr(ich) = i
            if (abs(ein - ein_prev) <= itrthr) exit
         enddo

      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_e1))
     &     call errquit('WFAT: Cannot deallocate e1', 1, MA_ERR)
      if (.not. ga_destroy(g_pwhi))
     &     call errquit('WFAT: Cannot destroy power_h_int', 1, GA_ERR)
      if (.not. ga_destroy(g_ms))
     &     call errquit('WFAT: Cannot destroy mask', 1, GA_ERR)
      if (.not. ga_destroy(g_it))
     &     call errquit('WFAT: Cannot destroy itermat', 1, GA_ERR)
         

      call ga_sync()
      
      end
c===============================================================
      
      
c===============================================================
c===============================================================
      subroutine wfat_bw_ergsum(geom, bas, nb, iv, ein, e0, 
     &           g_pwhi, g_ms, g_it, maxord, convthr, nord,
     &           max_c, max_c_i, eout)

      implicit none
#include "mafdecls.fh"      
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
      integer geom                       ! [input], geometry handle.
      integer bas                        ! [input], basis handle.
      integer nb                         ! [input], number of bases in bas.
      integer iv                         ! [input], index of MO whose correction is sought.
      double precision ein               ! [input], Exact or guess energy in Hartree of the iv-th orbital.
      double precision e0(nb)            ! [input], unperturbed energies in Hartree of all MO's.
      integer g_pwhi                     ! [in/output], on input, it must contain the matrix of interaction Hamiltonian
c                                                       represented in (orthonormal) MO basis. A scratch array having the size of [nb]x[nb].
      integer g_ms                       ! [output], mask vector. A scratch array having the size of [nb]x[1].
      integer g_it                       ! [in/output], on input, it must contain the matrix of interaction Hamiltonian
c                                                       represented in (orthonormal) MO basis. A scratch array having the size of [nb]x[nb].
      integer maxord                     ! [input], maximum correction order, must be larger than or equal to 2.
      double precision convthr           ! [input], threshold to control convergence of the sum.
      integer nord                       ! [input], actual number of correction terms to reach convthr.
      double precision max_c             ! [output], 
      integer max_c_i(2)                 ! [output], 
      double precision eout(1:maxord-1)  ! [output], corrections to energy starting from 2nd order, eout(1), onward.

      integer i, j, g_tmp, ndid, l_tmp, k_tmp, i_tmp
      double precision msk, elast1, elast2


      call ga_sync()
      ndid = ga_nodeid()
      

      if (maxord < 2)
     &     call errquit('WFAT: Argument ''maxord'' of '//
     &     '''wfat_bw_ergsum'' must be larger than or equal to 2',
     &     1, UERR)
      
      
c     ==== Construct mask vector ====
      call ga_zero(g_ms)
      do i = 1, nb     ! Loop over rows
         if (i == iv) then
            msk = 0.0d0
         else
            msk = 1.0d0 / (ein - e0(i))
         endif
         if (ndid==0) call ga_put(g_ms, i, i, 1, 1, msk, 1)
      enddo
      

c     ==== Construct iterator matrix ====
      do i = 1, nb
         call ga_elem_multiply_patch(
     &        g_ms, (/ 1, 1 /), (/ nb, 1 /),
     &        g_it, (/ 1, i /), (/ nb, i /),
     &        g_it, (/ 1, i /), (/ nb, i /))
      enddo


c     ==== Obtain convergence property ====
      if (.not. ma_alloc_get(mt_dbl, nb*nb, 'tmp', l_tmp, k_tmp))
     &     call errquit('WFAT: Cannot allocate tmp', 1, UERR)
      call ga_get(g_it, 1, nb, 1, nb, dbl_mb(k_tmp), nb)
      max_c = 0.0d0
      do j = 1, nb
         do i = 1, nb
            i_tmp = k_tmp + (j-1)*nb + i - 1
            if (abs(dbl_mb(i_tmp)) >= max_c) then
               max_c = abs(dbl_mb(i_tmp))
               max_c_i = (/ i, j /)
            endif
         enddo
      enddo
      if (.not. ma_free_heap(l_tmp))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, UERR)

cforcheck      if (ndid==0) write(luout, '("Iterator matrix: ", i3, es14.6)')
cforcheck     &     iv, ein
cforcheck      call ga_print(g_it)

      
c     ==== Calculate the correction orders ====
      if (.not. ga_duplicate(g_it, g_tmp, 'tmp'))
     &     call errquit('WFAT: Cannot create tmp', 0, GA_ERR)
      call ga_zero(g_tmp)
      eout = 0.0d0
      do i = 2, maxord         
         call ga_dgemm('n', 'n', nb, nb, nb, 1.0d0, g_pwhi, g_it,
     &        0.0d0, g_tmp)
         call ga_copy(g_tmp, g_pwhi)
         call ga_get(g_pwhi, iv, iv, iv, iv, eout(i-1), 1)

         nord = i - 1      ! no. of significant orders 

         if (i == 2) then
            elast1 = convthr + 1.0d2
            elast2 = eout(i-1)
         else
            elast1 = elast2
            elast2 = eout(i-1)
         endif

cforcheck         if (ga_nodeid() == 0) 
cforcheck     &        write(luout, '(3x, "eout= ", es14.6)', advance='no')
cforcheck     &        eout(i-1)
         

         if (abs(elast1) <= convthr .and. abs(elast2) <= convthr) then     ! 1)
            exit
         else
ccheck            if (i == maxord) 
         endif
      enddo

cforcheck      if (ndid==0) write(luout, '(x)')
      
c     ==== Deallocations ====
      if (.not. ga_destroy(g_tmp))
     &     call errquit('WFAT: Cannot destroy tmp', 1, GA_ERR)

      call ga_sync()
      
      end
c     NOTE:
c     1) Use current and previous energies to tell convergence, this is to get
c        around the situation where energy oscillates w.r.t orders.
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_bw_calcwfn(nb, ns, g_mo, erg0, nch, orbid,
     &     sid, erg, g_mohi, ord, g_cmo)

      implicit none
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"      
      integer nb
      integer ns
      integer g_mo(ns)                ! [input], all unperturbed MO's, IT'S HIGHLY RECOMMENDED that these MO's have been
c                                                degeneracy-based diagonalized.
      double precision erg0(nb,ns)    ! [input], unperturbed energies of all MO's.
      integer nch
      integer orbid(nch)
      integer sid(nch)
      double precision erg(nch)       ! [input], corrected energies in Hartree.
      integer g_mohi(ns)              ! [input], interaction Hamiltonian matrix represented in unperturbed (orthonormal) MO basis.
      integer ord                     ! [input], maximum correction order.
      integer g_cmo                   ! [output], coefficients of corrected orbitals. The size should be [nb]x[nch*(ord+1)]
      
      integer g_ms, g_it, g_tmp(2)
      integer i, ich, iv, is, icmo


      call ga_sync()

      
      if (ord < 0)
     &     call errquit('WFAT: Argument ''ord'' of '//
     &     '''wfat_bw_calcwfn'' must be zero or positive integer',
     &     1, UERR)
      

      if (.not. ga_create(mt_dbl, nb, 1, 'mask', 0, 0, g_ms))
     &     call errquit('WFAT: Cannot create mask', 1, GA_ERR)
      call ga_zero(g_ms)

      if (.not. ga_duplicate(g_mohi(1), g_it, 'itermat'))
     &     call errquit('WFAT: Cannot create itermat', 0, GA_ERR)

      do i = 1, 2
         if (.not. ga_duplicate(g_it, g_tmp(i), 'temp'))
     &        call errquit('WFAT: Cannot create temp', 0, GA_ERR)
      enddo
      

      do ich = 1, nch
         iv = orbid(ich)
         is = sid(ich)
         call ga_copy(g_mohi(is), g_it)
         call ga_zero(g_tmp(1))
         call ga_zero(g_tmp(2))

         icmo = (ich-1)*(ord+1) + 1
         call wfat_bw_calcwfn0(nb, iv, g_mo(is), erg(ich), erg0(:,is),
     &        ord, g_ms, g_it, g_tmp, icmo, g_cmo)
      enddo


      if (.not. ga_destroy(g_ms))
     &     call errquit('WFAT: Cannot destroy mask', 1, GA_ERR)
      if (.not. ga_destroy(g_it))
     &     call errquit('WFAT: Cannot destroy itermat', 1, GA_ERR)
      do i = 1, 2
         if (.not. ga_destroy(g_tmp(i)))
     &        call errquit('WFAT: Cannot destroy temp', 0, GA_ERR)
      enddo

      call ga_sync()

      end
c===============================================================
      

c===============================================================
c===============================================================
      subroutine wfat_bw_calcwfn0(nb, iv, g_v, erg, erg0, ord, g_ms,
     &           g_it, g_tmp, icv, g_cv)

      implicit none
#include "errquit.fh"      
#include "global.fh"
      integer nb
      integer iv
      integer g_v                ! [input], all unperturbed MO's. 
      double precision erg       ! [input], corrected energy in Hartree of the iv-th orbital.
      double precision erg0(nb)  ! [input], unperturbed energies of all MO's having spin 'is'.
      integer ord
      integer g_ms               ! [output], mask vector. A scratch array having the size of [nb]x[1].
      integer g_it               ! [in/output], on input, it must contain the matrix of interaction Hamiltonian
c                                               represented in unperturbed (orthonormal) MO basis. A scratch array having the size of [nb]x[nb].
      integer g_tmp(2)           ! [output], scratch array, each having size of [nb]x[nb]
      integer icv                ! [input], column of g_cv in which the 0-th order coefficient will be stored. (icv+1)-th to (icv+ord)-th columns will contain 1st, ..., ord-th correction orders.
      integer g_cv               ! [output], corrections to orbital with MO index iv and spin index is. The size is [nb]x[ord+1].

      integer i, j, ndid
      double precision msk


      call ga_sync()
      ndid = ga_nodeid()

      if (ord < 0)
     &     call errquit('WFAT: Argument ''ord'' of '//
     &     '''wfat_bw_calcwfn'' must be zero or positive integer',
     &     1, UERR)
      
      
c     ==== Construct mask vector ====
      call ga_zero(g_ms)
      do i = 1, nb     ! Loop over rows
         if (i == iv) then
            msk = 0.0d0
         else
            msk = 1.0d0 / (erg - erg0(i))
         endif
         if (ndid==0) call ga_put(g_ms, i, i, 1, 1, msk, 1)
      enddo
      

c     ==== Construct iterator matrix ====
      do i = 1, nb
         call ga_elem_multiply_patch(
     &        g_ms, (/ 1, 1 /), (/ nb, 1 /),
     &        g_it, (/ 1, i /), (/ nb, i /),
     &        g_it, (/ 1, i /), (/ nb, i /))
      enddo


c     ==== Get the 0th order wave function ====
      call ga_copy_patch('n', g_v, 1, nb, iv, iv, g_cv, 1, nb, icv, icv)


c     ==== Calculate higher order corrections ====
      do i = 1, ord
         if (i == 1) then
            call ga_copy_patch('n',
     &                         g_it,     1, nb, 1, nb,
     &                         g_tmp(1), 1, nb, 1, nb)
         else
            call ga_matmul_patch('n', 'n', 1.0d0, 0.0d0,
     &                           g_tmp(1), 1, nb, 1, nb,
     &                           g_it,     1, nb, 1, nb,
     &                           g_tmp(2), 1, nb, 1, nb)
            call ga_copy_patch('n',
     &                         g_tmp(2), 1, nb, 1, nb,
     &                         g_tmp(1), 1, nb, 1, nb)
         endif

         j = icv + i
         call ga_matmul_patch('n', 'n', 1.0d0, 0.0d0,
     &        g_v,      1, nb, 1,  nb,
     &        g_tmp(1), 1, nb, iv, iv,
     &        g_cv,     1, nb, j,  j)
         
      enddo

      call ga_sync()
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_perturb_analconv_bw(nb, ns, erg0, nch, sid, erg,
     &           numitr, numord, max_c, max_c_i)

      implicit none
#include "stdio.fh"
#include "global.fh"      
      integer nb
      integer ns
      double precision erg0(nb,ns)
      integer nch
      integer sid(nch)
      double precision erg(nch)
      integer numitr(nch)
      integer numord(nch)
      double precision max_c(nch)
      integer max_c_i(2, nch)

      character*15 lb
      character*4 r_s, c_s
      integer i, is, titpos, lwidth
      double precision derg


      if (ga_nodeid() /= 0) return

c     ==== Print title ====
      titpos = 5 + 124/2
      write(luout, '(2/)')
      call util_print_centered(luout, 'Convergence properties in '//
     &     'Brillouin-Wigner perturbation', titpos, .true.)
      write(luout, '(" ")')

      
      lwidth = 124
      call wfat_print_repeat(luout, lwidth, 5, '-')
      write(luout, '(1x)')
      write(luout, '(5x, a4)', advance='no') 'Ch.'
      call util_flush(luout)
      write(luout, '(2x, a25)', advance='no') 'Max. |<m1|H_int|m2>/dE|'
      call util_flush(luout)
      write(luout, '(2x, a14)', advance='no') 'dE = e_n-e_m1'
      call util_flush(luout)
      write(luout, '(2x, a15)', advance='no') '( m1 , m2 )'
      call util_flush(luout)
      write(luout, '(2x, a14)', advance='no') 'e_n'
      call util_flush(luout)
      write(luout, '(2x, a18)', advance='no') '#iter to solve e_n'
      call util_flush(luout)
      write(luout, '(2x, a22)') '#order to converge e_n'
      write(luout, '(5x, 102x, a22)') 'at the last iter.'
      call wfat_print_repeat(luout, lwidth, 5, '-')
      write(luout, '(1x)')
      

      do i = 1, nch
         is = sid(i)
         write(r_s, '(i4)') max_c_i(1, i)
         write(c_s, '(i4)') max_c_i(2, i)
         lb = '( '//trim(adjustl(r_s))//' , '//trim(adjustl(c_s))// ' )'
         derg = erg(i) - erg0( max_c_i(1, i) , is )
         write(luout, '(5x, i4)', advance='no') i
         call util_flush(luout)
         write(luout, '(2x, es25.6)', advance='no') max_c(i)
         call util_flush(luout)
         write(luout, '(2x, es14.6)', advance='no') derg
         call util_flush(luout)
         write(luout, '(2x, a15)', advance='no') trim(adjustl(lb))
         call util_flush(luout)
         write(luout, '(2x, es14.6)', advance='no') erg(i)
         call util_flush(luout)
         write(luout, '(2x, i18)', advance='no') numitr(i)
         call util_flush(luout)
         write(luout, '(2x, i22)') numord(i)
      enddo
      call wfat_print_repeat(luout, lwidth, 5, '-')
      write(luout, '(1x)')
      
      
      end
c===============================================================
      
