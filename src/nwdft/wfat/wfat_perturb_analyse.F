c===============================================================
c===============================================================
      subroutine wfat_perturb_analyse(rtdb, geom, bas, nb, ns, ne, nch, 
     &           sid, ord, maxord, g_hi, g_mo0, g_cmo, erg, tol2e, 
     &           oskel, ixcfit, iolgc, rmode)

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"      
      integer rtdb
      integer geom
      integer bas
      integer nb
      integer ns
      integer ne(ns)
      integer nch
      integer sid(nch)
      integer ord
      integer maxord
      integer g_hi
      integer g_mo0(ns)
      integer g_cmo
      double precision erg(0:maxord, nch)
      double precision tol2e
      logical oskel
      integer ixcfit
      logical iolgc
      logical rmode

      logical integxy, print_veff
      integer i, j, k
      integer g_s, g_tmp1, g_cv, l_cv, k_cv, l_vl, k_vl
      double precision norm, vpt(3)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      

      call ga_sync()

      
c     ==== Calculate overlap matrix ====
      g_s = ga_create_atom_blocked(geom, bas, 'overlapmat')
      call ga_zero(g_s)
      call int_1e_ga(bas, bas, g_s, 'overlap', oskel)
      
      
c     ==== Print energies ====
      call wfat_perturb_analerg(maxord, nch, erg)

      
c     ==== Print norms ====
      call wfat_perturb_analnorm(rtdb, geom, bas, nb, ns, ne, nch, sid, 
     &                           ord, g_hi, g_mo0, g_cmo, g_s, tol2e, 
     &                           oskel, ixcfit, iolgc)

      
c     ==== Print Mulliken populations ====
      call wfat_perturb_analmulliken(geom, bas, nb, nch, ord, g_cmo,
     &     g_s)


c     ==== Print overlaps ====
      call wfat_perturb_analovl(geom, bas, nb, nch, ord, g_cmo, g_s)
      

c     ==== Print xy-integrated density ====
      if (.not. rtdb_get(rtdb, 'wfat:perturb_dens_integxy', mt_log,
     &     1, integxy))
     &     integxy = .false.
      if (rmode .and. integxy) then
c        ==== Renormalize corrected ionized orbitals ====         
         if (.not. ga_create(mt_dbl, nb, nch, 'corct', 0, 0, g_cv))
     &        call errquit('WFAT: Cannot create corct', 1, GA_ERR)
         if (.not. ga_create(mt_dbl, nb, 1, 'tmp1', 0, 0, g_tmp1))
     &        call errquit('WFAT: Cannot create tmp1', 1, GA_ERR)
         do i = 1, nch
            call ga_zero(g_tmp1)
            do j = 0, ord
               k = (i-1)*(ord+1) + j + 1
               call ga_add_patch(
     &              1.0d0, g_tmp1, 1, nb, 1, 1,
     &              1.0d0, g_cmo,  1, nb, k, k,
     &              g_tmp1, 1, nb, 1, 1)
            enddo
            call wfat_ga_mat_element(g_tmp1, g_tmp1, g_s, norm)
            call ga_scale(g_tmp1, 1.0d0/sqrt(norm))
            call ga_copy_patch('n', g_tmp1, 1, nb, 1, 1,
     &                        g_cv, 1, nb, i, i)
         enddo


         

cforcheck         vpt = (/ 0.0d0, -0.0d0, 3.0d0 /)
cforcheck         if (.not. ma_alloc_get(mt_dbl, nch, 'val', l_vl, k_vl))
cforcheck     &        call errquit('WFAT: Cannot allocate val', 1, MA_ERR)
cforcheck         call wfat_eval_mo0_1pt(geom, bas, vpt, nch, g_cv, dbl_mb(k_vl))
cforcheck         if (ga_nodeid()==0) then
cforcheck            write(luout, '("Values of corrected MO''s at ", 3es14.4)')
cforcheck     &           vpt
cforcheck            do i = 1, nch            
cforcheck               write(luout, '(3x, i4, es14.6)') i, dbl_mb(k_vl+i-1)
cforcheck            enddo
cforcheck         endif
cforcheck         if (.not. ma_free_heap(l_vl))
cforcheck     &        call errquit('WFAT: Cannot deallocate val', 1, MA_ERR)


         

         if (.not. ma_alloc_get(mt_dbl, nb*nch, 'corct_MO', l_cv, k_cv))
     &        call errquit('WFAT: Cannot allocate corct_MO', 1, MA_ERR)
         call ga_get(g_cv, 1, nb, 1, nch, dbl_mb(k_cv), nb)
         call wfat_print_dens_integxy(rtdb, nch, dbl_mb(k_cv))

         if (.not. ga_destroy(g_cv))
     &        call errquit('WFAT: Cannot destroy corct', 1, GA_ERR)
         if (.not. ga_destroy(g_tmp1))
     &        call errquit('WFAT: Cannot destroy tmp1', 1, GA_ERR)
         if (.not. ma_free_heap(l_cv))
     &        call errquit('WFAT: Cannot deallocate corct', 1, MA_ERR)    
      endif


c     ==== Coefficients ====
cforcheck      do ich = 1, nch
c        ==== Print coefficients ====
cforcheck         if (ndid==0) then
cforcheck            write(luout, '(5x, "Coefficients of 1st order correction")')
cforcheck            call wfat_dprintmat(1, nb, 'e', dbl_mb(i_v1))
cforcheck            write(luout, '(5x, "Coefficients of normalized corrected")')
cforcheck            call wfat_dprintmat(1, nb, 'e', dbl_mb(i_v))
cforcheck         endif
cforcheck      enddo


c     ==== Print effective potential ====
      if (.not. rtdb_get(rtdb, 'wfat:perturb_veff_print', mt_log,
     &     1, print_veff))
     &     print_veff = .false.
      if (rmode .and. print_veff) then
         call wfat_perturb_potprint(rtdb, geom, bas, nb, ns, ne, 
     &        g_mo0, tol2e, oskel, ixcfit, iolgc, g_hi)
      endif
      

c     ==== Deallocations ====
      if (.not. ga_destroy(g_s))
     &     call errquit('WFAT: Cannot destroy overlapmat', 1, GA_ERR)

      call ga_sync()

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_perturb_analerg(ord, nch, erg)

      implicit none
#include "stdio.fh"
#include "global.fh"
      integer ord
      integer nch
      double precision erg(0:ord, nch)

      character*5 ic_s
      character*2 suf
      character*7 ord_s
      integer maxcol, nblock, titpos, lwidth
      integer ic1, ic2, ic, ib, ir

      parameter (maxcol=5)


      if (ga_nodeid() /= 0) return


c     ==== Print title ====
      titpos = 5 + (10 + maxcol*15)/2
      write(luout, '(2/)')
      call util_print_centered(luout, 'Perturbed Energies',
     &     titpos, .true.)
      write(luout, '(" ")')      

      if (mod(nch,maxcol) == 0) then
         nblock = nch/maxcol
      else
         nblock = nch/maxcol + 1
      endif


c     ==== Print tables ====
      ic1 = 1
      ic2 = maxcol
      do ib = 1, nblock
         lwidth = 7 + (min(ic2,nch) - ic1 + 1)*16
         call wfat_print_repeat(luout, lwidth, 5, '-')
         write(luout, '(1x)')

c        ==== Print column titles ====
         write(luout, '(5x, a7)', advance='no') 'Order'
         call util_flush(luout)
         do ic = ic1, ic2
            write(ic_s, '(i5)') ic
            write(luout, '(2x, a14)', advance='no')
     &           'Ch. '//trim(adjustl(ic_s))
            call util_flush(luout)
            if (ic == nch) exit
         enddo
         write(luout, '(" ")')

         call wfat_print_repeat(luout, lwidth, 5, '-')
         write(luout, '(1x)')

c        ==== Print energy of each order ====
         do ir = 0, ord
            call wfat_ordinal_suffix(ir, suf)
            write(luout, '(5x, i4, a3)', advance='no') ir, '-'//suf
            call util_flush(luout)
    
            do ic = ic1, ic2
               write(luout, '(2x, es14.6)', advance='no') erg(ir,ic)
               call util_flush(luout)
               if (ic == nch) exit
            enddo
            write(luout, '(" ")')
         enddo

         call wfat_print_repeat(luout, lwidth, 5, '-')
         write(luout, '(1x)')

         ic1 = ic1 + maxcol
         ic2 = ic2 + maxcol
         write(luout, '(1x)')
         
      enddo

      write(luout, '(5x, "Note :")')
      write(ord_s, '(i7)') ord      
      write(luout, 51) trim(adjustl(ord_s))
 51   format(7x,
     &     'The correction orders of energy shown above (= ', a, ') ',
     &     'is always larger than or equal ', /, 
     &     7x,
     &     'to the actual order used to compute the corrected wave ',
     &     'functions, which may be found in ', /,
     &     7x,
     &     '''Overlap with the unperturbed orbital'' section further ', 
     &     'below. In Brillouin-Wigner perturbation, ', /,
     &     7x,
     &     'the maximum order for energy is used to solve the ',
     &     'implicit equation for energy iteratively.')
      
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_perturb_analnorm(rtdb, geom, bas, nb, ns, ne, nch, 
     &                                 sid, ord, g_hi, g_mo0, g_cmo, 
     &                                 g_s, tol2e, oskel, ixcfit, iolgc)

      implicit none
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      integer geom
      integer bas
      integer nb
      integer ns
      integer ne(ns)
      integer nch
      integer sid(nch)
      integer ord
      integer g_hi
      integer g_mo0(ns)
      integer g_cmo
      integer g_s
      double precision tol2e
      logical oskel
      integer ixcfit
      logical iolgc

      character*5 ic_s
      character*2 suf
      integer maxcol, titpos, nblock, ndid, lwidth
      integer i, j, k, ib, ic, ic1, ic2, ir, is
      integer l_n, k_n, l_tn, k_tn, l_ev, k_ev, i_n, g_tmp, g_fk(ns)

      parameter (maxcol=5)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked


      ndid = ga_nodeid()


c     ==== Calculate norm of each orders ====
      if (.not. ga_create(mt_dbl, nb, 1, 'tmp', 0, 0, g_tmp))
     &     call errquit('WFAT: Cannot create tmp', 1, GA_ERR)
      if (.not. ma_alloc_get(mt_dbl, (ord+1)*nch, 'norms', l_n, k_n))
     &     call errquit('WFAT: Cannot allocate norms', 1, MA_ERR)
      k = 0
      do i = 1, nch
         do j = 0, ord
            k = k + 1
            call ga_copy_patch('n', g_cmo, 1, nb, k, k,
     &                         g_tmp, 1, nb, 1, 1)
            call wfat_ga_mat_element(g_tmp, g_tmp, g_s, dbl_mb(k_n+k-1))
         enddo
      enddo


      
c     ==== Calculate expectation value of perturbed Fock ====
      if (ndid==0) write(luout, '("here t0, geom = ", i4)') geom
      do is = 1, ns
         g_fk(is) = ga_create_atom_blocked(geom, bas, 'fock')
         call ga_zero(g_fk(is))
      enddo
cseelater      call wfat_get_fock(rtdb, geom, bas, nb, ns, ne, g_mo0,
cseelater     &                   tol2e, oskel, ixcfit, iolgc, g_fk)   chg, mult, and doxc have been added as new entries of this subroutine, update this call.
cseelater      if (ndid==0) write(luout, '("here t3")')
      do is = 1, ns
         call ga_dadd(1.0d0, g_fk(is), 1.0d0, g_hi, g_fk(is))
      enddo
      
      

c     ==== Calculate norm of total corrected MO's ====
      if (.not. ma_alloc_get(mt_dbl, nch, 'tot_norms', l_tn, k_tn))
     &     call errquit('WFAT: Cannot allocate tot_norms', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nch, 'expval', l_ev, k_ev))
     &     call errquit('WFAT: Cannot allocate expval', 1, MA_ERR)
      do i = 1, nch
         call ga_zero(g_tmp)
         do j = 0, ord
            k = (i-1)*(ord+1) + j + 1
            call ga_add_patch(
     &           1.0d0, g_tmp, 1, nb, 1, 1,
     &           1.0d0, g_cmo, 1, nb, k, k,
     &           g_tmp, 1, nb, 1, 1)
         enddo
         call wfat_ga_mat_element(g_tmp, g_tmp, g_s, dbl_mb(k_tn+i-1))
         is = sid(i)
         call wfat_ga_mat_element(g_tmp, g_tmp, g_fk(is),
     &        dbl_mb(k_ev+i-1))
      enddo
      if (.not. ga_destroy(g_tmp))
     &     call errquit('WFAT: Cannot destroy tmp', 1, MA_ERR)
      do is = 1, ns
         if (.not. ga_destroy(g_fk(is)))
     &        call errquit('WFAT: Cannot destroy fock', 1, MA_ERR)
      enddo
      

c     ==== Print title ====
      titpos = 5 + (10 + maxcol*15)/2
      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Norm of perturbation '//
     &        'corrections', titpos, .true.)
         write(luout, '(" ")')
      endif

      if (mod(nch,maxcol) == 0) then
         nblock = nch/maxcol
      else
         nblock = nch/maxcol + 1
      endif


c     ==== Print tables ====
      if (ndid == 0) then
         ic1 = 1
         ic2 = maxcol
         do ib = 1, nblock
            lwidth = 7 + (min(ic2,nch) - ic1 + 1)*16
            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
            
c           ==== Print column titles ====
            write(luout, '(5x, a7)', advance='no') 'Order'
            call util_flush(luout)
            do ic = ic1, ic2
               write(ic_s, '(i5)') ic
               write(luout, '(2x, a14)', advance='no')
     &              'Ch. '//trim(adjustl(ic_s))
               call util_flush(luout)
               if (ic == nch) exit
            enddo
            write(luout, '(" ")')

            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
            
c           ==== Print norm of each orders ====
            do ir = 0, ord
               call wfat_ordinal_suffix(ir, suf)
               write(luout, '(5x, i4, a3)', advance='no') ir,
     &              '-'//suf
               call util_flush(luout)
         
               do ic = ic1, ic2
                  i_n = k_n + (ic-1)*(ord+1) + ir
                  write(luout, '(2x, es14.6)', advance='no') dbl_mb(i_n)
                  call util_flush(luout)
                  if (ic == nch) exit
               enddo
               write(luout, '(" ")')
            enddo
         
c           ==== Print norm of corrected MO's ====
            write(luout, '(5x, a7)', advance='no') 'Total'
            call util_flush(luout)
            do ic = ic1, ic2
               write(luout, '(2x, es14.6)', advance='no')
     &              dbl_mb(k_tn+ic-1)
               call util_flush(luout)
               if (ic == nch) exit
            enddo
            write(luout, '(" ")')


c           ==== Print exp. value of total Fock ====
            write(luout, '(5x, a7)', advance='no') 'Expval.'
            call util_flush(luout)
            do ic = ic1, ic2
               write(luout, '(2x, es14.6)', advance='no')
     &              dbl_mb(k_ev+ic-1)/dbl_mb(k_tn+ic-1)
               call util_flush(luout)
               if (ic == nch) exit
            enddo
            write(luout, '(" ")')
            

            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
            
c           ==== Increment ====
            ic1 = ic1 + maxcol
            ic2 = ic2 + maxcol
            write(luout, '(/)')
         enddo

      endif


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_n))
     &     call errquit('WFAT: Cannot deallocate norms', 1, MA_ERR)
      if (.not. ma_free_heap(l_tn))
     &     call errquit('WFAT: Cannot deallocate tot_norms', 1, MA_ERR)
      if (.not. ma_free_heap(l_ev))
     &     call errquit('WFAT: Cannot deallocate expval', 1, MA_ERR)
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_perturb_analovl(geom, bas, nb, nch, ord, g_cmo,
     &           g_s)

      implicit none
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer geom
      integer bas
      integer nb
      integer nch
      integer ord
      integer g_cmo
      integer g_s

      character*5 ic_s
      character*2 suf
      integer maxcol, titpos, nblock, ndid, lwidth
      integer i, j, k, ii, ib, ic, ic1, ic2, ir
      integer l_ov, k_ov, i_ov, l_ovt, k_ovt, g_v0, g_cv, g_tot

      parameter (maxcol=5)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked


      ndid = ga_nodeid()
      

c     ==== Calculate overlap of each orders with 0th ====
      if (.not. ga_create(mt_dbl, nb, 1, '0th', 0, 0, g_v0))
     &     call errquit('WFAT: Cannot create 0th', 1, GA_ERR)
      if (.not. ga_duplicate(g_v0, g_cv, 'corct'))
     &     call errquit('WFAT: Cannot create corct', 1, GA_ERR)
      if (.not. ga_duplicate(g_v0, g_tot, 'tot'))
     &     call errquit('WFAT: Cannot create tot', 1, GA_ERR)
      
      if (.not. ma_alloc_get(mt_dbl, (ord+1)*nch, 'overlap',
     &     l_ov, k_ov))
     &     call errquit('WFAT: Cannot allocate overlap', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nch, 'overlap_tot',
     &     l_ovt, k_ovt))
     &     call errquit('WFAT: Cannot allocate overlap_tot', 1, MA_ERR)
      k = 0
      do i = 1, nch
         ii = (i-1)*(ord+1) + 1
         call ga_copy_patch('n', g_cmo, 1, nb, ii, ii,
     &                      g_v0, 1, nb, 1, 1)
         call ga_zero(g_tot)
         do j = 0, ord
            k = k + 1
            call ga_copy_patch('n', g_cmo, 1, nb, k, k,
     &                         g_cv, 1, nb, 1, 1)
            call wfat_ga_mat_element(g_cv, g_v0, g_s, dbl_mb(k_ov+k-1))

            call ga_add_patch(
     &           1.0d0, g_tot, 1, nb, 1, 1,
     &           1.0d0, g_cmo, 1, nb, k, k,
     &           g_tot, 1, nb, 1, 1)
         enddo
         call wfat_ga_mat_element(g_tot, g_v0, g_s, dbl_mb(k_ovt+i-1))
      enddo      
      if (.not. ga_destroy(g_v0))
     &     call errquit('WFAT: Cannot destroy 0th', 1, GA_ERR)
      if (.not. ga_destroy(g_cv))
     &     call errquit('WFAT: Cannot destroy corct', 1, GA_ERR)
      if (.not. ga_destroy(g_tot))
     &     call errquit('WFAT: Cannot destroy tot', 1, GA_ERR)

                  
c     ==== Print title ====
      titpos = 5 + (10 + maxcol*15)/2
      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Overlap with the '//
     &        'unperturbed orbital', titpos, .true.)
         write(luout, '(" ")')
      endif

      if (mod(nch,maxcol) == 0) then
         nblock = nch/maxcol
      else
         nblock = nch/maxcol + 1
      endif


c     ==== Print tables ====
      if (ndid == 0) then
         ic1 = 1
         ic2 = maxcol
         do ib = 1, nblock
            lwidth = 7 + (min(ic2,nch) - ic1 + 1)*16
            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
         
c           ==== Print column titles ====
            write(luout, '(5x, a7)', advance='no') 'Order'
            call util_flush(luout)
            do ic = ic1, ic2
               write(ic_s, '(i5)') ic
               write(luout, '(2x, a14)', advance='no')
     &              'Ch. '//trim(adjustl(ic_s))
               call util_flush(luout)
               if (ic == nch) exit
            enddo
            write(luout, '(" ")')

            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
            
c           ==== Print overlap from each orders ====
            do ir = 0, ord
               call wfat_ordinal_suffix(ir, suf)
               write(luout, '(5x, i4, a3)', advance='no') ir,
     &              '-'//suf
               call util_flush(luout)
         
               do ic = ic1, ic2
                  i_ov = k_ov + (ic-1)*(ord+1) + ir
                  write(luout, '(2x, es14.6)', advance='no')
     &                 dbl_mb(i_ov)
                  call util_flush(luout)
                  if (ic == nch) exit
               enddo
               write(luout, '(" ")')
            enddo

c           ==== Print norm of corrected MO's ====
            write(luout, '(5x, a7)', advance='no') 'Total'
            call util_flush(luout)
            do ic = ic1, ic2
               write(luout, '(2x, es14.6)', advance='no')
     &              dbl_mb(k_ovt+ic-1)
               call util_flush(luout)
               if (ic == nch) exit
            enddo
            write(luout, '(" ")')

            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
            
c           ==== Increment ====
            ic1 = ic1 + maxcol
            ic2 = ic2 + maxcol
            write(luout, '(x)')
         enddo

         write(luout, '(5x, "Note :")')
         write(luout, 51)
      endif

 51   format(7x,
     &     'The way the standard perturbation method is formulated ',
     &     'is to require that the zeroth ', /
     &     7x, 
     &     'order wavefunction be orthogonal with all higher order ',
     &     'corrections. So, the values in ', /
     &     7x, 
     &     'the rows between the second and the last ones in table ',
     &     'above should be close to zero. ',/
     &     7x, 
     &     'A non-zero value indicates either inappropriate inputs ',
     &     'or bugs in the code.')


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_ov))
     &     call errquit('WFAT: Cannot deallocate overlap', 1, MA_ERR)
      if (.not. ma_free_heap(l_ovt))
     &     call errquit('WFAT: Cannot deallocate overlap_tot',
     &     1, MA_ERR)
      

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_perturb_analmulliken(geom, bas, nb, nch, ord,
     &           g_cmo, g_s)

      implicit none
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
      integer geom
      integer bas
      integer nb
      integer nch
      integer ord
      integer g_cmo
      integer g_s
      
      character*16 aname
      character*5 ic_s
      integer nblock, maxcol, titpos, npo, nat, ndid, lwidth
      integer i, j, k, ich, ic, ic1, ic2, ia, ib
      integer g_cv, g_tmp1, g_tmp2
      integer l_po0, k_po0, l_po, k_po, i_po, i_po0
      double precision norm

      parameter (maxcol=3)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked


      ndid = ga_nodeid()
      

      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 51, GEOM_ERR)      


c     ==== Calculate total, normalized and corrected MO's ====
      if (.not. ga_create(mt_dbl, nb, nch, 'corct', 0, 0, g_cv))
     &     call errquit('WFAT: Cannot create corct', 1, GA_ERR)
      if (.not. ga_create(mt_dbl, nb, 1, 'tmp1', 0, 0, g_tmp1))
     &     call errquit('WFAT: Cannot create tmp1', 1, GA_ERR)
      do i = 1, nch
         call ga_zero(g_tmp1)
         do j = 0, ord
            k = (i-1)*(ord+1) + j + 1
            call ga_add_patch(
     &           1.0d0, g_tmp1, 1, nb, 1, 1,
     &           1.0d0, g_cmo,  1, nb, k, k,
     &           g_tmp1, 1, nb, 1, 1)
         enddo
         call wfat_ga_mat_element(g_tmp1, g_tmp1, g_s, norm)
         call ga_scale(g_tmp1, 1.0d0/sqrt(norm))
         call ga_copy_patch('n', g_tmp1, 1, nb, 1, 1,
     &                     g_cv, 1, nb, i, i)
      enddo
      

c     ==== Calculate Mulliken populations ====
      if (.not. ga_duplicate(g_tmp1, g_tmp2, 'tmp2'))
     &     call errquit('WFAT: Cannot create tmp2', 1, GA_ERR)
      call ga_zero(g_tmp2)
      
      npo = nat * nch
      if (.not. ma_alloc_get(mt_dbl, npo, 'pop0', l_po0, k_po0))
     &     call errquit('WFAT: Cannot allocate pop0', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, npo, 'ovl_pop', l_po, k_po))
     &     call errquit('WFAT: Cannot allocate ovl_pop', 1, MA_ERR)

      do ich = 1, nch
         i = (ich-1)*(ord+1) + 1
         call ga_copy_patch('n', g_cmo, 1, nb, i, i,
     &                      g_tmp1, 1, nb, 1, 1)           ! Unperturbed MO.
         call ga_copy_patch('n', g_cv, 1, nb, ich, ich,
     &                      g_tmp2, 1, nb, 1, 1)           ! Perturbed MO.

         i_po = k_po + (ich-1)*nat
         call wfat_ga_matelm_mulliken(geom, bas, g_tmp1, g_tmp2, g_s,
     &        dbl_mb(i_po))
         i_po0 = k_po0 + (ich-1)*nat
         call wfat_ga_matelm_mulliken(geom, bas, g_tmp1, g_tmp1, g_s,
     &        dbl_mb(i_po0))
      enddo

      if (.not. ga_destroy(g_tmp1))
     &     call errquit('WFAT: Cannot destroy tmp1', 1, GA_ERR)
      if (.not. ga_destroy(g_tmp2))
     &     call errquit('WFAT: Cannot destroy tmp2', 1, GA_ERR)
      if (.not. ga_destroy(g_cv))
     &     call errquit('WFAT: Cannot destroy corct', 1, GA_ERR)
      


c     ==== Print title ====
      titpos = 5 + (7 + maxcol*32)/2
      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Mulliken analysis of the '//
     &        'perturbation', titpos, .true.)
         write(luout, '(" ")')
      endif

      if (mod(nch,maxcol) == 0) then
         nblock = nch/maxcol
      else
         nblock = nch/maxcol + 1
      endif


c     ==== Print tables ====
      if (ndid == 0) then
         ic1 = 1
         ic2 = maxcol
         do ib = 1, nblock
            lwidth = 7 + (min(ic2,nch) - ic1 + 1)*32
            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
         
c           ==== Print column titles part 1 ====
            write(luout, '(5x, a7)', advance='no') 'Atom'
            call util_flush(luout)
            do ic = ic1, ic2
               write(ic_s, '(i5)') ic
               write(luout, '(4x, a28)', advance='no')
     &              'Ch. '//trim(adjustl(ic_s))//'         '
               call util_flush(luout)
               if (ic == nch) exit
            enddo
            write(luout, '(" ")')

c           ==== Print column titles part 2 ====
            write(luout, '(5x, a7)', advance='no') ' '
            call util_flush(luout)
            do ic = ic1, ic2
               write(luout, '(4x, a14, a14)', advance='no')
     &              '<v^(0)|v^(0)>', '<v^(0)|v>'
               call util_flush(luout)
               if (ic == nch) exit
            enddo
            write(luout, '(" ")')

            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
         
c           ==== Print norm of each orders ====
            do ia = 1, nat
               if (.not. geom_cent_tag(geom, ia, aname))
     &              call errquit('WFAT: Getting atomname failed',
     &              52, GEOM_ERR)
               write(luout, '(5x, a7)', advance='no')
     &              trim(aname)
               call util_flush(luout)
         
               do ic = ic1, ic2
                  i_po = k_po + (ic-1)*nat + ia - 1
                  i_po0 = k_po0 + (ic-1)*nat + ia - 1
                  write(luout, '(4x, 2es14.6)', advance='no')
     &                 dbl_mb(i_po0), dbl_mb(i_po)
                  call util_flush(luout)
                  if (ic == nch) exit
               enddo
               write(luout, '(" ")')
            enddo

            call wfat_print_repeat(luout, lwidth, 5, '-')
            write(luout, '(1x)')
                     
c           ==== Increment ====
            ic1 = ic1 + maxcol
            ic2 = ic2 + maxcol
            write(luout, '(x)')
         enddo

         write(luout, '(5x, "Note on the notation :")')
         write(luout, 51)
      endif

 51   format(7x, 
     &     '|v^(i)> is the i-th order correction to |v> such that',
     &     '|v> = |v^(0)> + |v^(1)> + |v^(2)> + ...')
      

      if (.not. ma_free_heap(l_po0))
     &     call errquit('WFAT: Cannot deallocate pop0', 1, MA_ERR)
      if (.not. ma_free_heap(l_po))
     &     call errquit('WFAT: Cannot deallocate pop', 1, MA_ERR)
      

      end
c===============================================================
