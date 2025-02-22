c===============================================================
c===============================================================
      subroutine wfat_perturb_rs(geom, bas, nb, ns, g_mo, erg0, g_hi,
     &     nch, orbid, sid, dgnthr, offthr, max_c, max_c_i, erg, g_cmo)

      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer geom
      integer bas
      integer nb
      integer ns
      integer g_mo(ns)
      double precision erg0(nb,ns)
      integer g_hi
      integer nch
      integer orbid(nch)
      integer sid(nch)
      double precision dgnthr
      double precision offthr
      double precision max_c(nch)     ! [output], maximum value of |<m|H_int|n>/(e_n-e_m)>| for each n in orbid array.                          
      integer max_c_i(nch)            ! [output], the orbital m corresponding to max_cc.                                                        
      double precision erg(0:1,nch)   ! [output], 0th and 1st orders correction of the energy of MO in g_cv. The size should be [nb]x[nch*2].   
      integer g_cmo                   ! [output], 0th and 1st orders correction of MO coefficient's with the MO ID's indicated by orbid and sid.

      character*9 movname(2)
      character*7 hiname(2)
      integer i, ich, is, ndid, ia
      integer g_mo2(ns), g_mohi(ns)

cOLD      parameter (dgnthr=1.0d-6, offthr=1.0d-10)


      call ga_sync()
      ndid = ga_nodeid()
      movname = (/ 'mov_alpha', 'mov_beta' /)      
      hiname = (/ 'h_int_a', 'h_int_b' /)
      

c     ==== Generate new MO's based on degeneracy ====
      do is = 1, ns
         if (.not. ga_duplicate(g_mo(is), g_mo2(is), movname(is)))
     &        call errquit('WFAT: Cannot duplicate '//movname(is)//'.',
     &        1, GA_ERR)
         call ga_copy(g_mo(is), g_mo2(is))
      enddo      
      call wfat_diag_perturb(nb, ns, g_mo2, erg0, g_hi, dgnthr, offthr)


c     ==== Calculate H_int matrix in the new MO basis ====
      do is = 1, ns
cforcheck         if (ndid==0)
cforcheck     &        write(luout, '("Adapted ", a, " coefficients :")')
cforcheck     &        trim(adjustl(movname(is)))    
cforcheck         call ga_print(g_mo2(is))
         
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

      
c     ==== Calculate 1st order correction to ionized orbitals ====      
      call wfat_rs_ergwfn(bas, nb, ns, g_mo2, erg0, nch, orbid, sid, 
     &                    g_mohi, dgnthr, max_c, max_c_i, erg, g_cmo)


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
      subroutine wfat_rs_ergwfn(bas, nb, ns, g_v, erg0, 
     &           nch, orbid, sid, g_mohi, dgnthr, max_cc, max_cc_i,
     &           erg, g_cv)  

      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer bas                       ! [input], basis handle in which g_v is represented.
      integer nb                        ! [input], no. of basis functions in bas.
      integer ns                        ! [input], no. of spin channels.
      integer g_v(ns)                   ! [input], all MO's after degeneracy-based mixing.
      double precision erg0(nb,ns)      ! [input], unperturbed energy associated with MO's in g_v.
      integer nch                       ! [input], number of MO to be corrected.
      integer orbid(nch)                ! [input], indices of MO's in g_v whose correction is sought.
      integer sid(nch)                  ! [input], spin values of the MO's in g_v whose correction is sought.
      integer g_mohi(ns)                ! [input], interaction Hamiltonian.
      double precision dgnthr           ! [input], threshold criteria for degeneracies.
      double precision max_cc(nch)      ! [output], maximum value of |<m|H_int|n>/(e_n-e_m)>| for each n in orbid array.
      integer max_cc_i(nch)             ! [output], the orbital m corresponding to max_cc.
      double precision erg(0:1,nch)     ! [output], 0th and 1st orders correction of the energy of MO in g_cv. The size should be [nb]x[nch*2].
      integer g_cv                      ! [output], 0th and 1st orders correction of MO coefficient's with the MO ID's indicated by orbid and sid.

      character*9 movname(2)
      integer i, ii, iv, is, ich, j, ndid
      integer l_sc, k_sc
      double precision cc, delta_erg, hithr


      call ga_sync()
      ndid = ga_nodeid()
      hithr = abs(dgnthr) * 1.0d-4


c     ==== Get 0th order MO's and energies ====
      do i = 1, nch
         iv = orbid(i)
         is = sid(i)
         j = (i-1)*2 + 1
         call ga_copy_patch('n',
     &                      g_v(is), 1, nb, iv, iv,
     &                      g_cv,    1, nb, j,  j)
         erg(0,i) = erg0(iv,is)
      enddo
      
      
      if (.not. ma_alloc_get(mt_dbl, nb, 'ma_scratch', l_sc, k_sc))
     &     call errquit('WFAT: Cannot allocate scratch', 1, MA_ERR)
      
c     ==== Loop over user-requested orbitals ====
      do ich = 1, nch
         is = sid(ich)
         i = orbid(ich)

         call ga_get(g_mohi(is), 1, nb, i, i, dbl_mb(k_sc), nb)
         erg(1,ich) = dbl_mb(k_sc+i-1)
         ii = (ich-1)*2 + 2
         call nga_zero_patch(g_cv, (/ 1, ii /), (/ nb, ii /))
         max_cc(ich) = 0.0d0
         max_cc_i(ich) = -1000

         
c        ==== Loop over all MO's ====
         do j = 1, nb

c           ==== Calculate coefficient of psi^(1) in MO representation ====
            delta_erg = erg0(i,is) - erg0(j,is)
            if (abs(delta_erg) < dgnthr .and. 
     &           abs(dbl_mb(k_sc+j-1)) < hithr) then
               cc = 0.0d0
            else if (abs(delta_erg) < dgnthr .and. 
     &              abs(dbl_mb(k_sc+j-1)) >= hithr) then
               if (i == j) then
                  cc = 0.0d0
               else
                  call errquit('WFAT: MO have not been adapted for '//
     &                 'degeneracies before calling wfat_perturb_1st.',
     &                 1, UERR)
               endif
            else
               cc = dbl_mb(k_sc+j-1) / delta_erg
            endif            

            if (abs(cc) > max_cc(ich)) then
               max_cc(ich) = abs(cc)
               max_cc_i(ich) = j
            endif

c           ==== Calculate psi^(1) in basis set representation ====
            ii = (ich-1)*2 + 2
            call ga_add_patch(
     &           1.0d0, g_cv,    1, nb, ii, ii,
     &           cc,    g_v(is), 1, nb, j,  j,
     &           g_cv, 1, nb, ii, ii)
         enddo
      enddo


      if (.not. ma_free_heap(l_sc))
     &     call errquit('WFAT: Cannot deallocate scratch.', 1, MA_ERR)

      call ga_sync()
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_perturb_analconv_rs(nb, ns, erg0, nch, orbid, sid,
     &           max_c, max_c_i)

      implicit none
#include "stdio.fh"
#include "global.fh"
      integer nb
      integer ns
      double precision erg0(nb,ns)
      integer nch
      integer orbid(nch)
      integer sid(nch)
      double precision max_c(nch)
      integer max_c_i(nch)

      integer i, iv, is, titpos, lwidth
      double precision derg


      if (ga_nodeid() /= 0) return

c     ==== Print title ====
      titpos = 5 + 53/2
      write(luout, '(2/)')
      call util_print_centered(luout, 'Convergence properties in '//
     &     'Rayleigh-Schroedinger perturbation', titpos, .true.)
      write(luout, '(" ")')

      
      lwidth = 53
      call wfat_print_repeat(luout, lwidth, 5, '-')
      write(luout, '(1x)')
      write(luout, '(5x, a4, 3x, a25, 1x, a14, 1x, a5)') 
     &     'Ch.', 'Max. | <m|H_int|n>/dE |', 'dE = e_n-e_m', 'm'
      call wfat_print_repeat(luout, lwidth, 5, '-')
      write(luout, '(1x)')
      

      do i = 1, nch
         iv = orbid(i)
         is = sid(i)
         derg = erg0(iv,is) - erg0( max_c_i(i) , is )
         write(luout, '(5x, i4, 3x, es25.6, 1x, es14.6, 1x, i5)')
     &        i, max_c(i), derg, max_c_i(i)
      enddo
      call wfat_print_repeat(luout, lwidth, 5, '-')
      write(luout, '(1x)')
      
      end
c===============================================================
 
