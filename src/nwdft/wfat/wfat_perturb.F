c===============================================================
c===============================================================
      subroutine wfat_calc_perturb(rtdb, ef, pr_iony, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom, rmode, nch, orbid, sid)

      implicit none
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
      integer rtdb
      double precision ef  ! [input],
      logical pr_iony      ! [input],
      logical pr_avbeta    ! [input],
      logical pr_avgamma   ! [input],
      logical pr_tmd       ! [input],
      logical pr_tmd_atom  ! [input],
      logical rmode        ! [input],
      integer nch          ! [input],
      integer orbid(nch)   ! [input],
      integer sid(nch)     ! [input],

      integer ndid
      double precision bt, gm


      ndid = ga_nodeid()
      if (ndid == 0) then
         write(luout, '(/)')
         call util_print_centered(
     &        luout, 'Perturbation Theory Analysis', 40, .true.)
      endif
      
      if (ef <= 1.0d-14 .or.
     &     .not. (pr_iony .or. pr_avbeta .or. pr_avgamma .or.
     &     pr_tmd .or. pr_tmd_atom)) then
         if (ndid == 0)
     &        call wfat_perturb_noinfo()
      else
         if (.not. rtdb_get(rtdb, 'wfat:perturb_beta', mt_dbl, 1, bt))
     &        call errquit('WFAT: Cannot read perturb_beta',
     &        1, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'wfat:perturb_gamma', mt_dbl, 1, gm))
     &        call errquit('WFAT: Cannot read perturb_gamma',
     &        1, RTDB_ERR)

         call wfat_perturb_info(rtdb, bt, gm, ef)
         call wfat_dip_perturb(rtdb, bt, gm, ef, nch, orbid, sid, rmode)
      endif

      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_dip_perturb(rtdb, bt, gm, ef, nch, orbid, sid,
     &           rmode)

      implicit none
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "wfat_qm_param.fh"
      integer rtdb
      double precision bt               ! [input], beta Euler angle in degree.
      double precision gm               ! [input], gamma Euler angle in degree.
      double precision ef               ! [input], field strength in a.u. of electric field.
      integer nch
      integer orbid(nch)
      integer sid(nch)
      logical rmode

      character*2 type
      character*8 dipname(3)
      integer i, i_e1, i_e2, dim1, dim2, typ
      integer ndid, ord, maxord, maxitr
      integer g_dp(3), g_hi, g_cmo
      integer l_mx, l_nitr, l_nord, l_e, l_mxi, l_ce
      integer k_mx, k_nitr, k_nord, k_e, k_mxi, k_ce
      double precision dgnthr, offthr, e1thr, itrthr, ordthr, btr, gmr
      double precision ef_vec(3), norm, deg2rad
      parameter (deg2rad=0.01745329251994d0)
      
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      ndid = ga_nodeid()
      
      
      dipname = (/ 'dipole_x', 'dipole_y', 'dipole_z' /)

      
c     ==== Calculate e. field vector in molecular frame ====
      btr = bt * deg2rad
      gmr = gm * deg2rad
      ef_vec(1) = - sin(btr) * cos(gmr) * ef
      ef_vec(2) = sin(btr) * sin(gmr) * ef
      ef_vec(3) = cos(btr) * ef
      
      if (ndid == 0) then
         write(luout, '("ef = ", es14.6)') ef
      endif
      if (ndid==0) write(luout, '("ef_vec = ", 3es14.6)') ef_vec


c     ==== Calculate dipole interaction Hamiltonian ====
      do i = 1, 3
         g_dp(i) = ga_create_atom_blocked(geom_hd, bas_hd, dipname(i))
         call ga_zero(g_dp(i))
      enddo
      call ga_inquire(g_dp(1), typ, dim1, dim2)
      if (ndid==0) write(luout, '("g_dp: dim1, dim2, nbas = ", 3i5)')
     &     dim1, dim2, nbas0
      call int_dip_ga(bas_hd, bas_hd, g_dp(1), g_dp(2), g_dp(3))

      g_hi = ga_create_atom_blocked(geom_hd, bas_hd, 'H_int')
      call ga_zero(g_hi)
      do i = 1, 3
cforcheck         if (ndid==0) write(luout, '("Dipole :", i1)') i
cforcheck         call ga_print(g_dp(i))
         call ga_add(1.0d0, g_hi, ef_vec(i), g_dp(i), g_hi)
      enddo
cforcheck      if (ndid==0) write(luout, '("h_int :")')
cforcheck      call ga_print(g_hi)

      do i = 1, 3
         if (.not. ga_destroy(g_dp(i)))
     &        call errquit('WFAT: Cannot destroy '//dipname(i)//'.',
     &        1, GA_ERR)
      enddo      
      

c     ==== RTDB reads ====
      ord = -100
      maxord = -100
      if (.not. rtdb_cget(rtdb, 'wfat:perturb_type', 1, type))
     &     call errquit('WFAT: Cannot read perturb_type', 1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:perturb_dgnthr', mt_dbl,
     &     1, dgnthr))
     &     call errquit('WFAT: Cannot read perturb_dgnthr', 1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:perturb_offthr', mt_dbl,
     &     1, offthr))
     &     call errquit('WFAT: Cannot read perturb_offthr', 1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:perturb_order', mt_int, 1, ord))
     &     call errquit('WFAT: Cannot read perturb_order', 1, RTDB_ERR)


c     ==== Allocations ====
      if (.not. ga_create(mt_dbl, nbas0, (ord+1)*nch, 'corct_MO', 0, 0,
     &     g_cmo))
     &     call errquit('WFAT: Cannot create corct_MO', 1, GA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nch, 'max_c', l_mx, k_mx))
     &     call errquit('WFAT: Cannot allocate max_c', 1, MA_ERR)

      
c     ==========================
c     ==== Brillouin-Wigner ====
c     ==========================      
      if (type == 'bw') then
c        == RTDB reads ==
         if (.not. rtdb_get(rtdb, 'wfat:perturb_e1thr', mt_dbl,
     &        1, e1thr))
     &        call errquit('WFAT: Cannot read perturb_e1thr',
     &        1, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'wfat:perturb_itrthr', mt_dbl,
     &        1, itrthr))
     &        call errquit('WFAT: Cannot read perturb_itrthr',
     &        1, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'wfat:perturb_maxitr', mt_int,
     &        1, maxitr))
     &        call errquit('WFAT: Cannot read perturb_maxitr',
     &        1, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'wfat:perturb_ordthr', mt_dbl,
     &        1, ordthr))
     &        call errquit('WFAT: Cannot read perturb_ordthr',
     &        1, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'wfat:perturb_maxord', mt_int,
     &        1, maxord))
     &        call errquit('WFAT: Cannot read perturb_maxord',
     &        1, RTDB_ERR)

c        == Allocations ==
         if (.not. ma_alloc_get(mt_int, nch, 'numitr', l_nitr, k_nitr))
     &        call errquit('WFAT: Cannot allocate numitr', 1, MA_ERR)
         if (.not. ma_alloc_get(mt_int, nch, 'numord', l_nord, k_nord))
     &        call errquit('WFAT: Cannot allocate numord', 1, MA_ERR)
         if (.not. ma_alloc_get(mt_dbl, (maxord+1)*nch, 'erg',
     &        l_e, k_e))
     &        call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)
         if (.not. ma_alloc_get(mt_int, 2*nch, 'max_c_i', l_mxi, k_mxi))
     &        call errquit('WFAT: Cannot allocate max_c_i', 1, MA_ERR)
         
         call wfat_perturb_bw(geom_hd, bas_hd, nbas0, nspin0, g_mov,
     &        dbl_mb(k_moerg), g_hi, nch, orbid, sid, dgnthr, offthr, 
     &        e1thr, itrthr, maxitr, ordthr, maxord, ord, 
     &        int_mb(k_nitr), int_mb(k_nord), dbl_mb(k_mx),
     &        int_mb(k_mxi), dbl_mb(k_e), g_cmo)

         if (.not. ma_alloc_get(mt_dbl, nch, 'corct_erg', l_ce, k_ce))
     &        call errquit('WFAT: Cannot allocate corct_erg', 1, MA_ERR)
         do i = 1, nch
            i_e1 = k_e + (i-1)*(maxord+1)
            i_e2 = i_e1 + (maxord+1) - 1
            dbl_mb(k_ce+i-1) = sum(dbl_mb(i_e1:i_e2))
         enddo
         call wfat_perturb_analconv_bw(nbas0, nspin0, dbl_mb(k_moerg),
     &        nch, sid, dbl_mb(k_ce), int_mb(k_nitr), int_mb(k_nord),
     &        dbl_mb(k_mx), int_mb(k_mxi))
         if (.not. ma_free_heap(l_ce))
     &        call errquit('WFAT: Cannot deallocate corct_erg',
     &        1, MA_ERR)

c        == Deallocations ==
         if (.not. ma_free_heap(l_nitr))
     &        call errquit('WFAT: Cannot deallocate numitr', 1, MA_ERR)
         if (.not. ma_free_heap(l_nord))
     &        call errquit('WFAT: Cannot deallocate numord', 1, MA_ERR)         

c     ===============================         
c     ==== Rayleigh-Schroedinger ====
c     ===============================
      else if (type == 'rs') then
c        == Allocations ==
         if (.not. ma_alloc_get(mt_dbl, 2*nch, 'erg', l_e, k_e))
     &        call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)
         if (.not. ma_alloc_get(mt_int, nch, 'max_c_i', l_mxi, k_mxi))
     &        call errquit('WFAT: Cannot allocate max_c_i', 1, MA_ERR)

         call wfat_perturb_rs(geom_hd, bas_hd, nbas0, nspin0, g_mov,
     &        dbl_mb(k_moerg), g_hi, nch, orbid, sid, dgnthr, offthr,
     &        dbl_mb(k_mx), int_mb(k_mxi), dbl_mb(k_e), g_cmo)

         call wfat_perturb_analconv_rs(nbas0, nspin0, dbl_mb(k_moerg),
     &        nch, orbid, sid, dbl_mb(k_mx), int_mb(k_mxi))
      endif

      
c     ==== Common analysis ====
      if (type == 'bw') then
         call wfat_perturb_analyse(rtdb, geom_hd, bas_hd, nbas0, nspin0, 
     &        nel, nch, sid, ord, maxord, g_hi, g_mov, g_cmo, 
     &        dbl_mb(k_e), tol2e, oskel, ixcfit, iolgc, rmode)
      else if (type == 'rs') then
         call wfat_perturb_analyse(rtdb, geom_hd, bas_hd, nbas0, nspin0, 
     &        nel, nch, sid, 1, 1, g_hi, g_mov, g_cmo, dbl_mb(k_e), 
     &        tol2e, oskel, ixcfit, iolgc, rmode)
      endif


c     ==== Deallocations ====
      if (.not. ga_destroy(g_cmo))
     &     call errquit('WFAT: Cannot destroy corct_MO', 1, GA_ERR)
      if (.not. ma_free_heap(l_e))
     &     call errquit('WFAT: Cannot deallocate erg', 1, MA_ERR)
      if (.not. ma_free_heap(l_mx))
     &     call errquit('WFAT: Cannot deallocate max_c', 1, MA_ERR)
      if (.not. ma_free_heap(l_mxi))
     &     call errquit('WFAT: Cannot deallocate max_c_i', 1, MA_ERR)
      
      end
c===============================================================      


c===============================================================
c===============================================================
      subroutine wfat_diag_perturb(nb, ns, g_v, erg, g_hi, thr, othr)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"      
      integer nb
      integer ns
      integer g_v(ns)               ! [in/output], vectors in which degeneracies are searched for. COLUMNS OF g_v(1:ns) MUST BE ORTHONORMAL. 
      double precision erg(nb,ns)   ! [input], unperturbed energies associated with g_v
      integer g_hi                  ! [input], interaction Hamiltonian in standard GTO basis
      double precision thr          ! [input], degeneracy threshold
      double precision othr         ! [input], off-diagonal threshold
      
      logical iscalc(nb, ns), isdiag
      integer g_v0, l_hi, k_hi, i_hi, l_ev, k_ev
      integer is, i, j, k1, k2, j1, j2
      integer ndg, dgid(nb)
      

      call ga_sync()

      iscalc = .false.
      if (.not. ga_duplicate(g_v(1), g_v0, 'scratch'))
     &     call errquit('WFAT: Cannot duplicate scratch.', 1, GA_ERR)
      call ga_zero(g_v0)
      


      do is = 1, ns
         do i = 1, nb
            if (iscalc(i, is)) cycle

c           ==== Search for degeneracies ====
            ndg = 0
            do j = 1, nb
               if (abs(erg(i,is) - erg(j,is)) < thr) then
                  ndg = ndg + 1
                  dgid(ndg) = j
               endif
            enddo

c           ==== If degenerate, construct H_int in degenerate subspace  ====
            if (ndg > 1) then
               if (.not. ma_alloc_get(mt_dbl, ndg*ndg, 'h_int',
     &              l_hi, k_hi))
     &              call errquit('WFAT: Cannot allocate h_int',
     &              1, MA_ERR)
               if (.not. ma_alloc_get(mt_dbl, ndg, 'h_eval',
     &              l_ev, k_ev))
     &              call errquit('WFAT: Cannot allocate h_eval',
     &              1, MA_ERR)
               
               isdiag = .true.
               do k2 = 1, ndg
                  j2 = dgid(k2)
                  do k1 = 1, ndg
                     j1 = dgid(k1)
                     call ga_matmul_patch(
     &                    'n', 'n', 1.0d0, 0.0d0,
     &                    g_hi,    1, nb, 1,  nb,
     &                    g_v(is), 1, nb, j2, j2,
     &                    g_v0,    1, nb, 1,  1)

                     i_hi = k_hi + (k2-1)*ndg + k1 - 1
                     dbl_mb(i_hi) = ga_ddot_patch(
     &                          g_v(is), 'n', 1, nb, j1, j1,
     &                          g_v0,    'n', 1, nb, 1,  1)
                     
                     if (k1 /= k2 .and. abs(dbl_mb(i_hi)) > othr) then
                        isdiag = .false.
                     endif
                  enddo
               enddo
               call ga_scale_patch(g_v0, 1, nb, 1, 1, 0.0d0)

c              ==== Diagonalize H_int if not diagonal ====
               if (.not. isdiag) then
                  call wfat_dsyev_wrap(ndg, dbl_mb(k_hi), dbl_mb(k_ev))

c                 ==== Construct new degenerate vectors in basis rep. ====
                  do k2 = 1, ndg        ! Loop over new vectors
                     j2 = dgid(k2) 
                     do k1 = 1, ndg     ! Sum over old vectors
                        j1 = dgid(k1)

                        i_hi = k_hi + (k2-1)*ndg + k1 - 1
                        call ga_add_patch(
     &                       1.0d0,        g_v0,    1, nb, j2, j2,
     &                       dbl_mb(i_hi), g_v(is), 1, nb, j1, j1,
     &                       g_v0, 1, nb, j2, j2)
                     enddo
                  enddo

c                 ==== Replace the old with the new degenerate vectors ====
                  do k1 = 1, ndg
                     j1 = dgid(k1)
                     call ga_copy_patch('n',
     &                    g_v0, 1, nb, j1, j1, g_v(is), 1, nb, j1, j1)
                  enddo
               endif

               if (.not. ma_free_heap(l_hi))
     &              call errquit('WFAT: Cannot deallocate h_int',
     &              1, MA_ERR)
               if (.not. ma_free_heap(l_ev))
     &              call errquit('WFAT: Cannot deallocate h_eval',
     &              1, MA_ERR)
               
               do k1 = 1, ndg
                  iscalc( dgid(k1), is ) = .true.
               enddo

c           ==== If non-degenerate, do nothing ====
            else
               iscalc(i, is) = .true. 
            endif
            
         enddo
      enddo
      

      if (.not. ga_destroy(g_v0))
     &     call errquit('WFAT: Cannot destroy scratch.', 1, GA_ERR)

      call ga_sync()
      
      end
c===============================================================
