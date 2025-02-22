c===============================================================
c===============================================================
      subroutine wfat_calc_ti(rtdb, geom, pr_srfc, pr_rate, fieldm,
     &           nch, orbid0, spinid, erg, chmmt, nqn, nid, mid, lmax, 
     &           zi, intgl, wfat1, nsub, maxn, dgs, ndgs, dps, pbs,
     &           wfat1_exc, intgl_j, intgl_k1, intgl_k2)

c     DESCRIPTION:
c     nsub = It is only used when wfat1=.true.. The number of degenerate subspaces,
c            a nondegenerate state is considered a degenerate state in a subspace of
c            dimension 1.
c     maxn = It is only used when wfat1=.true.. The largest dimension of all degenerate
c            subspaces, i.e. maxn=maxval(ndgs).
c     dgs = It is only used when wfat1=.true.. The channel ID arranged into degenerate
c           subspaces (columns). Each subspace (column) must only have states of the same
c           spin. However, note that dgs has no information about the spins of its
c           columns. This information is not needed because dgs contains the indices of the
c           channels which already encode the spin of that channel. The range of its 
c           elements is 1 <= dgs(i,j) <= nch. Its minimum size ismaxn*nsub.
c     ndgs = It is only used when wfat1=.true.. The dimension of each degenerate subspace.
c            Its minimum size is nsub.
c     dps = It is only used when wfat1=.true.. The electronic dipole moment matrix in
c           molecular frame in each degenerate subspace. The 3 indices in the 3rd
c           dimension correspond to x, y, z. Its minimum size is maxn*maxn*3*nsub.
c     pbs = It is only used when wfat1=.true.. The static polarizability matrix in 
c           molecular frame in each degenerate subspace. The 6 indices in the 3rd 
c           dimension correspond to xx, yy, zz, xy, xz, yz. Its minimum size is
c           maxn*maxn*6*nsub.

      implicit none
#include "global.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb                    ! [input]
      integer geom                    ! [input]
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_rate                 ! [input], if .true. then print ion yield.
      double precision fieldm
      integer nch                     ! [input], number of ionized channels.
      integer orbid0(nch)             ! [input], orbital ID of ionized channels.
      integer spinid(nch)             ! [input], spin value (1 or 2) of ionized channels.
      double precision erg(nch)       ! [input]
      double precision chmmt(3, nch)  ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax                    ! [input]
      double precision zi             ! [input]
      double complex intgl(1:*)       ! [input]
      logical wfat1                   ! [input]
      integer nsub                    ! [input]. Only used when wfat1=.true..
      integer maxn                    ! [input]. Only used when wfat1=.true..
      integer dgs(1:*)                ! [input]. Only used when wfat1=.true..
      integer ndgs(1:*)               ! [input]. Only used when wfat1=.true..
      double precision dps(1:*)       ! [input]. Only used when wfat1=.true..
      double precision pbs(1:*)       ! [input]. Only used when wfat1=.true..
      logical wfat1_exc(1:*)          ! [input]. Only used when wfat1=.true..
      double complex intgl_j(1:*)     ! [input]. Only used when wfat1=.true..
      double complex intgl_k1(1:*)    ! [input]. Only used when wfat1=.true..
      double complex intgl_k2(1:*)    ! [input]. Only used when wfat1=.true..
      
      logical gap, interm
      integer n_gn, nfc, nrt, nbt, ngm, nat, nsc, nsx
      integer npbd, n_jn, n_kn, n_on, ib, ig, ndid, i, iqn, i_on
      integer l_gn, l_fc, l_rt, l_sc, l_sx, l_pbd
      integer l_jn, l_kn, l_on
      integer k_gn, k_fc, k_rt, k_sc, k_sx, k_pbd
      integer k_jn, k_kn, k_on
      double precision bt_start, bt_last, dbt, dbtr, bt
      double precision gm_start, gm_last, dgm, dgmr, gm
      double precision rate, dthr, ooo(2), ooo_(2)
      double precision deg2rad, pi_const

      parameter (deg2rad = 0.01745329251994d0)
      parameter (pi_const = 3.14159265359d0)


      ndid = ga_nodeid()
      if (ndid /= 0) return


c     ==== Obtaining info on beta ====
      call wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)
      dbtr = dbt * deg2rad


c     ==== Obtaining info on gamma ====
      call wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)
      dgmr = dgm * deg2rad

      
c     ==== Allocation of structure factors and ion yields ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 31, GEOM_ERR)
      call wfat_prep_rate(nat, nqn, nch, n_gn, l_gn, k_gn, nfc,
     &     l_fc, k_fc, nrt, l_rt, k_rt)

      
c     ==== Obtain field factor ====
      call wfat_get_field_factor(rtdb, nch, nqn, nid, mid, zi, fieldm,
     &     erg, dbl_mb(k_fc))
                  

c     ==== When WFAT1 is active ====
      if (wfat1) then
         call wfat_prepare_an(maxn, nat, nch, nqn, nsc, l_sc, k_sc, 
     &        nsx, l_sx, k_sx, npbd, l_pbd, k_pbd, n_jn, l_jn, k_jn,
     &        n_kn, l_kn, k_kn, n_on, l_on, k_on)
         if (.not. rtdb_get(rtdb, 'wfat:dthr', mt_dbl, 1, dthr))
     &        call errquit('WFAT: Failed reading dthr from rtdb',
     &        1, RTDB_ERR)
         do i = 1, nch
            do iqn = 1, nqn
               i_on = k_on + 2*nqn*(i-1) + 2*(iqn-1)
               call wfat_calc_o(nid(iqn), mid(iqn), erg(i), zi,
     &              dbl_mb(i_on))
            enddo
         enddo
         call wfat_print_o(nch, nqn, nid, mid, dbl_mb(k_on))
         if (.not. rtdb_get(rtdb, 'wfat:wfat1_interm', mt_log, 1,
     &        interm)) interm = .false.
      endif

      
c     ==== Initializing the printings ====
      if (ndid == 0)
     &     call wfat_print_field_factor(rtdb, fieldm, nch, nqn, erg, 
     &     nid, mid, zi, dbl_mb(k_fc))
      if (ndid == 0 .and. pr_srfc) then
         call wfat_initprint_structuref(nch, orbid0, spinid, nqn,
     &        nid, mid)
         if (wfat1) call wfat_initprint_srfc1(nch, orbid0, spinid, nqn,
     &        nid, mid)
      endif
      if (ndid == 0 .and. pr_rate)
     &     call wfat_initprint_rate(nch, orbid0, spinid, nqn, nid, mid)
      
      
c     ==== Calculate and print structure factor and ion yield ====
      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)
         
         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)

            if (ib /= nbt) then
               gap = .false.
            else
               gap = .true.
            endif
            
            if (pr_srfc .or. pr_rate) then
c              ==== Calculate structure factors ====
               if (wfat1) then
                  call wfat_get_srfc1(bt, gm, nat, nch, nqn, lmax, mid, 
     &                 erg, chmmt, intgl, intgl_j, intgl_k1, intgl_k2, 
     &                 nsub, maxn, dgs, ndgs, dthr, dps, pbs, .false.,
     &                 dbl_mb(k_sc), dcpl_mb(k_sx), dbl_mb(k_pbd), 
     &                 dcpl_mb(k_gn), dcpl_mb(k_jn), dcpl_mb(k_kn))
               else
                  call wfat_get_srfc0(geom, bt, gm, nch, erg, nqn, mid, 
     &                 lmax, chmmt, intgl, .false., dcpl_mb(k_gn))
               endif
               

c              ==== Print structure factor ====
               if (ndid == 0 .and. pr_srfc) then
                  call wfat_print_structuref(nat, bt, gm, nch, nqn,
     &                 dcpl_mb(k_gn), gap)
                  if (wfat1) call wfat_print_srfc1(nat, bt, gm, nch, 
     &                 nqn, dcpl_mb(k_jn), dcpl_mb(k_kn), gap)
               endif
            endif

            if (pr_rate) then
c              ==== Calculate ionization rate ====
               if (wfat1) then
                  call wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid,
     &                 mid, fieldm, erg, zi, dcpl_mb(k_gn),
     &                 dcpl_mb(k_jn), dcpl_mb(k_kn), dbl_mb(k_on), 
     &                 chmmt, dbl_mb(k_pbd), dbl_mb(k_fc), wfat1_exc, 
     &                 interm, dbl_mb(k_rt), rate)
               else
                  call wfat_calc_rate(nat, nch, nqn, dcpl_mb(k_gn), 
     &                 dbl_mb(k_fc), dbl_mb(k_rt), rate)
               endif
               
c              ==== Print ionization rate ====
               if (ndid == 0 .and. pr_rate)
     &              call wfat_print_rate(bt, gm, nch, nqn, 
     &              dbl_mb(k_rt), rate, gap)
            endif
            
         enddo
      enddo
      

c     ==== Terminating printings ====      
      if (ndid==0 .and. pr_srfc) then
         call wfat_endprint_structuref()
         if (wfat1) call wfat_endprint_srfc1()
      endif
      if (ndid==0 .and. pr_rate) call wfat_endprint_rate()

      
c     ==== Deallocations ====
      call wfat_erase_rate(l_gn, l_fc, l_rt)
      if (wfat1) call wfat_erase_an(l_sc, l_sx, l_pbd, l_jn, l_kn, l_on)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb
      double precision bt_start
      double precision bt_last
      integer nbt
      double precision dbt


      if (.not. rtdb_get(rtdb, 'wfat:beta_start', mt_dbl, 1, bt_start))
     &     call errquit('WFAT: Failed reading beta_start from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:beta_last', mt_dbl, 1, bt_last))
     &     call errquit('WFAT: Failed reading beta_last from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:nbeta', mt_int, 1, nbt))
     &     call errquit('WFAT: Failed reading nbeta from rtdb',
     &     1, INPUT_ERR)

      if (nbt > 1) then
         dbt = (bt_last - bt_start) / (nbt - 1)
      else
         dbt = 0.0d0
      endif

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb
      double precision gm_start
      double precision gm_last
      integer ngm
      double precision dgm
      

      if (.not. rtdb_get(rtdb, 'wfat:gamma_start', mt_dbl, 1, gm_start))
     &     call errquit('WFAT: Failed reading gamma_start from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:gamma_last', mt_dbl, 1, gm_last))
     &     call errquit('WFAT: Failed reading gamma_last from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:ngamma', mt_int, 1, ngm))
     &     call errquit('WFAT: Failed reading ngamma from rtdb',
     &     1, INPUT_ERR)

      if (ngm > 1) then
         dgm = (gm_last - gm_start) / (ngm - 1)
      else
         dgm = 0.0d0
      endif
      
      end
c===============================================================
      
      
c===============================================================
c===============================================================
      subroutine wfat_prepare_an(maxn, nat, nch, nqn, nsc, l_sc, k_sc,
     &           nsx, l_sx, k_sx, npbd, l_pbd, k_pbd, n_jn, l_jn, k_jn, 
     &           n_kn, l_kn, k_kn, n_on, l_on, k_on)
      
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer maxn           ! [input], 
      integer nat            ! [input], 
      integer nch            ! [input], 
      integer nqn            ! [input], 
      integer nsc            ! [output],
      integer l_sc, k_sc     ! [output],
      integer nsx            ! [output],
      integer l_sx, k_sx     ! [output],
      integer npbd           ! [output],
      integer l_pbd, k_pbd   ! [output],
      integer n_jn           ! [output],
      integer l_jn, k_jn     ! [output],
      integer n_kn           ! [output],
      integer l_kn, k_kn     ! [output],
      integer n_on           ! [output],
      integer l_on, k_on     ! [output],

      
      nsc = maxn**2 + maxn
      if (.not. MA_alloc_get(MT_DBL, nsc, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 2, MA_ERR)
      nsx = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nsx, 'scr-cpx', l_sx, k_sx))
     &     call errquit('WFAT: scr-cpx allocation failed', 2, MA_ERR)
      npbd = 3 * 3 * nch
      if (.not. MA_alloc_get(MT_DBL, npbd, 'pbd', l_pbd, k_pbd))
     &     call errquit('WFAT: pbd allocation failed', 2, MA_ERR)
      n_jn = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, n_jn, 'jn', l_jn, k_jn))
     &     call errquit('WFAT: an allocation failed', 2, MA_ERR)
      n_kn = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, n_kn, 'kn', l_kn, k_kn))
     &     call errquit('WFAT: an allocation failed', 2, MA_ERR)
      n_on = 2 * nqn * nch
      if (.not. MA_alloc_get(MT_DBL, n_on, 'on', l_on, k_on))
     &     call errquit('WFAT: on allocation failed', 2, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_erase_an(l_sc, l_sx, l_pbd, l_jn, l_kn, l_on)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_sc
      integer l_sx
      integer l_pbd
      integer l_jn
      integer l_kn
      integer l_on

      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: scratch deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_sx))
     &     call errquit('WFAT: scr-cpx deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_pbd))
     &     call errquit('WFAT: pbd deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_jn))
     &     call errquit('WFAT: jn deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_kn))
     &     call errquit('WFAT: kn deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_on))
     &     call errquit('WFAT: on deallocation failed', 3, MA_ERR)

      end
c===============================================================
      
      
