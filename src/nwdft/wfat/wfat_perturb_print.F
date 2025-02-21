c===============================================================
c===============================================================
      subroutine wfat_perturb_info(rtdb, bt, gm, ef)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
      double precision bt
      double precision gm
      double precision ef

      character*21 type0
      character*2 type
      integer ord, maxord, maxitr, ndid
      double precision dgnthr, offthr, e1thr, itrthr, ordthr


      ndid = ga_nodeid()

      
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

      if (type == 'bw') then
         type0 = 'Brillouin-Wigner'
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
      else if (type == 'rs') then
         type0 = 'Rayleigh-Schroedinger'
      else
         call errquit('WFAT: The type of perturbation analysis '//
     &        'is invalid.', 1, UERR)         
      endif
         

      if (ndid == 0) then
         write(luout, '(" ")')
         write(luout, '(20x, a26, " = ", a)') "type", type0
         write(luout, '(20x, a26, " = ", i2)') "correction order", ord
         write(luout, '(20x, a26, " = ", f10.6, " deg.")') "beta", bt
         write(luout, '(20x, a26, " = ", f10.6, " deg.")') "gamma", gm
         write(luout, '(20x, a26, " = ", es10.3, " a.u.")')
     &        "field strength", ef
         write(luout, '(20x, a26, " = ", es14.6, " Hartree")')
     &        "degeneracy threshold", dgnthr
         write(luout, '(20x, a26, " = ", es14.6, " Hartree")')
     &        "off-diagonal threhsold", offthr
         if (type == 'bw') then
            write(luout, '(20x, a26, " = ", es14.6, " Hartree")')
     &        "energy iteration threshold", itrthr
            write(luout, '(20x, a26, " = ", i3)')
     &           "max. energy iterations", maxitr
            write(luout, '(20x, a26, " = ", es14.6, " Hartree")')
     &           "correction order threshold", ordthr
            write(luout, '(20x, a26, " = ", i3)')
     &           "max. correction order", maxord
            write(luout, '(20x, a26, " = ", es14.6, " Hartree")')
     &           "guess energy threshold", e1thr
         endif
      endif
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_perturb_noinfo()

      implicit none
#include "stdio.fh"      


      write(luout, '(" ")')
      write(luout, 50)
 50   format(5x,
     &     "Nothing to be done for perturbation analysis, either ",
     &     "because the field strength is ", /,
     &     5x, 
     &     "too small or only structure factor calculation is ",
     &     "requested, in which case the ", /,
     &     5x,
     &     "input for field strength is omitted.")
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_perturb_potprint(rtdb, geom, bas, nb, ns, ne, 
     &           g_mo, tol2e, oskel, ixcfit, iolgc, g_hi)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "wfat.fh"
      integer rtdb            ! [input]
      integer geom            ! [input]
      integer bas             ! [input]
      integer nb              ! [input]
      integer ns              ! [input]
      integer ne(ns)          ! [input]
      integer g_mo(ns)        ! [input]
      double precision tol2e  ! [input]
      logical oskel           ! [input]
      integer ixcfit          ! [input]
      logical iolgc           ! [input] 
      integer g_hi            ! [input]

      character*255 fname
      integer g_veff(ns), g_veff2(ns)
      integer l_veff, k_veff, l_mo, k_mo
      integer ny, nz, type0, m0, n0
      integer is, i_veff, i_mo, ndid, iy, iz
      double precision ylims(2), zlims(2), vpt(3), bohr2ang
      double precision v_val, y_low, y_upp, dy, z_low, z_upp, dz

      parameter (bohr2ang=0.52917721065d0,
     &     y_low=-15.0d0, y_upp=15.0d0, z_low=-15.0d0, z_upp=15.0d0)            ! Default values *_upp and *_low are in Bohr.


      call ga_sync()
      ndid = ga_nodeid()

c     ==== Construct perturbed potential matrix ====
      do is = 1, ns
         if (.not. ga_duplicate(g_hi, g_veff(is), 'v_eff'))
     &        call errquit('WFAT: Cannot create v_eff', 1, GA_ERR)
         if (.not. ga_duplicate(g_hi, g_veff2(is), 'v_eff_MO'))
     &        call errquit('WFAT: Cannot create v_eff_MO', 1, GA_ERR)
         call ga_zero(g_veff(is))
         call ga_zero(g_veff2(is))
      enddo
      call wfat_get_dft_veffmat(rtdb, geom, bas, nb, ns, ne, g_mo,
     &     tol2e, oskel, ixcfit, iolgc, g_veff)      
      do is = 1, ns
         call ga_dadd(1.0d0, g_veff(is), 1.0d0, g_hi, g_veff(is))
      enddo
      if (ndid==0) write(luout, '("Interaction Hamiltonian:")')
      call ga_print(g_hi)


c     ==== Represent potential matrix in MO basis ====
      if (.not.ma_alloc_get(mt_dbl, nb*nb*ns, 'veff', l_veff, k_veff))
     &     call errquit('WFAT: Cannot allocate veff', 1, MA_ERR)
      do is = 1, ns
         call wfat_ga_vrepresent(g_mo(is), g_veff(is), g_veff2(is))
         i_veff = k_veff + (is-1)*nb*nb
         call ga_get(g_veff2(is), 1, nb, 1, nb, dbl_mb(i_veff), nb)
      enddo

      
c     ==== Y range for printing ====
      if (.not. rtdb_get(rtdb, 'wfat:res_mode:ny_veffprint', mt_int,
     &     1, ny))
     &     ny = 150
      if (rtdb_get(rtdb, 'wfat:res_mode:ylims_veffprint', mt_dbl,
     &     2, ylims)) then
         ylims = ylims / bohr2ang
      else
         ylims = (/ y_low, y_upp /)
      endif

      
c     ==== Z range for printing ====
      if (.not. rtdb_get(rtdb, 'wfat:res_mode:nz_veffprint', mt_int,
     &     1, nz))
     &     nz = 150
      if (rtdb_get(rtdb, 'wfat:res_mode:zlims_veffprint', mt_dbl,
     &     2, zlims)) then
         zlims = zlims / bohr2ang
      else
         zlims = (/ z_low, z_upp /)
      endif
      
            
c     ==== Print in yz-plane ====
      if (.not.ma_alloc_get(mt_dbl, nb*nb*ns, 'MO', l_mo, k_mo))
     &     call errquit('WFAT: Cannot allocate MO', 1, MA_ERR)
      vpt(1) = 0.0d0
      dy = (ylims(2) - ylims(1)) / (ny-1)
      dz = (zlims(2) - zlims(1)) / (nz-1)

      call util_file_name('weff', .false., .false., fname)
      if (ndid==0) open(WF_PERTV_UNIT, file=fname, action='write')
      do iz = 1, nz
         vpt(3) = zlims(1) + (iz-1)*dz
         do iy = 1, ny
            vpt(2) = ylims(1) + (iy-1)*dy

            if (ndid==0) write(WF_PERTV_UNIT, '(2es14.5)', advance='no')
     &           vpt(2)*bohr2ang, vpt(3)*bohr2ang
            do is = 1, ns
               i_mo = k_mo + (is-1)*nb*nb
               call wfat_eval_mo0_1pt(geom, bas, vpt, nb, g_mo(is),
     &              dbl_mb(i_mo))
               i_veff = k_veff + (is-1)*nb*nb
               call wfat_1e_expval_1eopr(nb, dbl_mb(i_mo),
     &              dbl_mb(i_veff), v_val)

               if (ndid==0)
     &              write(WF_PERTV_UNIT, '(" ", es14.6)', advance='no')
     &              v_val
            enddo
            if (ndid==0) write(WF_PERTV_UNIT, '(" ")')
            
         enddo
         if (ndid==0) write(WF_PERTV_UNIT, '(" ")')
      enddo
      if (ndid==0) close(WF_PERTV_UNIT)


c     ==== Deallocations ====
      do is = 1, ns
         if (.not. ga_destroy(g_veff(is)))
     &        call errquit('WFAT: Cannot deallocate v_eff',
     &        1, GA_ERR)
         if (.not. ga_destroy(g_veff2(is)))
     &        call errquit('WFAT: Cannot deallocate v_eff_MO',
     &        1, GA_ERR)
      enddo
      if (.not. ma_free_heap(l_veff))
     &     call errquit('WFAT: Cannot deallocate veff', 1, MA_ERR)
      if (.not. ma_free_heap(l_mo))
     &     call errquit('WFAT: Cannot deallocate MO', 1, MA_ERR)


      call ga_sync()

      end
c===============================================================      
