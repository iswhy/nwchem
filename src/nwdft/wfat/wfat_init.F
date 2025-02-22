c=======================================================
c=======================================================
      subroutine wfat_init(rtdb)

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"      
#include "geom.fh"
#include "bas.fh"
#include "sym.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "wfat.fh"      
#include "wfat_qm_param.fh"
#include "wfat_rs.fh"
      integer rtdb

      logical use_nwxc, rtdb_stat
      integer is, i, ioc, ieg, basnormid, noio, ndid, offs
      character*20 basis_name, trans_name
      character*2 spin_name
      character*7 mov_name0, mov_name
      parameter(mov_name0='movec')
      double precision zn, rscent0(3), dgnthr
      
      logical movecs_read, nwxc_has_cam
      external movecs_read, nwxc_has_cam
      logical int_normalize
      external int_normalize
      integer ga_create_atom_blocked
      external ga_create_atom_blocked


      ndid = ga_nodeid()

cfordebug      if (.not. rtdb_print(rtdb, .false.))
cfordebug     &     call errquit('WFAT: Cannot print RTDB.', 1, RTDB_ERR)
      
      
c     ==== Check existence, format, and ====
c     ====     necessities of inputs    ====
      if (ga_nodeid() == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         call wfat_process_input(rtdb)
         rtdb_stat = rtdb_parallel(.true.)
      endif
      call ga_sync()


      if (ndid == 0) then
c        ==== WFAT module title ====
         call wfat_module_title(luout)
         
c        ==== WFAT level info ====
         rtdb_stat = rtdb_parallel(.false.)
         call wfat_print_level(rtdb)
         rtdb_stat = rtdb_parallel(.true.)
      endif

      
c     ==== Printing units ====
      WF_SF_PRINT_READY = .false.
      WF_IY_PRINT_READY = .false.
      WF_SF1_PRINT_READY(1:2) = .false.
      WF_SF1_UNIT = (/ WF_JN_UNIT, WF_KN_UNIT /)

      
c     ==== Geometry ====
      if (.not. geom_create(geom_hd, 'wfat_geom'))
     &     call errquit('WFAT: Error in creating geometry object.',
     &     1, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom_hd, 'geometry'))
     &     call errquit('WFAT: Error in loading geometry from rtdb',
     &     2, GEOM_ERR)
      if (.not. geom_ncent(geom_hd, natm))
     &     call errquit('WFAT: Getting natm failed', 51, GEOM_ERR)

c     ==== Geometry info ====
      if (ndid == 0) then
         write(luout, '(/)')
         if (.not. geom_print(geom_hd))
     &        call errquit('WFAT: Cannot print geometry info', 1,
     &        GEOM_ERR)
      endif


c      ==== No. of particles and spins ====
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, nspin0))
     &     call errquit('WFAT: Failed reading ipol from RTDB.', 1,
     &     RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, nel))
     &     call errquit('WFAT: Failed reading noc from RTDB.', 1,
     &     RTDB_ERR)      

      if (nspin0 == 1) then
         nelec0 = 2 * nel(1)
      else if (nspin0 == 2) then
         nelec0 = sum(nel)
      endif

      if (.not. geom_nuc_charge(geom_hd, zn))
     &     call errquit('WFAT: Cannot get nuclear charge', 1, GEOM_ERR)
      zion = zn - nelec0 + 1


c     ==== Basis ====
      if (.not. bas_create(bas_hd, 'wfat_bas'))
     &     call errquit('WFAT: Error in creating basis object.',
     &     1, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom_hd, bas_hd, 'ao basis'))
     &     call errquit('WFAT: Error in loading basis from rtdb',
     &     3, BASIS_ERR)
      if(.not. bas_name(bas_hd, basis_name, trans_name))
     &     call errquit('WFAT: Error in getting basis name.',
     &     1, BASIS_ERR)
      if (.not. bas_numbf(bas_hd, nbas0))
     &     call errquit('WFAT: Getting nbas0 failed in ', 53, UERR)

      if (.not. int_normalize(rtdb, bas_hd))
     &     call errquit('WFAT: Error in normalizing contraction '//
     &     'coefficients.', 1, BASIS_ERR)
      if (.not. bas_norm_get(bas_hd, basnormid))
     &     call errquit('WFAT: Error getting basis norm ID',
     &     1, BASIS_ERR)

c     ==== Basis info ====
      if (ndid == 0) then
         write(luout, '(/)')
         if (.not. bas_print(bas_hd))
     &        call errquit('WFAT: Cannot print basis info', 1,
     &        BASIS_ERR)
      endif
      

c     ==== Molecular orbitals ====
      nocc0 = nbas0 * nspin0
      if (.not. ma_alloc_get(mt_dbl, nocc0, 'moocc', l_moocc, k_moocc))
     &     call errquit('WFAT: Allocating occnum failed', 1, MA_ERR)
      nerg0 = nbas0 * nspin0
      if (.not. ma_alloc_get(mt_dbl, nerg0, 'moerg', l_moerg, k_moerg))
     &     call errquit('WFAT: Allocating erg failed', 2, MA_ERR)

      do is = 1, nspin0
         if (is == 1) then
            spin_name = '_a'
            mov_name = trim(adjustl(mov_name0))//spin_name
            if (.not. ga_create(mt_dbl, nbas0, nbas0, mov_name, 0, 0,
     &           g_mov(is)))
     &           call errquit ("WFAT: couldnt create "//mov_name,
     &           0, GA_ERR)
         else if (is == 2) then
            spin_name = '_b'
            mov_name = trim(adjustl(mov_name0))//spin_name
            if (.not. ga_duplicate(g_mov(1), g_mov(2), mov_name))
     &           call errquit('WFAT: Couldnt create '//mov_name,
     &           0, GA_ERR)
         endif
         
         call ga_zero(g_mov(is))
      enddo

      if (.not. rtdb_cget(rtdb, 'wfat:movec_name', 1, mov_in))
     &     call errquit('WFAT: Cannot read movec_name from RTDB',
     &     1, RTDB_ERR)
      do is = 1, nspin0
         ioc = k_moocc + (is-1) * nbas0
         ieg = k_moerg + (is-1) * nbas0
         if (.not. movecs_read(mov_in, is, dbl_mb(ioc),    !! 1)
     $             dbl_mb(ieg), g_mov(is)))
     $        call errquit("WFAT: Failed to read movecs data", 0, UERR)
      enddo

      
c     ==== Calculate charge state and multiplicity ====
      charge = zn - nelec0
      if (nspin0 == 1) then
         multip = 1
      else if (nspin0 == 2) then
         multip = nel(1) - nel(2) + 1
      endif
      

c     ==== Quantum system info ====
      call wfat_print_qmsys(rtdb)
      call ga_sync()
                  
      
c     ==== Detect degeneracies ====
      if (ga_nodeid() == 0) then
         write(luout, '(2/)')
         offs = 0
         if (nspin0 == 1) offs = 12
         call util_print_centered(luout, 'Degeneracy analysis', 40-offs,
     &        .true.)
         write(luout, '(" ")')
      endif
      if (.not. ma_alloc_get(mt_int, nbas0*nspin0, 'modgid',
     &     l_modgid, k_modgid))
     &     call errquit('WFAT: Allocating dgid failed', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_int, nbas0*nspin0, 'modgn',
     &     l_modgn, k_modgn))
     &     call errquit('WFAT: Allocating dgn failed', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_int, nspin0, 'norb_dg',
     &     l_norbdg, k_norbdg))
     &     call errquit('WFAT: Allocating norb_dg failed', 1, MA_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:dgnthr', mt_dbl, 1, dgnthr))
     &     call errquit('WFAT: Failed reading dgnthr from RTDB.', 1,
     &     RTDB_ERR)
      call wfat_detect_degen(nbas0, nspin0, dbl_mb(k_moerg), dgnthr,
     &     int_mb(k_modgid), int_mb(k_modgn), int_mb(k_norbdg))
      if (ga_nodeid()==0) call wfat_print_degen(nbas0, nspin0,
     &     dbl_mb(k_moerg), dbl_mb(k_moocc), dgnthr, int_mb(k_modgid),
     &     int_mb(k_modgn))
      call ga_sync()
      

c     ==== Molecular integrals ====
      call int_init(rtdb, 1, bas_hd)
      if (.not. rtdb_get(rtdb, 'dft:tol2e', mt_dbl, 1, tol2e))
     &     call errquit('WFAT: Failed reading tol2e from RTDB.', 1,
     &     RTDB_ERR)
      call scf_get_fock_param(rtdb, tol2e)
      call fock_force_direct(rtdb)          ! Not sure why but making this call speed up walltime.
      call schwarz_init(geom_hd, bas_hd)

      
c     ==== XC functionals ====
      call nwxc_rtdb_load(rtdb, 'dft', use_nwxc)
      if (.not. use_nwxc)
     &     call errquit('WFAT: NWXC library for the density '//
     &     'functionals must be used when WFAT task is requested.',
     &     1, UERR)


c     ==== Use of CAM functional ====
      if (nwxc_has_cam()) then
         call wfat_info_banner(luout, .true.)
         if (ndid==0) write(luout, '(5x, "Coulomb-attenuated method '//
     &        'will be used in the computation of the exact '//
     &        'exchange.")')
         if (.not. geom_create(rsgeom, rsgeom_name))
     &        call errquit('WFAT: Cannot create rsgeom', 1, GEOM_ERR)
         rscent0 = (/ 0.0d0, 5.0d0, 0.0d0 /)
         if (.not. geom_cart_set(rsgeom, 1, rstag, rscent0, rschg))
     &        call errquit('WFAT: Cannot set rsgeom', 1, GEOM_ERR)
         if (.not. geom_rtdb_store(rtdb, rsgeom, rsgeom_name))
     &        call errquit('WFAT: Cannot store rsgeom', 1, GEOM_ERR)
         if (.not. bas_create(rsbas, rsbas_name))
     &        call errquit('WFAT: Cannot create rsbas', 1, BASIS_ERR)
         if (.not. bas_rtdb_load(rtdb, rsgeom, rsbas, rsbas_name))
     &        call errquit('WFAT: Cannot load rsbas', 1, BASIS_ERR)
      endif
      

c     ==== Miscellaneous ====
cOLD      call wfat_set_fockreplicate(1)
      if (.not. rtdb_get(rtdb, 'dft:skeleton', mt_log, 1, oskel))
     &     oskel = sym_number_ops(geom_hd) .gt. 0
      if (.not. rtdb_get(rtdb, 'dft:noio', mt_int, 1, noio))
     &     call errquit('WFAT: Failed reading noio from RTDB.', 1,
     &     RTDB_ERR)
      iolgc = .true.
      if (noio.eq.1) iolgc = .false.
      ixcfit = 0
      call wfat_reset_lfqlib()        ! Reset Gauss-Laguerre nodes and weights for the calculations of the various U(a,b,z).
      

c     ==== Store common integrals ====
      gwf_ovl = ga_create_atom_blocked(geom_hd, bas_hd, 'overlap')
      call ga_zero(gwf_ovl)
      call int_1e_ga(bas_hd, bas_hd, gwf_ovl, 'overlap', oskel)
      gwf_kin = ga_create_atom_blocked(geom_hd, bas_hd, 'kinetic')
      call ga_zero(gwf_kin)
      call int_1e_ga(bas_hd, bas_hd, gwf_kin, 'kinetic', oskel)
      gwf_pot = ga_create_atom_blocked(geom_hd, bas_hd, 'potential')
      call ga_zero(gwf_pot)
      call int_1e_ga(bas_hd, bas_hd, gwf_pot, 'potential', oskel)
      

c     ==== Printing important information ====
      call wfat_print_info(rtdb)
      
      end
c=======================================================


      

c=======================================================
c=======================================================
      subroutine wfat_set_fockreplicate(b)

      implicit none
#include "cfock.fh"
#include "wfat_qm_param.fh"
#include "errquit.fh"
#include "stdio.fh"      
      integer b


c      write(lakkuout, '("b = ", a5)') b

      if (b == 1) then
         oreplicated_prev = oreplicated
         oreplicated = .false.
      else if (b == 0) then
         oreplicated = oreplicated_prev
      else
         call errquit('WFAT: Argument ''b'' of wfat_set_fockreplicate'//
     &        ' must be either ''1'' or ''0''.', 1, UERR)
      endif
      

      end
c=======================================================      


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_init_getcdft(rtdb, bas_hd, nspin, geom_hd, natm)
cOLD
cOLD      implicit none
cOLDc#include "cdft.fh"
cOLDc#include "geom.fh"
cOLD#include "errquit.fh"
cOLD#include "stdio.fh"
cOLD#include "rtdb.fh"
cOLD#include "mafdecls.fh"
cOLD      integer rtdb
cOLD      integer bas_hd
cOLD      integer nspin
cOLD      integer geom_hd
cOLD      integer natm
cOLD
cOLD      integer ipol
cOLD      
cOLD      
cOLDc      bas_hd = ao_bas_han
cOLDc      geom_hd = geom
cOLD      
cOLD
cOLD      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))
cOLD     &     call errquit('WFAT: Failed reading ipol from RTDB.', 1,
cOLD     &     RTDB_ERR)
cOLD      nspin = ipol
cOLD      
cOLDc OLD      if (.not. geom_check_handle(geom, 'not valid in init'))
cOLDc OLD     &     call errquit('WFAT: Getting geom failed', 1, GEOM_ERR)
cOLDc OLD      
cOLDc OLD      write(luout, '("geom handle in init = ", i4)') geom
cOLDc OLD      if (.not. geom_ncent(geom, natm))
cOLDc OLD     &     call errquit('WFAT: Getting natm failed', 51, GEOM_ERR)
cOLD      
cOLD      end
cOLDc=======================================================            


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_init_getcscf(rtdb, ch, nel, g_mov, mov_in)
cOLD
cOLD      implicit none
cOLD#include "cscf.fh"
cOLD#include "stdio.fh"
cOLD#include "rtdb.fh"
cOLD#include "inp.fh"      
cOLD#include "errquit.fh"
cOLD      integer rtdb
cOLD      double precision ch
cOLD      integer nel(2)
cOLD      integer g_mov(2)
cOLD      character*256 mov_in, rtdb_name
cOLD      integer i, pos1, pos2, rtlen, mvlen
cOLD      
cOLD
cOLD      ch = charge
cOLD      nel(1) = nalpha
cOLD      nel(2) = nbeta
cOLD      g_mov = g_movecs
cOLDc      mov_in = movecs_in
cOLDc      mov_in = movecs_out
cOLD
cOLD      if(.not. rtdb_getfname(rtdb, rtdb_name)) 
cOLD     &     call errquit()
cOLD      if (.not. rtdb_cget(rtdb, 'dft:output vectors', 1, mov_in))
cOLD     &     call errquit('WFAT: Failed reading tol2e from RTDB.', 1,
cOLD     &     RTDB_ERR)
cOLD
cOLD      rtlen = inp_strlen(rtdb_name)
cOLD      mvlen = inp_strlen(mov_in)
cOLD      pos1 = 0
cOLD      pos2 = 0
cOLD      do i = 1, max(rtlen, mvlen)
cOLD         if (i <= rtlen) then
cOLD            if (rtdb_name(i:i) == '/' ) pos1 = i
cOLD         endif
cOLD         if (i <= mvlen) then   
cOLD            if (mov_in(i:i) == '/') pos2 = i
cOLD         endif
cOLD      enddo
cOLD      mov_in = rtdb_name(1:pos1)//mov_in(pos2+1:mvlen)
cOLD      
cOLD      write(luout, '("Vector file : ", a256)') adjustl(mov_in)
cOLD      
cOLDcOLD      write(luout, '("movecs: ", i12)') g_movecs(1)
cOLD      
cOLD      end
cOLDc=======================================================


c=======================================================
c=======================================================
      subroutine wfat_terminate(rtdb)

      implicit none
#include "rtdb.fh"      
#include "errquit.fh"
#include "mafdecls.fh"      
#include "wfat_qm_param.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "wfat_rs.fh"
      integer rtdb
      
      character*2 spin_name
      character*7 mov_name0, mov_name
      parameter(mov_name0='movec')
      logical ll
      integer is

      logical nwxc_has_cam
      external nwxc_has_cam
      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_moocc))
     &     call errquit('WFAT: Deallocating occnum failed', 3, MA_ERR)
      if (.not. ma_free_heap(l_moerg))
     &     call errquit('WFAT: Deallocating erg failed', 3, MA_ERR)
      if (.not. ma_free_heap(l_modgid))
     &     call errquit('WFAT: Deallocating dgid failed', 3, MA_ERR)
      if (.not. ma_free_heap(l_modgn))
     &     call errquit('WFAT: Deallocating dgn failed', 3, MA_ERR)
      if (.not. ma_free_heap(l_norbdg))
     &     call errquit('WFAT: Deallocating norb_dg failed', 3, MA_ERR)

      
c     ==== Geometry ====
      if (.not. geom_destroy(geom_hd))
     &     call errquit('WFAT: Error in destroying geometry object.',
     &     1, GEOM_ERR)
      if (nwxc_has_cam()) then
         if (.not. geom_destroy(rsgeom)) call errquit('WFAT: Error'//
     &        ' in destroying geometry object for the RS functional.',
     &        1, GEOM_ERR)
      endif

      
c     ==== Basis ====
      if (.not. bas_destroy(bas_hd))
     &     call errquit('WFAT: Error in destroying basis object.',
     &     1, BASIS_ERR)
      if (nwxc_has_cam()) then
         if (.not. bas_destroy(rsbas)) call errquit('WFAT: Error'//
     &        ' in destroying basis object for the RS functional.',
     &        1, BASIS_ERR)
      endif

      
c     ==== Molecular orbitals ====
      do is = 1, nspin0
         if (is == 1) then
            spin_name = '_a'
         else if (is == 2) then
            spin_name = '_b'
         endif
   
         mov_name = trim(adjustl(mov_name0))//spin_name
         if (.not. ga_destroy(g_mov(is)))
     $        call errquit ("WFAT: couldnt destroy "//mov_name,
     &        0, GA_ERR)
      enddo

      
c     ==== Molecular integrals ====
      call int_terminate()
      call schwarz_tidy()


c     ==== Common integrals ====
      if (.not. ga_destroy(gwf_ovl))
     $     call errquit ("WFAT: couldnt destroy overlap", 0, GA_ERR)
      if (.not. ga_destroy(gwf_kin))
     $     call errquit ("WFAT: couldnt destroy kinetic", 0, GA_ERR)
      if (.not. ga_destroy(gwf_pot))
     $     call errquit ("WFAT: couldnt destroy potential", 0, GA_ERR)
      

c     ==== RTDB ====
      ll = rtdb_delete(rtdb, 'wfat:print_srfc')
      ll = rtdb_delete(rtdb, 'wfat:print_rate')
      ll = rtdb_delete(rtdb, 'wfat:print_yild')
      ll = rtdb_delete(rtdb, 'wfat:print_tmd')
      ll = rtdb_delete(rtdb, 'wfat:print_tmd_atom')
      ll = rtdb_delete(rtdb, 'wfat:beta_start')
      ll = rtdb_delete(rtdb, 'wfat:beta_last')
      ll = rtdb_delete(rtdb, 'wfat:nbeta')
      ll = rtdb_delete(rtdb, 'wfat:gamma_start')
      ll = rtdb_delete(rtdb, 'wfat:gamma_last')
      ll = rtdb_delete(rtdb, 'wfat:ngamma')
      ll = rtdb_delete(rtdb, 'wfat:nktmd')
      ll = rtdb_delete(rtdb, 'wfat:beta_tmd')
      ll = rtdb_delete(rtdb, 'wfat:gamma_tmd')
      ll = rtdb_delete(rtdb, 'wfat:intgl_name')
      ll = rtdb_delete(rtdb, 'wfat:intgl_nosave')
      ll = rtdb_delete(rtdb, 'wfat:movec_name')

      
     
      
c     ==== Miscellaneous ====
cOLD      call wfat_set_fockreplicate(0)

      
      end
c=======================================================
      
