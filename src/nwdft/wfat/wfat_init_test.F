c=======================================================
c=======================================================
      subroutine wfat_init_test()

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "errquit.fh"
#include "stdio.fh"      
#include "wfat.fh"      
#include "wfat_qm_param_test.fh"


      integer is, ioc, ieg
      integer k_oc, k_eg
      integer l_oc, l_eg, i

      
      
      call wfat_init_getcdft_test(bas_hd, nspin0, geom_hd, natm)
      call wfat_init_getcscf_test(zion, nel, g_mov, mov_in)

      if (nspin0 == 1) then
         nelec0 = 2 * nel(1)
      else if (nspin0 == 2) then
         nelec0 = sum(nel)
      endif

      zion = zion + 1


c      call wfat_init_print()
      
      
      WF_SF_PRINT_READY = .false.
      WF_IY_PRINT_READY = .false.

      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_init_getcdft_test(bas_hd, nspin, geom_hd, natm)

      implicit none
#include "cdft.fh"
#include "geom.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"      
      integer bas_hd
      integer nspin
      integer geom_hd
      integer natm

      character*20 basis_name, trans_name
      integer basnormid

      
      

      
      bas_hd = ao_bas_han
      geom_hd = geom
      nspin = ipol


      if(.not. bas_name(bas_hd, basis_name, trans_name))
     &     call errquit('WFAT: Error in getting basis name.',
     &     1, BASIS_ERR)
      write(luout, '("basis_name = ", a20)') basis_name
      write(luout, '("trans_name = ", a20)') trans_name

      if (.not. bas_norm_get(bas_hd, basnormid))
     &     call errquit('WFAT: Error getting basis norm ID',
     &     1, BASIS_ERR)
      write(luout, '("Basis norm ID test = ", i3)') basnormid



      

      if (.not. geom_check_handle(geom, 'not valid in init'))
     &     call errquit('WFAT: Getting geom failed', 1, GEOM_ERR)
      
      write(luout, '("geom handle in init = ", i4)') geom
      if (.not. geom_ncent(geom, natm))
     &     call errquit('WFAT: Getting natm failed', 51, GEOM_ERR)


      
      end
c=======================================================            


c=======================================================
c=======================================================
      subroutine wfat_init_getcscf_test(ch, nel, g_mov, mov_in)

      implicit none
#include "cscf.fh"
      double precision ch
      integer nel(2)
      integer g_mov(2)
      character*256 mov_in
      

      ch = charge
      nel(1) = nalpha
      nel(2) = nbeta
      g_mov = g_movecs
      mov_in = movecs_in
c      mov_in = movecs_out

      
      end
c=======================================================


