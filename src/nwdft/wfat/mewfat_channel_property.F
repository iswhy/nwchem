c=================================================
c=================================================
      subroutine mewfat_channel_property(rtdb, geom, bas, nb, gq, mcel, 
     &     n0, n1, g_mo0, g_mo1, chg0, chg1, mult0, mult1, gserg0, 
     &     gserg1, ionerg, edipm0, edipm1, chcent, chmmt, dipm1)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "mewfat_qm_param.fh"
      integer rtdb                  ! [input]
      integer geom                  ! [input]
      integer bas                   ! [input]
      integer nb                    ! [input]
      integer gq
      integer mcel
      integer n0(2)                 ! [input]
      integer n1(2)                 ! [input]
      integer g_mo0(2)              ! [input]
      integer g_mo1(2)              ! [input]
      double precision chg0         ! [input]
      double precision chg1         ! [input]
      integer mult0                 ! [input]
      integer mult1                 ! [input]
      double precision gserg0       ! [output], neutral GS energy
      double precision gserg1       ! [output], cation GS energy
      double precision ionerg       ! [output], ionization energy
      double precision edipm0(3)    ! [output], unshifted neutral e. dipole moment
      double precision edipm1(3)    ! [output], unshifted cation e. dipole moment
      double precision chcent(3)    ! [output], channel center.
      double precision chmmt(3)     ! [output], shifted channel moment.
      double precision dipm1(3)     ! [output], shifted cation dipole moment.

      integer l_mo0, k_mo0, l_mo1, k_mo1, dummy_int0, dummy_int1
      integer nocmo0, nocmo1, is, mo0, i_mo0, mo1, i_mo1
      double precision manip
      

c     ==== Ionization energy ====
      call mewfat_channel_ionize(rtdb, geom, bas, nb, gq, mcel, n0,
     &     n1, g_mo0, g_mo1, chg0, chg1, mult0, mult1, tol2e, oskel, 
     &     ixcfit, iolgc, gserg0, gserg1, ionerg)
      call wfat_manerg(rtdb, 1, ionerg)

      if (rtdb_get(rtdb, 'mewfat:manual_ip', mt_dbl, 1, manip))
     &     then
         ionerg = manip
         if (ga_nodeid()==0) write(luout, '(3x, "Manual IP is used.")')
      endif


c     ==== Center, 1st moment, and cation dipole ====
      nocmo0 = nb * sum(n0)
      if (.not. ma_alloc_get(mt_dbl, nocmo0, 'MO0', l_mo0, k_mo0))
     &     call errquit('MEWFAT: Cannot allocate MO0', 1, MA_ERR)
      nocmo1 = nb * sum(n1)
      if (.not. ma_alloc_get(mt_dbl, nocmo1, 'MO1', l_mo1, k_mo1))
     &     call errquit('MEWFAT: Cannot allocate MO1', 1, MA_ERR)
      do is = 1, 2
         mo0 = (is-1)*n0(1) + 1
         i_mo0 = k_mo0 + (mo0-1)*nb
         if (n0(is) > 0)
     &        call ga_get(g_mo0(is), 1, nb, 1, n0(is), dbl_mb(i_mo0),
     &        nb)
         mo1 = (is-1)*n1(1) + 1
         i_mo1 = k_mo1 + (mo1-1)*nb
         if (n1(is) > 0)
     &        call ga_get(g_mo1(is), 1, nb, 1, n1(is), dbl_mb(i_mo1),
     &        nb)
      enddo
      call mewfat_channel_moment(rtdb, geom, bas, nb, n0, n1, 
     &     dbl_mb(k_mo0), dbl_mb(k_mo1), edipm0, edipm1, chcent, chmmt,
     &     dipm1)


c     ==== Print Channel properties ====
      if (ga_nodeid() == 0) then
         dummy_int0 = 0
         dummy_int1 = 0
         call mewfat_print_channel_property(1, dummy_int0, dummy_int1, 
     &        spinx, gserg0, gserg1, ionerg, edipm0, edipm1, chcent,
     &        chmmt, dipm1)
      endif      

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_mo0))
     &     call errquit('MEWFAT: Cannot deallocate MO0', 1, MA_ERR)
      if (.not. ma_free_heap(l_mo1))
     &     call errquit('MEWFAT: Cannot deallocate MO1', 1, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_print_channel_property(nch, state0, state1, 
     &     spinid, gserg0, gserg1, ionerg, edipm0, edipm1, chcent, 
     &     chmmt, iondip)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      integer nch
      integer state0(nch)
      integer state1(nch)
      integer spinid(nch)
      integer gserg0(nch)
      integer gserg1(nch)
      double precision ionerg(nch)
      double precision edipm0(3,nch)    
      double precision edipm1(3,nch)
      double precision chcent(3,nch)
      double precision chmmt(3,nch)
      double precision iondip(3,nch)

      character*5 ss
      integer i, io0, io1
      

      write(luout, '(2/)')
      call util_print_centered(luout,
     &     'Properties of the user-requested ionization channels',
     &     40, .true.)
      write(luout, '(" ")')
      
      do i = 1, nch
         io0 = state0(i)
         io1 = state1(i)
         if (spinid(i) == 1) then
            ss = 'alpha'
         else if (spinid(i) == 2) then
            ss = 'beta'
         endif
         
         write(luout, '(5x, "Channel : ", i6)') i
         write(luout, '(7x, a15, " = ", a)')
     &        "* ionized spin ", trim(adjustl(ss))
         write(luout, '(7x, a26, " = ", i3, " -->", i3)')
     &        "* neutral --> cation states ", io0, io1
         write(luout, '(7x, a19, " = ", f18.12, " Hartree")')
     &        "* neutral GS energy", gserg0(i)
         write(luout, '(7x, a19, " = ", f18.12, " Hartree")')
     &        "* cation GS energy ", gserg1(i)
         write(luout, '(7x, a19, " = ", f18.12, " Hartree")')
     &        "* ionization energy", ionerg(i)
         write(luout, '(7x, a49, " = ", 3es14.5)')
     &        "* unshifted neutral e. dipole moment (xyz, Bohr) ",
     &        edipm0(:, i)
         write(luout, '(7x, a49, " = ", 3es14.5)')
     &        "* unshifted cation e. dipole moment (xyz, Bohr)  ",
     &        edipm1(:, i)
         write(luout, '(7x, a38, " = ", 3es14.5)')
     &        "* channel origin (xyz, Bohr)          ", chcent(:, i)
         write(luout, '(7x, a38, " = ", 3es14.5)') 
     &        "* unshifted channel moment (xyz, Bohr)",
     &        chmmt(:, i) + chcent(:, i)
         write(luout, '(7x, a38, " = ", 3es14.5)') 
     &        "* shifted channel moment (xyz, Bohr)  ", chmmt(:, i)
         write(luout, '(7x, a38, " = ", 3es14.5)')
     &        "* shifted ion dipole (xyz, e.Bohr)    ", iondip(:, i)
      enddo

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_channel_ionize(rtdb, geom, bas, nb, gq, mcel, 
     &     n0, n1, g_mo0, g_mo1, chg0, chg1, mult0, mult1, tol2e, oskel,
     &     ixcfit, iolgc, gserg0, gserg1, ionerg)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb              ! [input]
      integer geom              ! [input]
      integer bas               ! [input]
      integer nb                ! [input]
      integer gq
      integer mcel
      integer n0(2)             ! [input]
      integer n1(2)             ! [input]
      integer g_mo0(2)          ! [input]
      integer g_mo1(2)          ! [input]
      double precision chg0     ! [input]
      double precision chg1     ! [input]
      integer mult0             ! [input]
      integer mult1             ! [input]
      double precision tol2e    ! [input]
      logical oskel             ! [input]
      integer ixcfit            ! [input]
      logical iolgc             ! [input]
      double precision gserg0   ! [output], neutral GS energy
      double precision gserg1   ! [output], cation GS energy
      double precision ionerg   ! [output], ionization energy

      integer ndid, n1e, nfk, nocmo0, nocmo1, nrd, ggg, jjj
      integer l_1e, l_fk0, l_fk1, l_mo0, l_mo1, l_rd0, l_rd1
      integer k_1e, k_fk0, k_fk1, k_mo0, k_mo1, k_rd0, k_rd1
      integer is, i_fk, mo0, i_mo0, mo1, i_mo1
      integer g_1e, g_fk0(2), g_fk1(2)
      
      integer ga_create_atom_blocked
      external ga_create_atom_blocked


      ndid = ga_nodeid()
      
      
c     ==== 1e Hamiltonian ====
      g_1e = ga_create_atom_blocked(geom, bas, '1e_opr')
      call ga_zero(g_1e)
      call int_1e_ga(bas, bas, g_1e, 'kinetic', oskel)
      call int_1e_ga(bas, bas, g_1e, 'potential', oskel)
      n1e = nb * nb
      if (.not. ma_alloc_get(mt_dbl, n1e, '1e_opr', l_1e, k_1e))
     &     call errquit('MEWFAT: Cannot allocate 1e_opr', 1, MA_ERR)
      call ga_get(g_1e, 1, nb, 1, nb, dbl_mb(k_1e), nb)

      
c     ==== Fock matrices ====
      do is = 1, 2
         g_fk0(is) = ga_create_atom_blocked(geom, bas, 'fock0')
         g_fk1(is) = ga_create_atom_blocked(geom, bas, 'fock1')
      enddo
      call wfat_get_fock(rtdb, geom, bas, nb, 2, n0, g_mo0,
     &     chg0, mult0, tol2e, oskel, ixcfit, iolgc, .false., g_fk0)
      call wfat_get_fock(rtdb, geom, bas, nb, 2, n1, g_mo1,
     &     chg1, mult1, tol2e, oskel, ixcfit, iolgc, .false., g_fk1)
      nfk = nb * nb * 2
      if (.not. ma_alloc_get(mt_dbl, nfk, 'fock0', l_fk0, k_fk0))
     &     call errquit('MEWFAT: Cannot allocate fock0', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nfk, 'fock1', l_fk1, k_fk1))
     &     call errquit('MEWFAT: Cannot allocate fock1', 1, MA_ERR)
      do is = 1, 2
         i_fk = k_fk0 + (is-1)*nb*nb
         call ga_get(g_fk0(is), 1, nb, 1, nb, dbl_mb(i_fk), nb)
         i_fk = k_fk1 + (is-1)*nb*nb
         call ga_get(g_fk1(is), 1, nb, 1, nb, dbl_mb(i_fk), nb)
      enddo


c     ==== Calculate RDM's ====
      nocmo0 = nb * sum(n0)
      if (.not. ma_alloc_get(mt_dbl, nocmo0, 'MO0', l_mo0, k_mo0))
     &     call errquit('MEWFAT: Cannot allocate MO0', 1, MA_ERR)
      nocmo1 = nb * sum(n1)
      if (.not. ma_alloc_get(mt_dbl, nocmo1, 'MO1', l_mo1, k_mo1))
     &     call errquit('MEWFAT: Cannot allocate MO1', 1, MA_ERR)

      do is = 1, 2
         mo0 = (is-1)*n0(1) + 1
         i_mo0 = k_mo0 + (mo0-1)*nb
         if (n0(is) > 0)
     &        call ga_get(g_mo0(is), 1, nb, 1, n0(is), dbl_mb(i_mo0),
     &        nb)
         mo1 = (is-1)*n1(1) + 1
         i_mo1 = k_mo1 + (mo1-1)*nb
         if (n1(is) > 0)
     &        call ga_get(g_mo1(is), 1, nb, 1, n1(is), dbl_mb(i_mo1),
     &        nb)
      enddo

      nrd = nb * nb * 2
      if (.not. ma_alloc_get(mt_dbl, nrd, 'rdm0', l_rd0, k_rd0))
     &     call errquit('MEWFAT: Cannot allocate rdm0', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nrd, 'rdm1', l_rd1, k_rd1))
     &     call errquit('MEWFAT: Cannot allocate rdm1', 1, MA_ERR)
      call mewfat_calc_rdm(nb, n0, n1, dbl_mb(k_mo0), dbl_mb(k_mo1),
     &     dbl_mb(k_rd0), dbl_mb(k_rd1))


c     ==== GS and ionization energies ====
      call mewfat_calc_ionerg(geom, bas, nb, n0, n1, gq, mcel, 
     &     dbl_mb(k_1e), dbl_mb(k_rd0), dbl_mb(k_rd1), dbl_mb(k_fk0), 
     &     dbl_mb(k_fk1), gserg0, gserg1, ionerg)


c     ==== Deallocations ====
      if (.not. ga_destroy(g_1e))
     &     call errquit('MEWFAT: Cannot destroy 1e_opr', 18, GA_ERR)
      do is = 1, 2
         if (.not. ga_destroy(g_fk0(is)))
     &        call errquit('MEWFAT: Cannot destroy fock0', 13, GA_ERR)
         if (.not. ga_destroy(g_fk1(is)))
     &        call errquit('MEWFAT: Cannot destroy fock1', 13, GA_ERR)
      enddo
      if (.not. ma_free_heap(l_1e))
     &     call errquit('MEWFAT: Cannot deallocate 1e_opr', 1, MA_ERR)
      if (.not. ma_free_heap(l_fk0))
     &     call errquit('MEWFAT: Cannot deallocate fock0', 1, MA_ERR)
      if (.not. ma_free_heap(l_fk1))
     &     call errquit('MEWFAT: Cannot deallocate fock1', 1, MA_ERR)
      if (.not. ma_free_heap(l_mo0))
     &     call errquit('MEWFAT: Cannot deallocate MO0', 1, MA_ERR)
      if (.not. ma_free_heap(l_mo1))
     &     call errquit('MEWFAT: Cannot deallocate MO1', 1, MA_ERR)
      if (.not. ma_free_heap(l_rd0))
     &     call errquit('MEWFAT: Cannot deallocate rdm0', 5, MA_ERR)
      if (.not. ma_free_heap(l_rd1))
     &     call errquit('MEWFAT: Cannot deallocate rdm1', 1, MA_ERR)
      
      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_calc_ionerg(geom, bas, nb, n0, n1, gq, mcel, 
     &           h1e, rdm0, rdm1, fock0, fock1, gse0, gse1, ionerg)

      implicit none
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
      integer geom
      integer bas
      integer nb                ! [input]
      integer n0(2)
      integer n1(2)
      integer gq
      integer mcel
      double precision h1e(nb,nb)      ! [input]
      double precision rdm0(nb,nb,2)   ! [input]
      double precision rdm1(nb,nb,2)   ! [input]
      double precision fock0(nb,nb,2)  ! [input]
      double precision fock1(nb,nb,2)  ! [input]
      double precision gse0            ! [output]
      double precision gse1            ! [output]
      double precision ionerg          ! [output]

      logical ll
      integer ndid, i, j, ns0, ns1
      double precision g, xce0, xce1, enuc


      ndid = ga_nodeid()

      
c     ==== Calculate the DFT-XC part of total electronic energy ====
      if (n0(1)==n0(2)) then
         ns0 = 1
      else
         ns0 = 2
      endif
      call mewfat_dfterg(geom, bas, ns0, gq, mcel, rdm0, xce0)
      
      if (n1(1)==n1(2)) then
         ns1 = 1
      else
         ns1 = 2
      endif
      call mewfat_dfterg(geom, bas, ns1, gq, mcel, rdm1, xce1)
      

c     ==== Calculate the HF part of total electronic energy ====
      gse0 = 0.0d0
      gse1 = 0.0d0
      do j = 1, nb
         do i = 1, nb
            g = (rdm0(i,j,1) + rdm0(i,j,2)) * h1e(j,i) +
     &           rdm0(i,j,1) * fock0(j,i,1) +
     &           rdm0(i,j,2) * fock0(j,i,2)
            gse0 = gse0 + g

            g = (rdm1(i,j,1) + rdm1(i,j,2)) * h1e(j,i) +
     &           rdm1(i,j,1) * fock1(j,i,1) +
     &           rdm1(i,j,2) * fock1(j,i,2)
            gse1 = gse1 + g
         enddo
      enddo

      
c     ==== Calculate the total energy ====
      ll = geom_nuc_rep_energy(geom, enuc)
      gse0 = gse0 / 2 + xce0 + enuc
      gse1 = gse1 / 2 + xce1 + enuc
      ionerg = gse0 - gse1
      
      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_channel_moment(rtdb, geom, bas, nb, n0, n1, 
     &     ocmo0, ocmo1, edipm0, edipm1, chcent, chmmt, dipm1)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      integer geom                               ! [input],
      integer bas                                ! [input],
      integer nb                                 ! [input],
      integer n0(2)                              ! [input],
      integer n1(2)                              ! [input],
      double precision ocmo0(nb, n0(1)+n0(2))    ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))    ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision edipm0(3)                 ! [output], unshifted neutral e. dipole moment.
      double precision edipm1(3)                 ! [output], unshifted cation e. dipole moment.
      double precision chcent(3)                 ! [output], channel center.
      double precision chmmt(3)                  ! [output], shifted channel moment.
      double precision dipm1(3)                  ! [output], shifted cation dipole moment.

      logical umandip(1)
      double precision mdip(3)

      
      call mewfat_get_orbcent(rtdb, geom, bas, nb, n0, n1, ocmo0, ocmo1,
     &     edipm0, edipm1, chcent, dipm1)


c     ==== When using manual dipole ====
      if (rtdb_get(rtdb, 'wfat:use_mandip', mt_log, 1, umandip)) then
         if (umandip(1)) then
            if (.not. rtdb_get(rtdb, 'wfat:mandip', mt_dbl, 3*1, mdip))
     &           call errquit('MEWFAT: Manual dipole is used but '//
     &           'cannot find it in RTDB.', 1, RTDB_ERR)
c           mdip contains alternatives to
c              -(<Psi|\vec{r}|Psi> - <Psi^+|\vec{r}|Psi^+>)
         endif
      else
         umandip(1) = .false.
      endif

      
      if (umandip(1)) then
         chmmt = mdip - chcent               ! Shifted channel first moment
      else
         chmmt = -(edipm0 - edipm1) - chcent ! Shifted channel first moment
      endif
      
      end
c=================================================
      
      
c=================================================
c=================================================
      subroutine mewfat_get_orbcent(rtdb, geom, bas, nb, n0, n1, ocmo0,
     &     ocmo1, edipm0, edipm1, chcent, dipm1)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
      integer rtdb
      integer geom                               ! [input],
      integer bas                                ! [input],
      integer nb                                 ! [input],
      integer n0(2)                              ! [input],
      integer n1(2)                              ! [input],
      double precision ocmo0(nb, n0(1)+n0(2))    ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))    ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision edipm0(3)                 ! [output], unshifted neutral e. dipole moment.
      double precision edipm1(3)                 ! [output], unshifted cation e. dipole moment.
      double precision chcent(3)                 ! [output], channel center.
cOLD      double precision chmmt(3)                  ! [output], shifted channel moment.
      double precision dipm1(3)                  ! [output], shifted cation dipole moment.

      character*8 dipname(3)
      logical umanpos(1)
      integer NMMT, i, g_dp(3), i_dp, l_dp, k_dp, is, io, mo
      parameter (NMMT=4)
      double precision zt, nucmmt(NMMT), ndipm1(3), manpos(3,1)
      double precision moment(3)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      
      
      dipname = (/ 'dipole_x', 'dipole_y', 'dipole_z' /) 

      
c     ==== Calculating dipole operators ====
      do i = 1, 3
         g_dp(i) = ga_create_atom_blocked(geom, bas, dipname(i))
         call ga_zero(g_dp(i))
      enddo
      call int_dip_ga(bas, bas, g_dp(1), g_dp(2), g_dp(3))
      if (.not. ma_alloc_get(mt_dbl, 3*nb*nb, 'dip_opr', l_dp, k_dp))
     &     call errquit('MEWFAT: Cannot allocate dip_opr', 1, MA_ERR)
      do i = 1,3
         i_dp = k_dp + (i-1)*nb*nb
         call ga_get(g_dp(i), 1, nb, 1, nb, dbl_mb(i_dp), nb)
      enddo
      

c     ==== Calculating shifted channel first moment ====
cc      st = rtdb_get(rtdb, 'dft:dipole', mt_dbl, 3, edipm0)     ! This is for getting dipole moment from RTDB.
      edipm0 = 0.0d0
      edipm1 = 0.0d0
      do is = 1, 2
c        ==== Neutral dipole ====
         do io = 1, n0(is)
            mo = (is-1)*n0(1) + io
            do i = 1, 3
               i_dp = k_dp + (i-1)*nb*nb
               call wfat_1e_expval_1eopr(nb, ocmo0(:,mo),
     &              dbl_mb(i_dp), moment(i))
            enddo
            edipm0 = edipm0 - moment
         enddo
         
c        ==== Cation dipole ====
         do io = 1, n1(is)
            mo = (is-1)*n1(1) + io
            do i = 1, 3
               i_dp = k_dp + (i-1)*nb*nb
               call wfat_1e_expval_1eopr(nb, ocmo1(:,mo),
     &              dbl_mb(i_dp), moment(i))
            enddo
            edipm1 = edipm1 - moment
         enddo
      enddo

      
      if (.not. geom_nuc_charge(geom, zt))
     &     call errquit('MEWFAT: Getting zt failed', 5, UERR)
      call geom_nuc_mpole(geom, (/ 0.0d0, 0.0d0, 0.0d0 /), 1,
     &     nucmmt, NMMT)
      ndipm1 = nucmmt(2:4)

      
      call wfat_man_chpos(rtdb, 1, umanpos, manpos)
      if (umanpos(1)) then
         chcent = manpos(:,1)
      else
         chcent = (ndipm1 + edipm1) / (zt-sum(n1)) ! Channel center
      endif
cOLD      chmmt = -(edipm0 - edipm1) - chcent            ! Shifted channel first moment
      dipm1 = ndipm1 + edipm1 - (zt-sum(n1))*chcent  ! Shifted cation dipole moment, should be zero.


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dp))
     &     call errquit('MEWFAT: Cannot deallocate dip_opr', 1, MA_ERR)
      do i = 1, 3
         if (.not. ga_destroy(g_dp(i)))
     &        call errquit('MEWFAT: Cannot deallocate dipole',
     &        1, GA_ERR)
      enddo
         
      end
c=================================================
      
