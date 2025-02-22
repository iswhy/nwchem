c=======================================================
c=======================================================
      subroutine wfat_orb_property(rtdb, nch, orbid, spinid, 
     &           vcf, erg, orbcent, orbmmt, iondip)

      implicit none
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "wfat_qm_param.fh"
#include "global.fh"
      integer rtdb
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision vcf(nbas0, nch)
      double precision erg(nch)
      double precision orbcent(3, nch)
      double precision orbmmt(3, nch)
      double precision iondip(3, nch)

      logical ismandip
      
            
      call wfat_get_coef(nch, orbid, spinid, vcf, erg)
      call wfat_manerg(rtdb, nch, erg)
      call wfat_channel_moment(rtdb, nch, orbid, spinid, orbcent, 
     &     orbmmt, iondip)
      
      if (ga_nodeid() == 0) then
         call wfat_print_orb_property(nch, orbid, spinid, erg,
     &        orbcent, orbmmt, iondip)
      endif      
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_print_orb_property(nch, orbid, spinid, erg,
     &           orbcent, orbmmt, iondip)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision erg(nch)
      double precision orbcent(3,nch)
      double precision orbmmt(3,nch)
      double precision iondip(3,nch)

      character*5 ss
      integer i, io, is
      

      write(luout, '(2/)')
      call util_print_centered(luout,
     &     'Properties of the user-requested ionization channels',
     &     40, .true.)
      
      do i = 1, nch
         io = orbid(i)
         is = spinid(i)
         if (is == 1) then
            ss = 'alpha'
         else if (is == 2) then
            ss = 'beta'
         endif

         write(luout, '("")')
         write(luout, '(3x, "Channel : ", i6)') i
         write(luout, '(5x, "* ", a13, " = ", a)')
     &        "ionized spin ", trim(adjustl(ss))
         write(luout, '(5x, "* ", a13, " = ", i4)')
     &        "orbital no.  ", io
         write(luout, '(5x, "* ", a17, " = ", es14.5, " Hartree")')
     &        "ionization energy", erg(i)
         write(luout, '(5x, "* ", a36, " = ", 3es14.5)')
     &        "channel position (xyz, Bohr)        ", orbcent(:, i)
         write(luout, '(5x, "* ", a36, " = ", 3es14.5)') 
     &        "unshifted channel moment (xyz, Bohr)",
     &        orbmmt(:, i) + orbcent(:, i)
         write(luout, '(5x, "* ", a36, " = ", 3es14.5)') 
     &        "shifted channel moment (xyz, Bohr)  ", orbmmt(:, i)
         write(luout, '(5x, "* ", a36, " = ", 3es14.5)')
     &        "shifted ion dipole (xyz, e.Bohr)    ", iondip(:, i)
      enddo
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_get_coef(nch, orbid, spinid, vcf, erg)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "wfat_qm_param.fh"      
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision vcf(nbas0,nch)
      double precision erg(nch)

      integer is, io, i, ieg


      do i = 1, nch
         io = orbid(i)
         is = spinid(i)

         if (is > nspin0)
     &        call errquit('WFAT: An element of argument spinid of '//
     &        'wfat_get_coef exceeds the number of spin channels.',
     &        1, UERR)
         call ga_get(g_mov(is), 1, nbas0, io, io, vcf(1,i), nbas0)

         ieg = k_moerg + (is-1)*nbas0 + io - 1
         erg(i) = dbl_mb(ieg)
      enddo

      end
c     Notes:
c       1) Note that this call to movecs_read(...) assumes that the basis set and thus
c          the number of basis functions is equal to nbas0 (argument of this routine).
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_get_channel_coef(nch, orbid, spinid, chcf)
c This routine is OLD.
      implicit none
#include "errquit.fh"
#include "stdio.fh"      
#include "global.fh"
c     #include "cscf.fh"
#include "wfat_qm_param.fh"
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision chcf(nbas0, nch)

      integer ndid
      integer i, is, io, it, nb, nn


      call ga_inquire(g_mov(1), it, nb, nn)
      if (nbas0 /= nb) then
         call errquit('WFAT: Argument nbas0 of wfat_get_channel_coef '//
     &        'must be equal to the number of basis.', 1, UERR)
      endif
      if (maxval(orbid) > nn .or. minval(orbid) < 1) then
         write(luout, '("In wfat_get_channel_coef, norb = ", i5)') nn
         call errquit('WFAT: The elements of argument array orbid '//
     &        'must be positive and less than norb (see above).',
     &        2, UERR)
      endif
      
      
      do i = 1, nch
         is = spinid(i)
         io = orbid(i)
         call ga_get(g_mov(is), 1, nbas0, io, io, chcf(1,i), nbas0)
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_channel_moment(rtdb, nch, orbid, spinid, 
     &           orbcent, orbmmt, iondip)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision orbcent(3, nch)   ! [output], channel center (shifting amount).
      double precision orbmmt(3, nch)    ! [output], shifted channel dipole moment.
      double precision iondip(3, nch)    ! [output], shifted total dipole of cation.

      logical umandip(nch)
      integer l_om0, k_om0, l_omm, k_omm, iom0, iomm
      integer i, id


      if (.not. MA_alloc_get(MT_DBL, 3*nch, 'chmmt0', l_om0, k_om0))
     &     call errquit('WFAT: chmmt0 allocation failed', 1, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, 3*nch, 'chmmt0', l_omm, k_omm))
     &     call errquit('WFAT: chmmt_man allocation failed', 1, MA_ERR)
      call wfat_get_orbcent(rtdb, nch, orbid, spinid, orbcent,
     &     dbl_mb(k_om0), iondip)

      
c     ==== When using manual dipole ====
      if (rtdb_get(rtdb, 'wfat:use_mandip', mt_log, nch, umandip)) then
         do i = 1, nch
            if (umandip(i)) then
               if (rtdb_get(rtdb, 'wfat:mandip', mt_dbl, 3*nch,
     &              dbl_mb(k_omm))) then
c                 dbl_mb(k_omm) contains alternatives to <psi|\vec{r}|psi>,
c                 note that there is no negative sign associated with electron charge.
                  exit
               else
                  call errquit('WFAT: Manual dipole is used but '//
     &                 'cannot find it in RTDB.', 1, RTDB_ERR)
               endif
            endif
         enddo
      else
         umandip(:) = .false.
      endif

      
c     ==== Compute orb. moments relative to the orb. center ====
      do i = 1, nch
         do id = 1, 3
            if (umandip(i)) then
               iomm = k_omm + 3*(i-1) + (id-1)
               orbmmt(id,i) = dbl_mb(iomm) - orbcent(id,i)
            else
               iom0 = k_om0 + 3*(i-1) + (id-1)
               orbmmt(id,i) = dbl_mb(iom0) - orbcent(id,i)
            endif
         enddo
      enddo

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_om0))
     &     call errquit('WFAT: chmmt0 deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_omm))
     &     call errquit('WFAT: chmmt_man deallocation failed', 1,
     &     MA_ERR)
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_get_orbcent(rtdb, nch, orbid, spinid, 
     &           orbcent, orbmmt0, iondip)

      implicit none
#include "stdio.fh"
#include "geom.fh"
#include "global.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "wfat_qm_param.fh"
      integer rtdb
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision orbcent(3, nch)   ! [output], channel center (shifting amount).
      double precision orbmmt0(3, nch)   ! [output], unshifted channel dipole moment, <psi|\vec{r}|psi>. Note that there is no negative sign.
      double precision iondip(3, nch)    ! [output], shifted total dipole of cation.

      logical umanpos(nch)
      integer ndid, nmp, NMMT, ncf, nom0, niondip
      integer l_mp, l_cf, l_om0
      integer k_mp, k_cf, k_om0
      integer is, icf, id, iq, imp, iom0, i, io, k
      parameter (NMMT = 4)
      double precision vorig(3), ndip0(NMMT), edip0(3), totdip0(3)
      double precision totzn, manpos(3,nch), ndip(3), edip(3), totdip(3)


      ndid = ga_nodeid()
      vorig = 0.0d0

c     ==== Nuclear moment relative to origin ====
      call geom_nuc_mpole(geom_hd, vorig, 1, ndip0, NMMT)

c     ==== Dipole matrices relative to origin ====
      nmp = nbas0 * nbas0 * NMMT     !! 2)
      if (.not. MA_alloc_get(MT_DBL, nmp, 'mpmat', l_mp, k_mp))
     &     call errquit('WFAT: mpmat allocation failed', 1, MA_ERR)
      call int_mpole_all(bas_hd, nbas0, 1, vorig, dbl_mb(k_mp))
      

c     ==== Get orb. coeffs to calculate electronic moments ====
      ncf = nbas0 * nbas0 * nspin0
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vcf', l_cf, k_cf))
     &     call errquit('WFAT: vcf allocation failed', 1, MA_ERR)
      do is = 1, nspin0
         do io = 1, nbas0
            icf = k_cf + (is-1)*nbas0*nbas0 + (io-1)*nbas0
            call ga_get(g_mov(is), 1, nbas0, io, io,
     &                  dbl_mb(icf), nbas0)
         enddo
      enddo

c     ==== Electronic moments relative to origin ====
      nom0 = nbas0 * nspin0 * 3     !! 1)
      if (.not. MA_alloc_get(MT_DBL, nom0, 'orbmmt0', l_om0, k_om0))
     &     call errquit('WFAT: orbmmt0 allocation failed', 1, MA_ERR)

      do id = 1, 3
         iq = id + 1
         imp = k_mp + (iq-1)*nbas0*nbas0
         iom0 = k_om0 + (id-1)*nbas0*nspin0
         call wfat_ne_expval_1eopr(nbas0, nspin0, nel, dbl_mb(k_cf),
     &        dbl_mb(imp), dbl_mb(iom0), edip0(id))
      enddo
      
c     ==== Total dipole relative to origin ====
      totdip0 = ndip0(NMMT-2:NMMT) - edip0


c     ==== Orb. center and orb. moments relative to orb. center ====
      call wfat_man_chpos(rtdb, nch, umanpos, manpos)
      if (.not. geom_nuc_charge(geom_hd, totzn))
     &     call errquit('WFAT: Getting totzn failed', 5, UERR)
      do i = 1, nch
         io = orbid(i)
         is = spinid(i)
         do id = 1, 3
c           ==== Compute orb. center ====
            iom0 = k_om0 + nbas0*nspin0*(id-1) + nbas0*(is-1) + (io-1)
            if (umanpos(i)) then
               orbcent(id,i) = manpos(id,i)
            else
               orbcent(id,i) = (totdip0(id) + dbl_mb(iom0)) /
     &                         (totzn - nelec0 + 1)
            endif

c           ==== Get orb. moments relative to origin ====
            orbmmt0(id,i) = dbl_mb(iom0)
         enddo
         
c        ==== Shifted ion dip. moment (should be vanishing) ====
         ndip = ndip0(2:4) - totzn * orbcent(:,i)
         edip = edip0(1:3) - nelec0 * orbcent(:,i)
         totdip = ndip - edip
         iondip(:,i) = totdip + (orbmmt0(:,i) - orbcent(:,i))
      enddo

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_mp))
     &     call errquit('WFAT: mpmat deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: vcf deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_om0))
     &     call errquit('WFAT: orbmmt0 deallocation failed', 1, MA_ERR)

      end
c     Notes:
c       1) orbmmt0 is conceptually a 3D array of size [nbas0, nspin0, 3].
c       2) mpmat is conceptually a 3D array of size [nbas0, nbas0, NMMT].
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_manerg(rtdb, nch, erg)
c     WARNING:
c     erg may be overwritten by this subroutine.
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer rtdb                ! [input]
      integer nch                 ! [input]
      double precision erg(nch)   ! [in/output]

      logical umanerg(nch)
      integer i, l_eg, k_eg
      

      if (.not. MA_alloc_get(MT_DBL, nch, 'manerg', l_eg, k_eg))
     &     call errquit('WFAT: manerg allocation failed', 1, MA_ERR)
      
      
c     ==== When using manual dipole ====
      if (rtdb_get(rtdb, 'wfat:use_manerg', mt_log, nch, umanerg)) then
         do i = 1, nch
            if (umanerg(i)) then
               if (rtdb_get(rtdb, 'wfat:manerg', mt_dbl, nch,
     &              dbl_mb(k_eg))) then
c                 dbl_mb(k_eg) contains alternatives to orbital energies.
                  exit
               else
                  call errquit('WFAT: Manual energy is used but '//
     &                 'cannot find it in RTDB.', 1, RTDB_ERR)
               endif
            endif
         enddo
      else
         umanerg(:) = .false.
      endif

      
c     ==== Replace proper orb. energy with the custom ones, if any ====
      do i = 1, nch
         if (umanerg(i)) erg(i) = dbl_mb(k_eg + i - 1)
      enddo
      

      if (.not. MA_free_heap(l_eg))
     &     call errquit('WFAT: manerg deallocation failed', 1, MA_ERR)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_man_chpos(rtdb, nch, umanpos, chpos)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      integer nch
      logical umanpos(nch)
      double precision chpos(3,nch)

      integer i

      
c     ==== When using manual dipole ====
      if (rtdb_get(rtdb, 'wfat:use_manpos', mt_log, nch, umanpos)) then
         do i = 1, nch
            if (umanpos(i)) then
               if (rtdb_get(rtdb, 'wfat:manpos', mt_dbl, 3*nch,
     &              chpos)) then
                  exit
               else
                  call errquit('WFAT: Manual chpos is used but '//
     &                 'cannot find it in RTDB.', 1, RTDB_ERR)
               endif
            endif
         enddo
      else
         umanpos(:) = .false.
      endif


      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_get_man_orbcent_OLD(rtdb, nch, orbcent, orbmmt,
     &     iondip)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      integer nch
      double precision orbcent(3, nch)
      double precision orbmmt(3, nch)
      double precision iondip(3, nch)

      integer i

      
      do i = 1, nch
         if (.not. rtdb_get(rtdb, 'wfat:man_orbcent', mt_dbl, 3,
     &        orbcent(:,i)))
     &        call errquit('WFAT: Cannot get man_orbcent from RTDB',
     &        1, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'wfat:man_orbmmt', mt_dbl, 3,
     &        orbmmt(:,i)))
     &        call errquit('WFAT: Cannot get man_orbmmt from RTDB',
     &        1, RTDB_ERR)
         iondip(:,i) = (/ 1.0d3, 1.0d3, 1.0d3 /)
      enddo

      end
c=======================================================
