c=======================================================
c=======================================================
      subroutine mewfat_hf_vcpsi_1pt_c(icent, geom, bas, nb, sx, sy, 
     &     n0, n1, vpt, cgrid, ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, 
     &     ccent, smh, scr, vcpsi)

      implicit none
      integer icent
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      double precision vpt(3)                         ! [input], evaluation point.
      logical cgrid          
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input],
      double precision ccent(3)                       ! [input],
      double precision smh                            ! [input],
      double precision scr(1:*)                       ! [input], scratch array. 3)
      double precision vcpsi                          ! [output],


      if (icent /= -1 .or. (icent == -1 .and. cgrid)) then
         call mewfat_hf_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &        vpt, ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, ccent, smh,
     &        scr, vcpsi)
      else
         vcpsi = 0.0d0
      endif
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine mewfat_hf_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &     vpt, ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, ccent, smh, 
     &     scr, vcpsi)

      implicit none
#include "global.fh"
#include "stdio.fh"      
#include "errquit.fh"
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      double precision vpt(3)                         ! [input], evaluation point.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input],
      double precision ccent(3)                       ! [input],
      double precision smh                            ! [input],
      double precision scr(1:*)                       ! [input], scratch array. 3)
      double precision vcpsi                          ! [output],

c     Notes:
c     3) The minimum number of elements should be 
c           21*nb + n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)

      integer nb0, pb0, pb0n, nb1, pb1, pb1n, nb2, pb2, pb2n
      integer nb3, pb3, pb3n, nsc, psc, pscn, pscn_, nrf, prf, prfn
      integer ndid
      double precision vhfpsi, v1psi

      ndid = ga_nodeid()


c     ==== pb0 = pointer to evaluated bases ====
      nb0 = nb
      pb0 = 1
      pb0n = pb0 + nb0 - 1

c     ==== pb1 = pointer to evaluated 1st derivative bases ====
      nb1 = 3 * nb
      pb1 = pb0n + 1
      pb1n = pb1 + nb1 - 1
      
c     ==== pb2 = pointer to evaluated 2nd derivative bases ====
      nb2 = 6 * nb
      pb2 = pb1n + 1
      pb2n = pb2 + nb2 - 1

c     ==== pb3 = pointer to evaluated 3rd derivative bases ====
      nb3 = 10 * nb
      pb3 = pb2n + 1
      pb3n = pb3 + nb3 - 1

c     ==== psc = pointer to scratch array ====
      nsc = n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)
      psc = pb3n + 1
      pscn = psc + nsc - 1

c     ==== prf = pointer to the radial part of basis functions ====
      nrf = nb
      prf = pscn + 1
      prfn = prf + nrf - 1
      

c     ==== Calculate basis and its derivatives ====
      call wfat_eval_basis_radf(geom, bas, 0, vpt, scr(prf : prfn),
     &                         scr(pb0 : pb0n), scr(pb1 : pb1n),
     &                         scr(pb2 : pb2n), scr(pb3 : pb3n))

      
c     ==== Calculate vhf_psi ====
      call mewfat_hf_veepsi_1ptbas(bas, nb, n0, n1, sx, sy, vpt,
     &     scr(pb0:pb0n), ocmo0, ocmo1, pcf, rcf, qcf, scf, 
     &     scr(psc:pscn), .false., 'none', vhfpsi)

      
c     ==== Calculate v1_psi = vnuc_psi + vs_psi ====
      pscn_ = psc + n0(sx) - 1
      call mewfat_v1psi_1ptbas(geom, nb, sx, n0, vpt, scr(pb0:pb0n),
     &     ocmo0, pcf, rcf, zi, ccent, smh, scr(psc:pscn_), v1psi)


c     ==== Calculate vc_psi ====
      vcpsi = vhfpsi + v1psi

cfordebug      if (ndid==0) write(luout, '("    v= ", 4es14.5)')
cfordebug     &     vhfpsi, vhfpsi/rcf, v1psi, v1psi/rcf
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_dft_vcpsi_1pt_c(icent, geom, bas, nb, sx, sy, 
     &     n0, n1, vpt, cgrid, ocmo0, ocmo1, rdm0, pcf, rcf, qcf, scf, 
     &     zi, ccent, smh, rs, si_type, scr, vcpsi)

      implicit none
      integer icent
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      double precision vpt(3)                         ! [input], evaluation point.
      logical cgrid
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision rdm0(nb,nb,2)                  ! [input], RDM of the neutral (that corresponds to ocmo0).
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input],
      double precision ccent(3)                       ! [input],
      double precision smh                            ! [input],
      logical rs                                      ! [input],
      character*5 si_type                             ! [input], self-interaction type.
      double precision scr(1:*)                       ! [input], scratch array.  3)
      double precision vcpsi                          ! [output],

      
      if (icent /= -1 .or. (icent == -1 .and. cgrid)) then
         call mewfat_dft_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &     vpt, ocmo0, ocmo1, rdm0, pcf, rcf, qcf, scf, zi, ccent, smh,
     &     rs, si_type, scr, vcpsi)
      else
         vcpsi = 0.0d0
      endif

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine mewfat_dft_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &     vpt, ocmo0, ocmo1, rdm0, pcf, rcf, qcf, scf, zi, ccent, smh,
     &     rs, si_type, scr, vcpsi)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      double precision vpt(3)                         ! [input], evaluation point.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision rdm0(nb,nb,2)                  ! [input], RDM of the neutral (that corresponds to ocmo0).
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input],
      double precision ccent(3)                       ! [input],
      double precision smh                            ! [input],
      logical rs                                      ! [input],
      character*5 si_type                             ! [input], self-interaction type.
      double precision scr(1:*)                       ! [input], scratch array.  3)
      double precision vcpsi                          ! [output],

c     Notes:
c     3) The minimum number of elements should be 
c           21*nb + n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0) + n_rs
c        where
c           n_rs = nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)
c        when rs = .true. and n_rs = 0 otherwise.
      
      character*4 xclvl
      integer nb0, pb0, pb0n, nb1, pb1, pb1n, nb2, pb2, pb2n
      integer nb3, pb3, pb3n, nsc, psc, pscn, pscn_, nrf, prf, prfn
      integer dbas, nspin
      double precision vhfpsi, vxcpsi, v1psi

      logical xc_chktau
      external xc_chktau
      

c     ==== pb0 = pointer to evaluated bases ====
      nb0 = nb
      pb0 = 1
      pb0n = pb0 + nb0 - 1

c     ==== pb1 = pointer to evaluated 1st derivative bases ====
      nb1 = 3 * nb
      pb1 = pb0n + 1
      pb1n = pb1 + nb1 - 1
      
c     ==== pb2 = pointer to evaluated 2nd derivative bases ====
      nb2 = 6 * nb
      pb2 = pb1n + 1
      pb2n = pb2 + nb2 - 1

c     ==== pb3 = pointer to evaluated 3rd derivative bases ====
      nb3 = 10 * nb
      pb3 = pb2n + 1
      pb3n = pb3 + nb3 - 1

c     ==== psc = pointer to scratch array ====
      nsc = n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)    change to include rs
      psc = pb3n + 1
      pscn = psc + nsc - 1

c     ==== prf = pointer to the radial part of basis functions ====
      nrf = nb
      prf = pscn + 1
      prfn = prf + nrf - 1

      
c     ==== Determine DFT XC level ====
      if (xc_chktau()) then
         call errquit('MEWFAT: At the moment, mGGA functional is not '//
     &                'yet supported in WFAT.', 11, UERR)
      endif
      call wfat_get_xclvl(dbas, xclvl)
      

c     ==== Calculate basis and its derivatives ====
      call wfat_eval_basis_radf(geom, bas, dbas, vpt, scr(prf : prfn),
     &                         scr(pb0 : pb0n), scr(pb1 : pb1n),
     &                         scr(pb2 : pb2n), scr(pb3 : pb3n))

      
c     ==== Calculate v_psi ====
      call mewfat_hf_veepsi_1ptbas(bas, nb, n0, n1, sx, sy, vpt,
     &     scr(pb0:pb0n), ocmo0, ocmo1, pcf, rcf, qcf, scf, 
     &     scr(psc:pscn), rs, si_type, vhfpsi)


c     ==== Calculate vxc_psi ====
      pscn_ = psc + n0(sx) - 1
cOLD      call mewfat_vxcpsi_1ptbas_OLD(nb, sx, n0, scr(pb0), scr(pb1), 
cOLD     &     scr(pb2), ocmo0, rdm0, pcf, rcf, smh, xclvl, scr(psc:pscn_),
cOLD     &     vxcpsi)
      if (n0(1) == n0(2)) then
         nspin = 1
      else
         nspin = 2
      endif
      call mewfat_vxcpsi_1ptbas(nb, nspin, sx, n0, scr(pb0), 
     &     scr(pb1), scr(pb2), ocmo0, rdm0, pcf, rcf, smh, xclvl, 
     &     scr(psc:pscn_), vxcpsi)
      
      
c     ==== Calculate v1_psi = vnuc_psi + vs_psi ====
      pscn_ = psc + n0(sx) - 1
      call mewfat_v1psi_1ptbas(geom, nb, sx, n0, vpt, scr(pb0:pb0n),
     &     ocmo0, pcf, rcf, zi, ccent, smh, scr(psc:pscn_), v1psi)


c     ==== Calculate vc_psi ====
      vcpsi = vhfpsi + vxcpsi + v1psi

      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_collect_vcpsi(geom, bas, nb, sx, sy, n0, n1,
     &     nch, nptp, vpt, ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, ccent,
     &     smh, vc_psi)
      
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      integer nch
      integer nptp
      double precision vpt(3,nptp)                    ! [input], evaluation point.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input],
      double precision ccent(3,nch)                   ! [input],
      double precision smh                            ! [input],
      double precision vc_psi(nptp,nch)               ! [output],

      integer nscr, i
      integer l_sc, k_sc, l_vc, k_vc

      
      if (nch /= 1) call errquit('MEWFAT: At the moment, the argument'//
     &     ' nch of mewfat_collect_vcpsi can only be equal to 1.',
     &     32, UERR)
      
c     ==== Prepare scratch ====
      nscr = 21*nb + n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) +
     &      sum(n1*n0)
      if (.not. MA_alloc_get(MT_DBL, nscr, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 59, MA_ERR)
      
      if (.not. MA_alloc_get(MT_DBL, nch, 'vc_psi', l_vc, k_vc))
     &     call errquit('WFAT: vc_psi allocation failed', 25, MA_ERR)
      do i = 1, nptp
         call mewfat_hf_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &        vpt(:,i), ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, ccent, 
     &        smh, dbl_mb(k_sc), dbl_mb(k_vc))
         vc_psi(i,:) = dbl_mb(k_vc:k_vc+nch-1)
      enddo

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: scratch deallocation failed', 49, MA_ERR)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc_psi deallocation failed', 49, MA_ERR)

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine mewfat_collect_dft_vcpsi(geom, bas, nb, sx, sy, n0, n1,
     &     nch, nptp, vpt, ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, ccent,
     &     smh, rs, si_type, vc_psi)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      integer nch
      integer nptp
      double precision vpt(3,nptp)                    ! [input], evaluation point.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input],
      double precision ccent(3,nch)                   ! [input],
      double precision smh                            ! [input],
      logical rs                                      ! [input],
      character*5 si_type                             ! [input], self-interaction type.
      double precision vc_psi(nptp,nch)               ! [output],

      integer nsc_rs, nscr, i
      integer l_sc, k_sc, l_rd0, k_rd0, l_vc, k_vc

      
      if (nch /= 1) call errquit('MEWFAT: At the moment, the argument'//
     &     ' nch of mewfat_collect_dft_vcpsi can only be equal to 1.',
     &     32, UERR)
      
c     ==== Prepare scratch ====
      if (rs) then
         nsc_rs = nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)
      else
         nsc_rs = 0
      endif
      nscr = 21*nb + n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) +
     &       sum(n1*n0) + nsc_rs
      if (.not. MA_alloc_get(MT_DBL, nscr, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 59, MA_ERR)

c     ==== Obtain RDM of the neutral ====
      if (.not. ma_alloc_get(mt_dbl, nb*nb*2, 'rdm0', l_rd0, k_rd0))
     &     call errquit('MEWFAT: Cannot allocate rdm0', 61, MA_ERR)
      call mewfat_calc_rdm0(nb, n0, ocmo0, dbl_mb(k_rd0))
      
      if (.not. MA_alloc_get(MT_DBL, nch, 'vc_psi', l_vc, k_vc))
     &     call errquit('WFAT: vc_psi allocation failed', 25, MA_ERR)
      do i = 1, nptp
         call mewfat_dft_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &        vpt(:,i), ocmo0, ocmo1, dbl_mb(k_rd0), pcf, rcf, qcf,
     &        scf, zi, ccent, smh, rs, si_type, dbl_mb(k_sc),
     &        dbl_mb(k_vc))
         vc_psi(i,:) = dbl_mb(k_vc:k_vc+nch-1)
      enddo

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: scratch deallocation failed', 49, MA_ERR)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc_psi deallocation failed', 49, MA_ERR)
      if (.not. MA_free_heap(l_rd0))
     &     call errquit('WFAT: rdm0 deallocation failed', 49, MA_ERR)

      end
c=======================================================
      
