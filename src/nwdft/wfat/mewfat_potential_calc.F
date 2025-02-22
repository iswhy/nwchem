c=======================================================
c=======================================================
      subroutine mewfat_hf_veepsi_1ptbas(bas, nb, n0, n1, sx, sy, vpt,
     &     bval, ocmo0, ocmo1, pcf, rcf, qcf, scf, scr, rs, si_tp,
     &     vhfpsi)

      implicit none
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      double precision vpt(3)                         ! [input], evaluation point.
      double precision bval(nb)                       ! [input], basis functions evaluated at vpt.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision scr(1:*)                       ! [input], scratch array.   1)
      logical rs                                      ! [input], controls whether using range-separated (rs=.true.) exchange or not.
      character*5 si_tp                               ! [input], self-interaction type.
      double precision vhfpsi                         ! [output], one-determinant, e-e contribution to ionization.

      logical srhfx, isdft
      integer pmpt0, pmptn, pc0, pcn, ptm0, ptmn, ptmn_, pcmo0, pcmon
      integer per0, pern, ptrs0, ptrsn, ptrsn_, pcrs0, pcrsn
      integer nmpt, nc, ntm, ncmo, ner, ntrs, ncrs, ndid
      integer i, io, jo, is, c1, c2, j1, j2, icmo0, icmon, icrs0, icrsn
      integer j, dl, ntot0, l_ds, k_ds, i_ds
      double precision crg, same_sp, comp_sp, chfx, alpha, beta, gamma
      double precision si, si0, vmo0 ,vrs0, BETA_THR, ds, fp, f0

      parameter (BETA_THR=1.0d-13)

      logical nwxc_has_xc
      double precision nwxc_wght_hfx
      external nwxc_wght_hfx, nwxc_has_xc


      ndid = ga_nodeid()
      

c     ==== Check no. of electrons ====
      call mewfat_check_nelec(sx, sy, n0, n1, 100)

      
c     ==== pmpt0 = pointer to neutral MO in sx channel ====
      nmpt = n0(sx)
      pmpt0 = 1
      pmptn = pmpt0 + nmpt - 1

c     ==== pc0 = pointer to Coulomb matrix in basis ====
      nc = nb**2
      pc0 = pmptn + 1
      pcn = pc0 + nc - 1

c     ==== ptm0 = pointer to a temporary array ====
      ntm = nb * max(maxval(n0),maxval(n1))
      ptm0 = pcn + 1
      ptmn = ptm0 + ntm - 1

c     ==== pcmo0 = pointer to Coulomb matrix in occ. MO's ====
      ncmo = sum(n1*n0)
      pcmo0 = ptmn + 1
      pcmon = pcmo0 + ncmo - 1

c     ==== per0 = pointer to erf Coulomb matrix in basis ====
      ner = nb**2
      per0 = pcmon + 1
      pern = per0 + ner - 1
      
c     ==== ptrs0 = pointer to a temporary array for RS ====
      ntrs = nb * max(maxval(n0),maxval(n1))
      ptrs0 = pern + 1
      ptrsn = ptrs0 + ntrs - 1      

c     ==== pcrs0 = pointer to RS Coulomb matrix in occ. MO's ====
      ncrs = sum(n1*n0)
      pcrs0 = ptrsn + 1
      pcrsn = pcrs0 + ncrs - 1
      
      
c     ==== Calculate occupied neutral MO's ====
c     ====   in the ionized spin channel   ====
      do i = 1, n0(sx)
         io = pmpt0 + i - 1
         jo = (sx-1)*n0(1) + i
         scr(io) = sum(bval*ocmo0(:,jo))
      enddo
      

c     ==== Calculate Coulomb matrix in basis ====
      crg = 1.0d0
      call int_1ehyd_all(bas, nb, vpt, crg, scr(pc0:pcn)) !! scr = 1e coulomb interaction matrix

      if (rs) then
         call nwxc_get_cam(alpha, beta, gamma, srhfx)
         if (srhfx) call errquit('WFAT: Use of range-separated '//
     &        'exchange in WFAT at the moment does not allow '//
     &        'short range HF.', 1, UERR)
         if (beta > BETA_THR) then
            call wfat_1ehyd_erf_all(bas, nb, gamma, vpt, crg,
     &        scr(per0:pern))
         else
            scr(per0:pern) = 0.0d0
         endif
      endif
      

c     ==== Calculate Coulomb matrix between ====
c     ====     cation and neutral MO's      ====
      do is = 1, 2
c        ==== U * MO1 ====
         c1 = (is-1)*n1(1) + 1
         c2 = c1 + n1(is) - 1
         ptmn_ = ptm0 + nb*n1(is) - 1
         call wfat_dsymm_wrap(nb, n1(is), 1.0d0, scr(pc0:pcn), 
     &        ocmo1(:, c1:c2), 0.0d0, scr(ptm0:ptmn_))
         if (rs) then                         ! CHECK: check correctness here
            ptrsn_ = ptrs0 + nb*n1(is) - 1
            call wfat_dsymm_wrap(nb, n1(is), 1.0d0, scr(per0:pern),
     &           ocmo1(:, c1:c2), 0.0d0, scr(ptrs0:ptrsn_))
         endif
         
c        ==== MO0^T * U ====
         c1 = (is-1)*n0(1) + 1
         c2 = c1 + n0(is) - 1
         j1 = (is-1)*n0(1)*n1(1) + pcmo0
         j2 = j1 + n0(is)*n1(is) - 1
         call wfat_dgemm_wrap('t', 'n', nb, n0(is), nb, n1(is), 1.0d0, 
     &        ocmo0(:, c1:c2), scr(ptm0:ptmn_), 0.0d0, scr(j1:j2))
         if (rs) then
            j1 = (is-1)*n0(1)*n1(1) + pcrs0
            j2 = j1 + n0(is)*n1(is) - 1
            call wfat_dgemm_wrap('t', 'n', nb, n0(is), nb, n1(is), 
     &           1.0d0, ocmo0(:, c1:c2), scr(ptrs0:ptrsn_), 0.0d0,
     &           scr(j1:j2))
         endif
      enddo
      
            
c     ==== Calculate same spin Hartree+exchange term ====
      chfx = nwxc_wght_hfx()
      icmo0 = (sx-1)*n0(1)*n1(1) + pcmo0
      icmon = icmo0 + n0(sx)*n1(sx) - 1
      if (rs) then
         icrs0 = (sx-1)*n0(1)*n1(1) + pcrs0          ! CHECK: check correctness here
         icrsn = icrs0 + n0(sx)*n1(sx) - 1
         call mewfat_hf_samespin_rs(n0(sx), n1(sx), 
     &        scr(icmo0:icmon), scr(icrs0:icrsn), scr(pmpt0:pmptn),
     &        rcf, qcf, chfx, alpha, beta, same_sp)
      else
         call mewfat_hf_samespin(n0(sx), n1(sx), scr(icmo0:icmon),     ! CHECK: add hfexch coefficient as input
     &        scr(pmpt0:pmptn), rcf, qcf, chfx, same_sp)
      endif

      
c     ==== Calculate different spin Hartree term ====
      icmo0 = (sy-1)*n0(1)*n1(1) + pcmo0
      icmon = icmo0 + n0(sy)*n1(sy) - 1
      call mewfat_hf_compspin(sx, sy, n0, n1(sy), scr(icmo0:icmon), 
     &     scr(pmpt0:pmptn), scf, pcf, comp_sp)
      

c     ==== DFT self-interaction ====
      isdft = nwxc_has_xc()
      if (isdft) then
         if (si_tp == 'none') then
            si = 0.0d0
         else if (si_tp == 'mo') then
            call mewfat_siterm_ocmo(nb, n0, sx, chfx, rs, alpha, beta,
     &           ocmo0, scr(pmpt0:pmptn), scr(pc0:pcn), scr(per0:pern),
     &           pcf, rcf, si)
         else if (si_tp == 'dyson') then
            call mewfat_siterm_dyson(nb, n0, sx, chfx, rs, alpha, beta,
     &           ocmo0, scr(pmpt0:pmptn), scr(pc0:pcn), scr(per0:pern),
     &           pcf, rcf, si)
         endif
      else
         si = 0.0d0
      endif

      
c     ==== Calculate HF ee contribution ====
      dl = 0
      if (sx==2) dl = 1
      ntot0 = sum(n0)
      vhfpsi = (-1)**(ntot0+dl*n0(1)) / sqrt(dble(ntot0)) *
     &     (same_sp + si + comp_sp)       ! check here if ntot0 is independent of choice of which electron is removed from neutral.

cfordebug      if (ndid==0) write(luout, '("    vee= ", 3es14.5)')
cfordebug     &     same_sp, same_sp/rcf, comp_sp
      
      
      end
c     NOTES:
c     1) Minimum number of elements should be n0(sx) + nb**2 + nb * max(maxval(n0),maxval(n1))
c        + sum(n1*n0) + nrs, where nrs = 0 when rs=.false. and nrs = nb**2 + 
c        nb*max(maxval(n0),maxval(n1)) + sum(n1*n0) when rs=.true..
c=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_hf_samespin(n0, n1, vmo, phi, rcf, qcf, chfx, h)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"      
      integer n0                         ! [input], no. of electrons in the ionized spin channel in the neutral.
      integer n1                         ! [input], no. of electrons in the ionized spin channel in the cation.
      double precision vmo(n0, n1)       ! [input], the ionized spin block of Coulomb matrix. The rows index neutral MO's while the columns index cation MO's.
      double precision phi(n0)           ! [input], neutral MO's in the ionized spin channel evaluated at the desired point.
      double precision rcf               ! [input], R det. coefficient.
      double precision qcf(n0, n0, n1)   ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision chfx              ! [input],
      double precision h                 ! [output], contribution of the ionized spin channel to ionization.

      integer kk, j, k
      double precision h0

      
      if (n1 /= n0-1)
     &     call errquit('mewfat_hf_samespin: Inconsistent values of '//
     &     'n0 and n1', 1, UERR)

      h = 0.0d0
      do kk = 1, n1
         do j = 1, n0
            do k = 1, j-1
               h0 = qcf(k,j,kk) * (vmo(k,kk)*phi(j) - 
     &              chfx*vmo(j,kk)*phi(k)) * (-1)**(j+k+kk)
               h = h + h0
            enddo
         enddo
      enddo

cdebug      h = h + vmo(n0,n0)*phi(n0) - chfx*vmo(n0,n0)*phi(n0)    this line is for trying to reproduce OEWFAT with DFT, it turns out vmo(n0,n0) is nonexistent.
      
      h = h * rcf

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_hf_samespin_rs(n0, n1, vmo, vrs, phi, rcf, qcf, 
     &     chfx, alpha, beta, h)

      implicit none
#include "errquit.fh"
      integer n0                         ! [input], no. of electrons in the ionized spin channel in the neutral.
      integer n1                         ! [input], no. of electrons in the ionized spin channel in the cation.
      double precision vmo(n0, n1)       ! [input], the ionized spin block of Coulomb matrix. The rows index neutral MO's while the columns index cation MO's.
      double precision vrs(n0, n1)       ! [input], the ionized spin block of Coulomb matrix. The rows index neutral MO's while the columns index cation MO's.
      double precision phi(n0)           ! [input], neutral MO's in the ionized spin channel evaluated at the desired point.
      double precision rcf               ! [input], R det. coefficient.
      double precision qcf(n0, n0, n1)   ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision chfx              ! [input],
      double precision alpha             ! [input],
      double precision beta              ! [input],
      double precision h                 ! [output], contribution of the ionized spin channel to ionization.

      integer kk, j, k
      double precision h0

      
      if (n1 /= n0-1)
     &     call errquit('mewfat_hf_samespin: Inconsistent values of '//
     &     'n0 and n1', 1, UERR)
      
      h = 0.0d0
      do kk = 1, n1
         do j = 1, n0
            do k = 1, j-1
               h0 = qcf(k,j,kk) * (vmo(k,kk)*phi(j) - 
     &              chfx*(alpha*vmo(j,kk) + beta*vrs(j,kk))*phi(k))
     &              * (-1)**(j+k+kk)
               h = h + h0
            enddo
         enddo
      enddo
      h = h * rcf

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine mewfat_hf_compspin(sx, sy, n0, n1, vmo, phi, scf,
     &     pcf, h)

      implicit none
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer sx                        ! [input], sx=1 or 2, the ionized electron spin.
      integer sy                        ! [input], if sx=1 then sy=2, and vice versa, the complementary spin.
      integer n0(2)                     ! [input], no. of electrons in both spin channels in the neutral.
      integer n1                        ! [input], no. of electrons in the complementary spin channel, sy, in the cation.
      double precision vmo(n0(sy), n1)  ! [input], the complementary spin (sy) block of Coulomb matrix. The rows index neutral MO's while the columnds index cation MO's.
      double precision phi(n0(sx))      ! [input], neutral MO's in the ionized spin channel evaluated at the desired point.
      double precision scf(n0(sy), n1)  ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision pcf(n0(sx))      ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision h                ! [output], contribution of the complementary spin channel to ionization.

      integer kk, k, j, ndid
      double precision g, g0, f, f0, h0


      ndid = ga_nodeid()
      call mewfat_check_spin(sx, sy, 101)
      

      if (n1 /= n0(sy))
     &     call errquit('mewfat_hf_compspin: Inconsistent values of '//
     &     'n0(sy) and n1', 1, UERR)

      g = 0.0d0
      do kk = 1, n1
         do k = 1, n0(sy)
            g0 = (-1)**(kk+k) * scf(k,kk) * vmo(k,kk)
            g = g + g0
         enddo
      enddo

      f = 0.0d0
      do j = 1, n0(sx)
         f0 = (-1)**j * pcf(j) * phi(j)
         f = f + f0
      enddo

      h = g * f

cfordebug      if (ndid==0) write(luout, '("    vcomp= ", 4es14.5)')
cfordebug     &     scf(1,1), vmo(1,1), pcf(1), phi(1)
      
      end
c=======================================================


cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine mewfat_vxcpsi_1ptbas_OLD(nb, sx, n0, bval, d1bval, d2bval,
cOLD     &           ocmo0, rdm0, pcf, rcf, smh, xclvl, scr, vxcpsi)
cOLD
cOLD      implicit none
cOLD#include "stdio.fh"
cOLD#include "wfat_indexing.fh"
cOLD#include "errquit.fh"
cOLD      integer nb                               ! [input], no. of basis.
cOLD      integer sx                               ! [input], ionized spin channel.
cOLD      integer n0(2)                            ! [input], no. of electrons in the neutral.
cOLD      double precision bval(nb)                ! [input], basis functions evaluated at the desired point.
cOLD      double precision d1bval(3,nb)            ! [input], 1st derivative of basis functions evaluated at the desired point.   1)
cOLD      double precision d2bval(6,nb)            ! [input], 2nd derivative of basis functions evaluated at the desired point.   1)
cOLD      double precision ocmo0(nb, n0(1)+n0(2))  ! [input], occupied MO's of the neutral.
cOLD      double precision rdm0(nb,nb,2)           ! [input], RDM of the neutral (that corresponds to ocmo0).
cOLD      double precision pcf(n0(sx))             ! [input], 
cOLD      double precision rcf                     ! [input], 
cOLD      double precision smh                     ! [input], 
cOLD      character*4 xclvl                        ! [input],      1)
cOLD      double precision scr(1:*)                ! [output/scratch], min. no. of elements should be n0(sx).
cOLD      double precision vxcpsi                  ! [output],
cOLD
cOLD      integer nspin
cOLD      parameter (nspin=2)
cOLD      integer nmpt, pmpt0, pmptn, i, io, jo, maxd
cOLD      double precision vxc(nspin), f
cOLD      double precision dens(nspin), d1dens(3,nspin), d2dens(6,nspin)
cOLD      double precision gamma(3), d1gamma(3,3), tau(nspin), lapl(nspin)
cOLD      double precision fn
cOLD      double precision dfndr(nspin), dfndr2(NCOL_AMAT2)
cOLD      double precision dfndg(3), dfndg2(NCOL_CMAT2)
cOLD      double precision dfndt(nspin), dfndt2(NCOL_MMAT2)
cOLD
cOLD
cOLD      if (xclvl /= 'lda' .and. xclvl /= 'gga' .and. xclvl /= 'mgga'
cOLD     &    .and. xclvl /= 'none')
cOLD     &     call errquit('MEWFAT: Argument xclvl of mewfat_vxcpsi_'//
cOLD     &                  '1ptbas must be either lda, gga, mgga, or '//
cOLD     &                  'none.', 1, UERR)
cOLD      if (xclvl == 'none' .or. sum(n0)== 1) then
cOLD         vxcpsi = 0.0d0
cOLD         return
cOLD      endif
cOLD
cOLD      
cOLD      if (xclvl == 'lda') maxd = 0
cOLD      if (xclvl == 'gga' .or. xclvl == 'mgga') maxd = 2
cOLD      call wfat_dens_1ptbas(nb, nspin, rdm0, maxd,
cOLD     &                      bval, d1bval, d2bval,
cOLD     &                      dens, d1dens, d2dens)
cOLD
cOLD      gamma = 0.0d0
cOLD      tau = 0.0d0
cOLD      fn = 0.0d0
cOLD      dfndr = 0.0d0
cOLD      dfndr2 = 0.0d0
cOLD      dfndg = 0.0d0
cOLD      dfndg2 = 0.0d0
cOLD      dfndt = 0.0d0
cOLD      dfndt2 = 0.0d0
cOLD      if (xclvl == 'lda') then
cOLD         call nwxc_eval_df(nspin, 1, dens, gamma, tau, fn, dfndr,
cOLD     &                     dfndg, dfndt)
cOLD
cOLDc        ==== Calculate LDA potential ====
cOLD         call wfat_vxclda_1pt(nspin, dfndr, vxc)
cOLD         
cOLD      else if (xclvl == 'gga' .or. xclvl == 'mgga') then
cOLD         call wfat_gamma_1pt(nspin, d1dens, d2dens, gamma, d1gamma)
cOLDccheck         if (xclvl == 'mgga') 
cOLDccheck  &           call wfat_tau_1pt(nb, nspin, rdm, d1bval, tau)
cOLD         call wfat_tau_1pt(nb, nspin, rdm0, d1bval, tau)
cOLD         
cOLD         call nwxc_eval_df2(nspin, 1, dens, gamma, tau, fn, dfndr, 
cOLD     &                      dfndr2, dfndg, dfndg2, dfndt, dfndt2)
cOLD
cOLD         call wfat_lapl_1pt(nb, nspin, d2dens, lapl)
cOLD
cOLD         
cOLDc        ==== Calculate GGA or mGGA potential ====
cOLD         if (xclvl == 'gga') then
cOLD            call wfat_vxcgga_1pt(nspin, d1dens, lapl, gamma, d1gamma,
cOLD     &                           dfndr, dfndg, dfndg2, vxc)
cOLD         else if(xclvl == 'mgga') then
cOLD            call wfat_vxcmgga_1pt()
cOLD         endif
cOLD      endif
cOLD
cOLD      
cOLDc     ==== Calculate occupied neutral MO's ====
cOLDc     ====   in the ionized spin channel   ====
cOLD      nmpt = n0(sx)
cOLD      pmpt0 = 1
cOLD      pmptn = pmpt0 + nmpt - 1
cOLD      do i = 1, n0(sx)
cOLD         io = pmpt0 + i - 1
cOLD         jo = (sx-1)*n0(1) + i
cOLD         scr(io) = sum(bval*ocmo0(:,jo))
cOLD      enddo      
cOLD      call mewfat_eval_dyson(sx, n0, pcf, rcf, scr(pmpt0:pmptn), f)
cOLD
cOLD
cOLDc     ==== Calculate v1psi ====
cOLD      vxcpsi = vxc(sx) * f
cOLD      
cOLD      end
cOLDc     Notes:
cOLDc       1) d1bval and d2bval are used only when xclvl is equal to 'gga' or 'mgga'.
cOLDc=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_vxcpsi_1ptbas(nb, nspin, sx, n0, bval, 
     &           d1bval, d2bval, ocmo0, rdm0, pcf, rcf, smh, xclvl, scr,
     &           vxcpsi)

      implicit none
#include "stdio.fh"
#include "wfat_indexing.fh"
#include "errquit.fh"
      integer nb                ! [input], no. of basis.
      integer nspin
      integer sx                               ! [input], ionized spin channel.
      integer n0(2)                            ! [input], no. of electrons in the neutral.
      double precision bval(nb)                ! [input], basis functions evaluated at the desired point.
      double precision d1bval(3,nb)            ! [input], 1st derivative of basis functions evaluated at the desired point.   1)
      double precision d2bval(6,nb)            ! [input], 2nd derivative of basis functions evaluated at the desired point.   1)
      double precision ocmo0(nb, n0(1)+n0(2))  ! [input], occupied MO's of the neutral.
      double precision rdm0(nb,nb,2)           ! [input], RDM of the neutral (that corresponds to ocmo0).
      double precision pcf(n0(sx))             ! [input], 
      double precision rcf                     ! [input], 
      double precision smh                     ! [input], 
      character*4 xclvl                        ! [input],      1)
      double precision scr(1:*)                ! [output/scratch], min. no. of elements should be n0(sx).
      double precision vxcpsi                  ! [output],

      integer nmpt, pmpt0, pmptn, i, io, jo, maxd
      double precision vxc(nspin), f
      double precision dens(nspin), d1dens(3,nspin), d2dens(6,nspin)
      double precision gamma(3), d1gamma(3,3), tau(nspin), lapl(nspin)
      double precision fn
      double precision dfndr(nspin), dfndr2(NCOL_AMAT2)
      double precision dfndg(3), dfndg2(NCOL_CMAT2)
      double precision dfndt(nspin), dfndt2(NCOL_MMAT2)


      if (xclvl /= 'lda' .and. xclvl /= 'gga' .and. xclvl /= 'mgga'
     &    .and. xclvl /= 'none')
     &     call errquit('MEWFAT: Argument xclvl of mewfat_vxcpsi_'//
     &                  '1ptbas must be either lda, gga, mgga, or '//
     &                  'none.', 1, UERR)
      if (xclvl == 'none' .or. sum(n0)== 1) then
         vxcpsi = 0.0d0
         return
      endif

      
      if (xclvl == 'lda') maxd = 0
      if (xclvl == 'gga' .or. xclvl == 'mgga') maxd = 2
      call wfat_dens_1ptbas(nb, nspin, rdm0, maxd,
     &                      bval, d1bval, d2bval,
     &                      dens, d1dens, d2dens)

      gamma = 0.0d0
      tau = 0.0d0
      fn = 0.0d0
      dfndr = 0.0d0
      dfndr2 = 0.0d0
      dfndg = 0.0d0
      dfndg2 = 0.0d0
      dfndt = 0.0d0
      dfndt2 = 0.0d0
      if (xclvl == 'lda') then
         call nwxc_eval_df(nspin, 1, dens, gamma, tau, fn, dfndr,
     &                     dfndg, dfndt)

c        ==== Calculate LDA potential ====
         call wfat_vxclda_1pt(nspin, dfndr, vxc)
         
      else if (xclvl == 'gga' .or. xclvl == 'mgga') then
         call wfat_gamma_1pt(nspin, d1dens, d2dens, gamma, d1gamma)
ccheck         if (xclvl == 'mgga') 
ccheck  &           call wfat_tau_1pt(nb, nspin, rdm, d1bval, tau)
         call wfat_tau_1pt(nb, nspin, rdm0, d1bval, tau)
         
         call nwxc_eval_df2(nspin, 1, dens, gamma, tau, fn, dfndr, 
     &                      dfndr2, dfndg, dfndg2, dfndt, dfndt2)

         call wfat_lapl_1pt(nb, nspin, d2dens, lapl)

         
c        ==== Calculate GGA or mGGA potential ====
         if (xclvl == 'gga') then
            call wfat_vxcgga_1pt(nspin, d1dens, lapl, gamma, d1gamma,
     &                           dfndr, dfndg, dfndg2, vxc)
         else if(xclvl == 'mgga') then
            call wfat_vxcmgga_1pt()
         endif
      endif

      
c     ==== Calculate occupied neutral MO's ====
c     ====   in the ionized spin channel   ====
      nmpt = n0(sx)
      pmpt0 = 1
      pmptn = pmpt0 + nmpt - 1
      do i = 1, n0(sx)
         io = pmpt0 + i - 1
         jo = (sx-1)*n0(1) + i
         scr(io) = sum(bval*ocmo0(:,jo))
      enddo      
      call mewfat_eval_dyson(sx, n0, pcf, rcf, scr(pmpt0:pmptn), f)


c     ==== Calculate v1psi ====
      if (nspin == 1) then
         vxcpsi = vxc(1) * f
      else
         vxcpsi = vxc(sx) * f
      endif
      
      end
c     Notes:
c       1) d1bval and d2bval are used only when xclvl is equal to 'gga' or 'mgga'.
c=======================================================
      

c=======================================================
c=======================================================
      subroutine mewfat_v1psi_1ptbas(geom, nb, sx, n0, vpt, bval,
     &     ocmo0, pcf, rcf, zi, ccent, smh, scr, v1psi)

      implicit none
#include "errquit.fh"
#include "geom.fh"
      integer geom
      integer nb
      integer sx
      integer n0(2)
      double precision vpt(3)
      double precision bval(nb)
      double precision ocmo0(nb, n0(1)+n0(2))
      double precision pcf(n0(sx))
      double precision rcf
      double precision zi
      double precision ccent(3)
      double precision smh
      double precision scr(1:*)                  ! [output/scratch], min. no. of elements should be n0(sx).
      double precision v1psi

      character*16 tag
      integer nat, ia, nmpt, pmpt0, pmptn, i, io, jo
      double precision vn, vh, atp(3), zn, dv(3), dist, f


      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 3, GEOM_ERR)


c     ==== Nuclear attraction potential ====
      vn = 0.0d0
      do ia = 1, nat
         if (.not. geom_cent_get(geom, ia, tag, atp, zn))
     &        call errquit('WFAT: Getting zn and atp failed',
     &             4, GEOM_ERR)
         dv = vpt - atp
         dist = sqrt(sum(dv*dv) + smh)
         if (abs(dist) < 1.0d-12) then
            smh = 1.0d-8
            dist = sqrt(sum(dv*dv) + smh)
         endif
         vn = vn + (-zn/dist)
      enddo


c     ==== H-like potential ====
      dv = vpt - ccent
      dist = sqrt(sum(dv*dv) + smh)
      if (abs(dist) < 1.0d-12) then
         smh = 1.0d-8
         dist = sqrt(sum(dv*dv) + smh)
      endif
      vh = zi/dist


c     ==== Calculate occupied neutral MO's ====
c     ====   in the ionized spin channel   ====
      nmpt = n0(sx)
      pmpt0 = 1
      pmptn = pmpt0 + nmpt - 1
      do i = 1, n0(sx)
         io = pmpt0 + i - 1
         jo = (sx-1)*n0(1) + i
         scr(io) = sum(bval*ocmo0(:,jo))
      enddo      
      call mewfat_eval_dyson(sx, n0, pcf, rcf, scr(pmpt0:pmptn), f)


c     ==== Calculate v1psi ====
      v1psi = (vn + vh) * f

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_siterm_ocmo(nb, n0, sx, chfx, rs, alpha, beta,
     &           ocmo0, phi, u_cl, u_er, pcf, rcf, si)

      implicit none
      integer nb                                 ! [input]
      integer n0(2)                              ! [input]
      integer sx                                 ! [input]
      double precision chfx                      ! [input]
      logical rs                                 ! [input]
      double precision alpha                     ! [input], only used when rs=.true..
      double precision beta                      ! [input], only used when rs=.true..
      double precision ocmo0(nb, n0(1)+n0(2))    ! [input]
      double precision phi(n0(sx))               ! [input]
      double precision u_cl(nb,nb)               ! [input]
      double precision u_er(nb,nb)               ! [input], only used when rs=.true..
      double precision pcf(n0(sx))               ! [input]
      double precision rcf                       ! [input]
      double precision si                        ! [output] 

      integer i, jo
      double precision si0, vmo0, vrs0
      

      do i = 1, n0(sx)
c        ==== Calculate the Self-Interaction potential ====
         jo = (sx-1)*n0(1) + i
         call wfat_1e_expval_1eopr(nb, ocmo0(:, jo), u_cl, vmo0)
         if (rs) then
            call wfat_1e_expval_1eopr(nb, ocmo0(:, jo), u_er, vrs0)
         endif

c        ==== Calculate the self-interaction term ====
         if (rs) then
            si0 = pcf(i) * (-1)**i *
     &           (vmo0 - chfx*(alpha*vmo0 + beta*vrs0)) * phi(i)
         else               
            si0 = pcf(i) * (-1)**i * (1.0d0-chfx) * vmo0 * phi(i)
         endif
         si = si + si0
      enddo
      si = rcf * si


      end
c=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_siterm_dyson(nb, n0, sx, chfx, rs, alpha, beta,
     &           ocmo0, phi, u_cl, u_er, pcf, rcf, si)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nb                                 ! [input]
      integer n0(2)                              ! [input]
      integer sx                                 ! [input]
      double precision chfx                      ! [input]
      logical rs                                 ! [input]
      double precision alpha                     ! [input], only used when rs=.true..
      double precision beta                      ! [input], only used when rs=.true..
      double precision ocmo0(nb, n0(1)+n0(2))    ! [input]
      double precision phi(n0(sx))               ! [input]
      double precision u_cl(nb,nb)               ! [input]
      double precision u_er(nb,nb)               ! [input], only used when rs=.true..
      double precision pcf(n0(sx))               ! [input]
      double precision rcf                       ! [input]
      double precision si                        ! [output] 

      integer j, i, jo, i_ds, l_ds, k_ds
      double precision f0, fp, ds, vmo0, vrs0
      

c     ==== Evaluate imperfect Dyson orbital ====
      fp = 0.0d0
      do j = 1, n0(sx)
         f0 = (-1)**j * pcf(j) * phi(j)
         fp = fp + f0
      enddo


c     ==== Compute the coefficients of imperfect Dyson orbitals ====
      if (.not. ma_alloc_get(mt_dbl, nb, 'ds2', l_ds, k_ds))
     &     call errquit('MEWFAT: Allocating ds2 failed', 1, MA_ERR)
      do i = 1, nb
         i_ds = k_ds + i - 1
         dbl_mb(i_ds) = 0.0d0
         do j = 1, n0(sx)
            jo = (sx-1)*n0(1) + j
            ds = (-1)**j * pcf(j) * ocmo0(i,jo)
            dbl_mb(i_ds) = dbl_mb(i_ds) + ds
         enddo
      enddo


c     ==== Calculate the Self-Interaction potential ====
      call wfat_1e_expval_1eopr(nb, dbl_mb(k_ds), u_cl, vmo0)
      if (rs) then
         call wfat_1e_expval_1eopr(nb, dbl_mb(k_ds), u_er, vrs0)
      endif


c     ==== Calculate the self-interaction term ====
      if (rs) then
         si = rcf * (vmo0 - chfx*(alpha*vmo0 + beta*vrs0)) * fp
      else               
         si = rcf * (1.0d0-chfx) * vmo0 * fp
      endif
      

c     ==== Deallocation ====
      if (.not. ma_free_heap(l_ds))
     &     call errquit('MEWFAT: Failed to deallocate ds2',
     &     1, MA_ERR)

      end
c=======================================================
