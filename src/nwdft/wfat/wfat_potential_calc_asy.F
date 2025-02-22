c=======================================================
c=======================================================
      subroutine wfat_dft_veeasy_1ptbas(bhan, nbas, nspin, vpt, bval,
     &                          d1bval, d2bval, bvalas, norb, vcf, 
     &                          vcfas, spinid, rdm, rdmas, xclvl, scr, 
     &                          veepsi)
      
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "wfat_qm_param.fh"      
      integer bhan                               ! [input] basis handle
      integer nbas                               ! [input] no. of basis corresponding to bhan
      integer nspin                              ! [input] number of spin (1=RHF/unpol., 2=UHF/pol.)
      double precision vpt(3)                    ! [input] evaluation point
      double precision bval(nbas)                ! [input] basis functions evaluated at vpt
      double precision d1bval(3, nbas)           ! [input] 1st derivative of basis functions evaluated at vpt   1)
      double precision d2bval(6, nbas)           ! [input] 2nd derivative of basis functions evaluated at vpt   1) 
      double precision bvalas(nbas)              ! [input] basis functions evaluated at vpt
      integer norb                               ! [input] no. of input orbitals
      double precision vcf(nbas, norb)           ! [input] coefficients of input orbitals
      double precision vcfas(nbas, norb)         ! [input] coefficients of input orbitals
      integer spinid(norb)                       ! [input] spin of input orbitals
      double precision rdm(nbas, nbas, nspin)    ! [input] reduced density matrix
      double precision rdmas(nbas, nbas, nspin)  ! [input] reduced density matrix
      character*4 xclvl                          ! [input] level of XC functional      1)
      double precision scr(nbas, *)              ! [output] scratch array, minimum no. of columns is 3*nbas+1+nspin.
      double precision veepsi(1:norb)            ! [output] output orbitals evaluated at vpt

      logical rs, srhfx
      integer l_hf, l_xc
      integer k_hf, k_xc
      double precision chfx, alpha, beta, gamma

      logical nwxc_has_cam
      double precision xc_hfexch, nwxc_wght_hfx
      external nwxc_has_cam, xc_hfexch, nwxc_wght_hfx

      
      if (nelec0 == 1) then
         veepsi = 0.0d0
         return
      endif

      if (.not. MA_alloc_get(MT_DBL, norb, 'vhfpsi', l_hf, k_hf))
     &     call errquit('WFAT: vhfpsi allocation failed', 1, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, norb, 'vxcpsi', l_xc, k_xc))
     &     call errquit('WFAT: vxcpsi allocation failed', 2, MA_ERR)
      

c     ==== Calculate vhf_psi ====
      chfx = nwxc_wght_hfx()
      rs = nwxc_has_cam()
      if (rs) then
         call nwxc_get_cam(alpha, beta, gamma, srhfx)
         if (srhfx) call errquit('WFAT: Use of range-separated '//
     &        'exchange in WFAT at the moment does not allow '//
     &        'short range HF.', 1, UERR)
      else
         alpha = 1.0d0
         beta = 0.0d0
         gamma = 0.0d0
         srhfx = .false.
      endif
      call wfat_vhfasy_1ptbas(bhan, nbas, nspin, vpt, bvalas, 
     &                        chfx, rs, alpha, beta, gamma, norb, 
     &                        vcf, vcfas, spinid, rdm, rdmas, scr, 
     &                        dbl_mb(k_hf : k_hf+norb-1))
      
      
c     ==== Calculate vxc_psi ====
      call wfat_vxcasy_1ptbas(nbas, nspin, bval, d1bval, d2bval,
     &                        bvalas, norb, vcfas, spinid, rdm, xclvl, 
     &                        dbl_mb(k_xc : k_xc+norb-1))

      
c     ==== Calculate vee_psi ====
      veepsi = dbl_mb(k_hf : k_hf+norb-1) + dbl_mb(k_xc : k_xc+norb-1)


      if (.not. MA_free_heap(l_hf))
     &     call errquit('WFAT: vhfpsi deallocation failed', 3, UERR)
      if (.not. MA_free_heap(l_xc))
     &     call errquit('WFAT: vxcpsi deallocation failed', 4, UERR)

      end

c     Notes:
c       1) d1bval and d2bval are used only when xclvl is equal to 'gga' or 'mgga'.
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_vhfasy_1ptbas(bhan, nbas, nspin, vpt, bvalas, 
     &                              chfx, rs, alpha, beta, gamma, norb,
     &                              vcf, vcfas, spinid, rdm, rdmas, scr,
     &                              vhfpsi)
      
      implicit none
#include "stdio.fh"
#include "global.fh"      
#include "errquit.fh"      
#include "wfat_qm_param.fh"      
      integer bhan                               ! [input] basis handle
      integer nbas                               ! [input] no. of basis corresponding to bhand
      integer nspin                              ! [input] number of spin (1=RHF/unpol., 2=UHF/pol.)
      double precision vpt(3)                    ! [input] evaluation point
      double precision bvalas(nbas)              ! [input] basis functions evaluated at vpt
      double precision chfx                      ! [input] coefficient of exact exchange
      logical rs                                 ! [input] rs=.true. means range-separation is employed
      double precision alpha                     ! [input] alpha coefficient of range-separated exact exchange (only used when rs=.true.)
      double precision beta                      ! [input] beta coefficient of range-separated exact exchange (only used when rs=.true.)
      double precision gamma                     ! [input] attenuation paramaterfor range-separated exact exchange (only used when rs=.true.)
      integer norb                               ! [input] no. of input orbitals
      double precision vcf(nbas, norb)           ! [input] coefficients of input orbitals
      double precision vcfas(nbas, norb)         ! [input] coefficients of input orbitals
      integer spinid(norb)                       ! [input] spin of input orbitals
      double precision rdm(nbas, nbas, nspin)    ! [input] reduced density matrix
      double precision rdmas(nbas, nbas, nspin)  ! [input] reduced density matrix
      double precision scr(nbas, *)              ! [output] scratch array, if rs=.false. then minimum no. of colmuns is 2*nbas+1+nspin, otherwise 3*nbas+1+nspin.
      double precision vhfpsi(1:norb)            ! [output] output orbitals evaluated at vpt

      integer nc, nr, nh, nx, ner
      integer pc0, pcn, pr0, prn, ph0, phn, px0, pxn, per0, pern
      integer io, is
      double precision crg, trhar, BETA_THR, CHFX_THR

      parameter (BETA_THR=1.0d-13, CHFX_THR=1.0d-13)


c     ==== Initial checks ====
      if (minval(spinid) < 1 .or. maxval(spinid) > 2) then
         call errquit('WFAT: Elements of spinid should be either 1 '//
     &                'or 2.', 1, UERR)
      endif

      if (chfx < 0.0d0 .or. chfx > 1.0d0)
     &     call errquit('WFAT: The argument chfx of '//
     &     'wfat_vhfpsi_1ptbas must be between 0 and 1', 1, UERR)
      if ((alpha < 0.0d0 .or. alpha > 1.0d0) .and. rs)
     &     call errquit('WFAT: The argument alpha of '//
     &     'wfat_vhfpsi_1ptbas must be between 0 and 1', 1, UERR)
      if ((beta < 0.0d0 .or. beta > 1.0d0) .and. rs)
     &     call errquit('WFAT: The argument beta of '//
     &     'wfat_vhfpsi_1ptbas must be between 0 and 1', 1, UERR)
      if (gamma < 0.0d0 .and. rs)
     &     call errquit('WFAT: The argument gamma of '//
     &     'wfat_vhfpsi_1ptbas must be positive or zero', 1, UERR)      
      if (alpha+beta > 1.0d0 .and. rs)
     &     call errquit('WFAT: The sum of argument alpha and beta '//
     &     'of wfat_vhfpsi_1ptbas must be between 0 and 1', 1, UERR)
      



c     clmat : nbas x nbas
c     rdm.clmat :  nbas x nbas
      
c     ==== pc0 = pointer to Coulomb matrix ====
      nc = nbas
      pc0 = 1
      pcn = pc0 + nc - 1

c     ==== pr0 = pointer to the resulting matrix ====
c     ====       from rdm.coulomb_matrix         ====
      nr = nbas
      pr0 = pc0 + nc
      prn = pr0 + nr - 1

c     ==== ph0 = pointer to Hartree operator at vpt ====
      nh = 1
      ph0 = pr0 + nr
      phn = ph0 + nh - 1      !! i.e. phn = ph0 always

c     ==== px0 = pointer to HF exchange operator at vpt ====
      nx = nspin
      px0 = ph0 + nh
      pxn = px0 + nx - 1

c     ==== per0 = pointer to erf-modified COulomb matrix ====
      ner = nbas
      per0 = px0 + nx
      pern = per0 + ner - 1
      

      if (nelec0 > 1) then
c        ==== When multi-electron, ... ====
         crg = 1.0d0
         call int_1ehyd_all(bhan, nbas, vpt, crg, scr(:, pc0:pcn)) !! scr = 1e coulomb interaction matrix

         if (rs .and. beta > BETA_THR) then
            call wfat_1ehyd_erf_all(bhan, nbas, gamma, vpt, crg,
     &           scr(:, per0:pern))
         else if (rs .and. beta <= BETA_THR) then
            scr(:, per0:pern) = 0.0d0
         endif
      else if (nelec0 == 1) then
c        ==== When single electron, ... ====
         vhfpsi = 0.0d0
         return
      endif


      
c     ++++++ START: CALCULATE HARTREE OPERATOR AT VPT ++++++
c     ==== Calculate rdm*clmat for Hartree operator ====
      if (nspin == 1) then
         scr(:, pr0:prn) = 2 * matmul(rdm(:,:,1), scr(:, pc0:pcn))
      else if (nspin == 2) then
         scr(:, pr0:prn) = matmul(sum(rdm, 3), scr(:, pc0:pcn))
      endif

c     ==== Calculate e-e Hartree operator at vpt ====
      call wfat_trace(nbas, scr(:, pr0:prn), trhar)
      scr(:, ph0) = trhar * bvalas
c     ++++++ END: CALCULATE HARTREE OPERATOR AT VPT ++++++
      

      
c     ++++++ START: CALCULATE HF EXCHANGE OPERATOR AT VPT ++++++
      do is = 1, nspin
         if (rs) then
            if (nel(is) == 0 .or. chfx*(alpha+beta) < CHFX_THR) then
               scr(:, px0+is-1) = 0.0d0
               cycle
            endif
         else
            if (nel(is) == 0 .or. chfx < CHFX_THR) then
               scr(:, px0+is-1) = 0.0d0
               cycle
            endif
         endif
         
c        ==== Calculate rdm*clmat for HF exchange operator ====
         if (rs) then
ccc            scr(:, pr0:prn) = matmul(rdm(:,:,is),
ccc     &           alpha*scr(:, pc0:pcn) + beta*scr(:, per0:pern))
            scr(:, pr0:prn) = matmul(rdmas(:,:,is),
     &           alpha*scr(:, pc0:pcn) + beta*scr(:, per0:pern))
         else            
            scr(:, pr0:prn) = matmul(rdmas(:,:,is), scr(:, pc0:pcn))
         endif

c        ==== Calculate HF exchange operator at vpt ====
         scr(:, px0+is-1) = - matmul(bvalas, scr(:, pr0:prn))
      enddo
c     ++++++ END: CALCULATE HF EXCHANGE OPERATOR AT VPT ++++++

      

c     ==== Calculate vhfpsi ====
      do io = 1, norb
         is = spinid(io)

         vhfpsi(io) = sum(scr(:, ph0) * vcfas(:,io)) +
     &                sum(chfx * scr(:, px0+is-1) * vcf(:,io))
      enddo
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_vxcasy_1ptbas(nbas, nspin, bval, d1bval, d2bval, 
     &                              bvalas, norb, vcfas, spinid, rdm, 
     &                              xclvl, vxcpsi)


      implicit none
#include "stdio.fh"            
#include "wfat_indexing.fh"      
#include "errquit.fh"
#include "wfat_qm_param.fh"      
      integer nbas
      integer nspin
      double precision bval(nbas)             ! [input], value of the bases at the evaluation point.
      double precision d1bval(3, nbas)        ! [input], value of the 1st derivative bases at the evaluation point.    1)
      double precision d2bval(6, nbas)        ! [input], value of the 2nd derivative bases at the evaluation point.    1)
      double precision bvalas(nbas)           ! [input], value of the bases at a point at the asymptotic sphere in the same radial line as the evaluation point of bval
      integer norb
      double precision vcfas(nbas, norb)      ! [input], asymptotically modified coefficients of the channels.
      integer spinid(norb)
      double precision rdm(nbas, nbas, nspin)
      character*4 xclvl        ! 1)
      double precision vxcpsi(norb)

      integer maxd, io, is
      double precision dens(nspin), d1dens(3,nspin), d2dens(6,nspin)
      double precision gamma(3), d1gamma(3,3), tau(nspin), lapl(nspin)
      double precision vxc(nspin), fn
      double precision dfndr(nspin), dfndr2(NCOL_AMAT2)
      double precision dfndg(3), dfndg2(NCOL_CMAT2)
      double precision dfndt(nspin), dfndt2(NCOL_MMAT2)
      


      if (xclvl /= 'lda' .and. xclvl /= 'gga' .and. xclvl /= 'mgga'
     &    .and. xclvl /= 'none')
     &     call errquit('WFAT: Argument xclvl of wfat_vxcpsi_1ptbas '//
     &                  'must be either lda, gga, mgga, or none.',
     &                  1, UERR)

      if (xclvl == 'none') then
         vxcpsi = 0.0d0
         return
      endif
      
      
      if (nelec0 > 1) then
c        ==== When multi-electron, ... ====
         if (xclvl == 'lda') maxd = 0
         if (xclvl == 'gga' .or. xclvl == 'mgga') maxd = 2
         call wfat_dens_1ptbas(nbas, nspin, rdm, maxd,
     &                         bval, d1bval, d2bval,
     &                         dens, d1dens, d2dens)

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
     &                        dfndg, dfndt)

c           ==== Calculate LDA potential ====
            call wfat_vxclda_1pt(nspin, dfndr, vxc)
            
         else if (xclvl == 'gga' .or. xclvl == 'mgga') then
            call wfat_gamma_1pt(nspin, d1dens, d2dens, gamma, d1gamma)
ccheck            if (xclvl == 'mgga') 
ccheck     &           call wfat_tau_1pt(nbas, nspin, rdm, d1bval, tau)
            call wfat_tau_1pt(nbas, nspin, rdm, d1bval, tau)
            
            call nwxc_eval_df2(nspin, 1, dens, gamma, tau, fn, dfndr, 
     &                         dfndr2, dfndg, dfndg2, dfndt, dfndt2)

            call wfat_lapl_1pt(nbas, nspin, d2dens, lapl)

            
c           ==== Calculate GGA or mGGA potential ====
            if (xclvl == 'gga') then
               call wfat_vxcgga_1pt(nspin, d1dens, lapl, gamma, d1gamma,
     &                              dfndr, dfndg, dfndg2, vxc)
            else if(xclvl == 'mgga') then
               call wfat_vxcmgga_1pt()
            endif
            
         endif


         do io = 1, norb
            is = spinid(io)            
c            moas = sum(vcf(:,io) * bvalas)
c            rr = sqrt(sum((vpt-asc)**2))
c            kappa = sqrt(2*abs(erg(io)))
c            vxcpsi(io) = vxc(is) * moas * exp(-kappa*(rr-ras))
            vxcpsi(io) = vxc(is) * sum(vcfas(:,io) * bvalas)
         enddo
      else
c        ==== When single electron, ... ====
         vxcpsi = 0.0d0
         return
      endif

      end
c     Notes:
c       1) d1bval and d2bval are used only when xclvl is equal to 'gga' or 'mgga'.
c=======================================================
