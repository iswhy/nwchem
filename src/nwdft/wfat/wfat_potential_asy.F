c=======================================================
c=======================================================
      subroutine wfat_dft_vcasy_1pt(ghan, bhan, nbas, nspin,
     &           nch, vch, erg, spinid, chpos, vpt, asc, ras,
     &           rdm, rdmas, zi, smh, scr, vc_psi)

      implicit none
#include "stdio.fh"      
#include "errquit.fh"      
cOLD #include "cdft.fh"
      integer ghan
      integer bhan
      integer nbas
      integer nspin
      integer nch
      double precision vch(nbas, nch)    ! [input], nbasis*nch
      double precision erg(nch)
      integer spinid(nch)
      double precision chpos(3, nch)
      double precision vpt(3)
      double precision asc(3)            ! [input], coordinate of asymptotic center.
      double precision ras
      double precision rdm(nbas, nbas, nspin)
      double precision rdmas(nbas, nbas, nspin)
      double precision zi
      double precision smh
      double precision scr(1:*) ! [input], scratch array with min. no. of elements
c                                 being 20*nbas+4*nspin + nbas * (3*nbas+1+nspin) + 4*nch + nbas + nbas*nch
      double precision vc_psi(nch)

      character*4 xclvl
      integer nb0, nb1, nb2, nb3, nd0, nd1, nsc1, nsc2, nve, nvn, nvs
      integer nba, ncha
      integer pb0, pb1, pb2, pb3, pd0, pd1, psc1, psc2, pve, pvn, pvs
      integer pba, pcha
      integer io, io1, io2
      integer dbas
      double precision vn_new, rc(3), bb(3), dummy1, dummy2, dummy3, rr
      double precision kappa

      logical xc_gotxc, xc_chkgrad, xc_chktau, xc_hascam, xc_hassrhf
      external xc_gotxc, xc_chkgrad, xc_chktau, xc_hascam, xc_hassrhf

      
      
      if (xc_chktau()) then
         call errquit('WFAT: At the moment, mGGA functional is not '//
     &                'yet supported in WFAT.', 11, UERR)
      endif

      dbas = 0
      xclvl = 'none'
      if (xc_gotxc()) then
         dbas = 0
         xclvl = 'lda'
      endif
      if (xc_chkgrad()) then
         dbas = 2
         xclvl = 'gga'
      endif      
      if (xc_chktau()) then
         dbas = 2
         xclvl = 'mgga'
      endif
      if (xc_hascam()) then
cOLD         call errquit('WFAT: At the moment, WFAT cannot be used with '//
cOLD     &                'CAM functionals.', 1, UERR)
      endif
      if (xc_hassrhf()) then
         call errquit('WFAT: At the moment, WFAT cannot be used with '//
     &                'short-range Hartree-Fock functionals.', 1, UERR)
      endif


c     ==== Assign indexes for scratch array ====
      nb0 = nbas
      nb1 =  3 * nbas
      nb2 =  6 * nbas
      nb3 = 10 * nbas
      nd0 = nspin
      nd1 = 3 * nspin
cOLD      nsc1 = nbas * (2*nbas+1+nspin)
      nsc1 = nbas * (3*nbas+1+nspin)      
      nsc2 = nch
      nve = nch
      nvn = nch
      nvs = nch
      nba = nbas
      ncha = nbas * nch
      
      pb0 = 1             !! basis functions
      pb1 = pb0 + nb0     !! 1st derivative of basis functions
      pb2 = pb1 + nb1     !! 2nd derivative of basis functions
      pb3 = pb2 + nb2     !! 3rd derivative of basis functions
      pd0 = pb3 + nb3     !! density functions
      pd1 = pd0 + nd0     !! 1st derivative of density functions
      psc1 = pd1 + nd1    !! scratch array for wfat_dft_veepsi_1ptbas
      psc2 = psc1 + nsc1  !! scratch array for wfat_vnucpsi_1ptbas
      pve = psc2 + nsc2   !! vee|psi> array
      pvn = pve + nve     !! vnuc|psi> array
      pvs = pvn + nvn     !! vs|psi> array
      pba = pvs + nvs     !! basis functions at the asymptotic sphere
      pcha = pba + nba    !! asymptotically modified coefficients

      
c     ==== Calculate basis and its derivatives ====
      call wfat_eval_basis_1pt(ghan, bhan, dbas, vpt,
     &                         scr(pb0 : pb0+nb0-1),
     &                         scr(pb1 : pb1+nb1-1),
     &                         scr(pb2 : pb2+nb2-1),
     &                         scr(pb3 : pb3+nb3-1))     !! 1)

      bb = vpt - asc
      rc = asc + bb / sqrt(sum(bb*bb)) * ras
      call wfat_eval_basis_1pt(ghan, bhan, 0, rc,
     &                         scr(pba : pba+nba-1),
     &                         dummy1, dummy2, dummy3)
      
      
c     ==== Calculate the asymptotically modified coefficients and RDM ====
      do io = 1, nch
         rr = sqrt(sum((vpt-asc)**2))
         kappa = sqrt(2*abs(erg(io)))
         io1 = pcha + (io-1)*nbas
         io2 = io1 + nbas - 1
         scr(io1:io2) = vch(:,io) * exp(-kappa*(rr-ras))
         if (rr-ras < 0.0d0)
     &        call errquit('WFAT: Evaluation point is still not '//
     &        'within asymptotic criteria', 1, UERR)
      enddo
      

c     ==== Calculate vee_psi ====
      call wfat_dft_veeasy_1ptbas(bhan, nbas, nspin, vpt,
     &     scr(pb0 : pb0+nb0-1),
     &     scr(pb1 : pb1+nb1-1),
     &     scr(pb2 : pb2+nb2-1),
     &     scr(pba : pba+nba-1), 
     &     nch,
     &     vch, scr(pcha : pcha+ncha-1),
     &     spinid, rdm, rdmas, xclvl,
     &     scr(psc1 : psc1+nsc1-1),
     &     scr(pve : pve+nve-1))


c     ==== Calculate vn_psi ====
      call wfat_vnucpsi_1ptbas(ghan, smh, vpt, nbas,
     &     scr(pba : pba+nba-1), nch, scr(pcha : pcha+ncha-1),
     &     scr(psc2 : psc2+nsc2-1),
     &     scr(pvn : pvn+nvn-1))

      
c     ==== Calculate vs_psi ====
cchecklater      zi = 1.0d0
      do io = 1, nch
         io1 = pcha + (io-1)*nbas
         io2 = io1 + nbas - 1
         call wfat_vhpsi_1ptbas(chpos(:,io), zi, smh, vpt, nbas,
     &        scr(pba : pba+nba-1), 1, scr(io1:io2),
     &        scr(pvs+io-1 : pvs+io-1))
      enddo


c     ==== Calculate vc_psi ====
      vc_psi = -(scr(pve : pve+nve-1) +
     &           scr(pvn : pvn+nvn-1) +
     &           scr(pvs : pvs+nvs-1))      


      end
c Notes:      
c   1) The 3rd (pb3) derivative is currently not used.
c=======================================================
