c=======================================================
c=======================================================
      subroutine wfat_dft_vcpsi_1pt_c(icent, ghan, bhan, nbas, nspin,
     &           nch, wfat1, vch, spinid, chpos, vpt, cgrid, rdm, zi, 
     &           smh, scr, mask, ergmask, amsk, bmsk, vc_psi)
      
c     DESCRIPTION:
c     When wfat is false, vch, spinid, chpos, and vc_psi must have nch columns
c     containing the GTO coefficients of the unperturbed orbital. Otherwise, they must
c     have 3*nch more columns (thus 4*nch columns in total) containing the Cartesian
c     components of the correction of each orbital arranged in this way,
c       orb1_x orb1_y orb1_z orb2_x orb2_y orb2_z ... orb<nch>_x orb<nch>_y orb<nch>_z.
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer icent
      integer ghan
      integer bhan
      integer nbas
      integer nspin
      integer nch                     ! [input], the number of unperturbed channels.
      logical wfat1
      double precision vch(nbas,*)    ! [input], no. of column is nch if wfat1 is false, otherwise, 4*nch.
      integer spinid(*)               ! [input], no. of elements is nch if wfat1 is false, otherwise, 4*nch.
      double precision chpos(3,*)     ! [input], no. of column is nch if wfat1 is false, otherwise, 4*nch.
      double precision vpt(3, nch)    ! [input], Cartesian coordinate of the evaluation point in the molecular frame.
      logical cgrid(nch)
      double precision rdm(nbas, nbas, nspin)
      double precision zi
      double precision smh
      double precision scr(1:*)    ! [input], scratch array with min. no. of elements being 21*nbas+4*nspin + nbas * (3*nbas+1+nspin) + 4*nch
      logical mask
      double precision ergmask      !< [Input]. Only used when mask=.true.
      double precision amsk         !< [Input]. Only used when mask=.true.
      double precision bmsk         !< [Input]. Only used when mask=.true.
      double precision vc_psi(*)    ! [output], no. of elements is nch if wfat1 is false, otherwise, 4*nch.

      integer i, nch1, l_cf, k_cf, i_cf, b, i1, i2
      double precision vcp(4)

      
      nch1 = 0
      b = 0
      if (wfat1) then
         nch1 = 3*nch
         b = 3
      endif
      
      if (icent /= -1) then
c        When icent /= -1, the columns of vpt should all be identical. In this
c        case, the first column will be used for the following calculation.
         call wfat_dft_vcpsi_1pt(ghan, bhan, nbas, nspin, nch+nch1, vch,
     &        spinid, chpos, vpt(:,1), rdm, zi, smh, scr, mask, ergmask,
     &        amsk, bmsk, vc_psi)
      else
         if (.not. MA_push_get(MT_DBL, nbas*(1+b), 'vc', l_cf, k_cf))
     &        call errquit('WFAT: vc allocation failed', 24, MA_ERR)

         do i = 1, nch
            if (wfat1) then
               i1 = 1 + nch + (i-1)*3
               i2 = i1 + 2
            endif
            
            if (cgrid(i)) then
               dbl_mb(k_cf:k_cf+nbas-1) = vch(:,i)
               if (wfat1) then
                  i_cf = k_cf + nbas
                  call dcopy(nbas*3, vch(1,i1), 1, dbl_mb(i_cf), 1)
               endif
               call wfat_dft_vcpsi_1pt(ghan, bhan, nbas, nspin, 1+b, 
     &              dbl_mb(k_cf), spinid(i), chpos(:,i), vpt(:,i), rdm, 
     &              zi, smh, scr, mask, ergmask, amsk, bmsk, vcp)
               vc_psi(i) = vcp(1)
               if (wfat1) vc_psi(i1:i2) = vcp(2:4)
            else
               vc_psi(i) = 0.0d0
               if (wfat1) vc_psi(i1:i2) = 0.0d0
            endif
         enddo

         if (.not. MA_pop_stack(l_cf))
     &        call errquit('WFAT: vc deallocation failed', 24, GA_ERR)
      endif
      

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dft_vcpsi_1pt(ghan, bhan, nbas, nspin,
     &           nch, vch, spinid, chpos, vpt, rdm, zi, smh,
     &           scr, mask, ergmask, amsk, bmsk, vc_psi)


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
      integer spinid(nch)
      double precision chpos(3, nch)
      double precision vpt(3)      ! [input], Cartesian coordinate of the evaluation point in the molecular frame.
      double precision rdm(nbas, nbas, nspin)
      double precision zi
      double precision smh
      double precision scr(1:*)    ! [input], scratch array with min. no. of elements being 21*nbas+4*nspin + nbas * (3*nbas+1+nspin) + 4*nch
      logical mask
      double precision ergmask      !< [Input]. Only used when mask=.true.
      double precision amsk         !< [Input]. Only used when mask=.true.
      double precision bmsk         !< [Input]. Only used when mask=.true.
      double precision vc_psi(nch)

      character*4 xclvl
      integer nb0, nb1, nb2, nb3, nd0, nd1, nsc1, nsc2, nve, nvn, nvs
      integer nrf
      integer pb0, pb1, pb2, pb3, pd0, pd1, psc1, psc2, pve, pvn, pvs
      integer prf
      integer io, i, ib
      integer dbas

      double precision vn_new, kp, maskf

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
      nrf = nbas

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
      prf = pvs + nvs     !! radial part of basis functions

      
c     ==== Calculate basis and its derivatives ====
      call wfat_eval_basis_radf(ghan, bhan, dbas, vpt,
     &                         scr(prf : prf+nrf-1),
     &                         scr(pb0 : pb0+nb0-1),
     &                         scr(pb1 : pb1+nb1-1),
     &                         scr(pb2 : pb2+nb2-1),
     &                         scr(pb3 : pb3+nb3-1))     !! 1)

      if (mask) then
         do ib = 1, nbas
            kp = sqrt(2*abs(ergmask))
            call wfat_mask_exp(1, vpt, chpos(:,1), zi, kp, amsk, bmsk,
     &           scr(prf+ib-1), maskf)
            scr(pb0+ib-1) = scr(pb0+ib-1) * maskf
         enddo
      endif
         
      
c     ==== Calculate vee_psi ====
      call wfat_dft_veepsi_1ptbas(bhan, nbas, nspin, vpt,
     &     scr(pb0 : pb0+nb0-1),
     &     scr(pb1 : pb1+nb1-1),
     &     scr(pb2 : pb2+nb2-1),
     &     nch, vch, spinid, rdm, xclvl,
     &     scr(psc1 : psc1+nsc1-1),
     &     scr(pve : pve+nve-1))


c     ==== Calculate vn_psi ====
      call wfat_vnucpsi_1ptbas(ghan, smh, vpt, nbas,
     &     scr(pb0 : pb0+nb0-1), nch, vch,
     &     scr(psc2 : psc2+nsc2-1),
     &     scr(pvn : pvn+nvn-1))

      
c     ==== Calculate vs_psi ====
cchecklater      zi = 1.0d0
      do io = 1, nch
         call wfat_vhpsi_1ptbas(chpos(:,io), zi, smh, vpt, nbas,
     &        scr(pb0 : pb0+nb0-1), 1, vch(:,io:io),
     &        scr(pvs+io-1 : pvs+io-1))
      enddo


c     ==== Calculate vc_psi ====
c      vc_psi = scr(pve : pve+nve-1) +
c     &         scr(pvn : pvn+nvn-1) -
c     &         scr(pvs : pvs+nvs-1)      

      vc_psi = -(scr(pve : pve+nve-1) +
     &           scr(pvn : pvn+nvn-1) +
     &           scr(pvs : pvs+nvs-1))
      
      end
c Notes:      
c   1) The 3rd (pb3) derivative is currently not used.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dft_vcpsi_getcoef(rtdb, geom, bas, nb, ns, rdm,
     &           nv, vin, sid, rc, vnucvin, vhfvin, vxcvin, vsvin)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "wfat_qm_param.fh"
      integer rtdb                    ! [input] 
      integer geom                    ! [input]
      integer bas                     ! [input]
      integer nb                      ! [input]
      integer ns                      ! [input]
      double precision rdm(nb,nb,ns)  ! [input]
      integer nv                      ! [input]
      double precision vin(nb,nv)     ! [input]
      integer sid(nv)                 ! [input]
      double precision rc(3,nv)       ! [input]
      double precision vnucvin(nb,nv) ! [output]
      double precision vhfvin(nb,nv)  ! [output]
      double precision vxcvin(nb,nv)  ! [output]
      double precision vsvin(nb,nv)   ! [output]
      
      character*16 atomname
      character*9 rdm_name(2)
      integer g_rd(2)
      integer is, i, ia, ndid
      double precision chfx, numnel, exc(2)

      double precision xc_hfexch
      external xc_hfexch


cforcheck      call wfat_test_fock(rtdb, geom, bas, nb, ns, nel, g_mov, 
cforcheck     &     charge, multip, tol2e, oskel, ixcfit, iolgc)
      

      rdm_name = (/ 'rdm_alpha', 'rdm_beta' /)
      ndid = ga_nodeid()
      

c     ==== vnuc|vin> Mulliken populations ====
      call wfat_vnucvec_coef(bas, nb, nv, vin, vnucvin)


c     ==== Create RDM GA ====
      do is = 1, ns
         if (.not. ga_create(mt_dbl, nb, nb, rdm_name(is), 0, 0,
     &        g_rd(is)))
     &        call errquit('WFAT: Cannot create '//rdm_name(is)//
     &        ' global array.', 0, GA_ERR)
         call ga_zero(g_rd(is))

         if (ndid == 0) then
            call ga_put(g_rd(is), 1, nb, 1, nb, rdm(:,:,is), nb)
         endif
      enddo
      call ga_sync()


c     ==== (J + cexch*K)|vin> Mulliken populations ====      
      chfx = xc_hfexch()      
      call wfat_vhfvec_coef(geom, bas, nb, ns, g_rd(1:ns), nv, vin,
     &                      sid, chfx, tol2e, oskel, vhfvin)

  
c     ==== Vxc|vin> Mulliken populations ====
      call wfat_vxcvec_coef(rtdb, geom, bas, nb, ns, g_rd, nv, vin, sid,
     &                      ixcfit, iolgc, numnel, exc, vxcvin)           
cforcheck      if (ndid==0) then
cforcheck         write(luout, '(10x, "Numerically integrated no. of '//
cforcheck     &         'electrons = ", f14.6)') numnel
cforcheck         write(luout, '(10x, "Exchange energy = ", es14.6)') exc(1)
cforcheck         write(luout, '(10x, "Correlation energy = ", es14.6)') exc(2)
cforcheck      endif
      call ga_sync()

      
c     ==== Destroy RDM GA ====
      do is = 1, ns
         if (.not. ga_destroy(g_rd(is)))
     &        call errquit('WFAT: Cannot destroy '//rdm_name(is)//'.',
     &        1, GA_ERR)
      enddo

      
c     ==== V_s|vin> ====
      do i = 1, nv
         call wfat_vhydvec_coef(bas, nb, rc(:,i), 1.0d0, 1, vin(:,i),
     &        vsvin(:,i))
      enddo
      

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dft_vcpsi_fromcoef(ghan, bhan, nbas, nch,
     &           vnucpsi, vhfpsi, vxcpsi, vspsi, vpt, scr, vc_psi)

      implicit none
#include "stdio.fh"      
#include "errquit.fh"
      integer ghan
      integer bhan
      integer nbas
      integer nch
      double precision vnucpsi(nbas,nch)
      double precision vhfpsi(nbas,nch)
      double precision vxcpsi(nbas,nch)
      double precision vspsi(nbas,nch)
      double precision vpt(3)
      double precision scr(1:*)
      double precision vc_psi(nch)

      integer nb0, nb1, nb2, nb3, nvn, nve, nvs
      integer pb0, pb1, pb2, pb3, pvn, pve, pvs
      integer io
      integer dbas
      
      
c     ==== Assign indexes for scratch array ====
      dbas = 0
      
      nb0 = nbas
      nb1 =  3 * nbas
      nb2 =  6 * nbas
      nb3 = 10 * nbas
      nve = nch
      nvn = nch
      nvs = nch
      
      pb0 = 1             !! basis functions
      pb1 = pb0 + nb0     !! 1st derivative of basis functions
      pb2 = pb1 + nb1     !! 2nd derivative of basis functions
      pb3 = pb2 + nb2     !! 3rd derivative of basis functions
      pve = pb3 + nb3     !! vee|psi> array
      pvn = pve + nve     !! vnuc|psi> array
      pvs = pvn + nvn     !! vs|psi> array 

      
c     ==== Calculate basis and its derivatives ====
      call errquit('Think of whether better using '//
     &     'wfat_eval_basis_radf instead of wfat_eval_basis_1pt below',
     &     46, UERR)
      call wfat_eval_basis_1pt(ghan, bhan, dbas, vpt,
     &                         scr(pb0 : pb0+nb0-1),
     &                         scr(pb1 : pb1+nb1-1),
     &                         scr(pb2 : pb2+nb2-1),
     &                         scr(pb3 : pb3+nb3-1))     !! 1)


c     ==== Calculate vc_psi ====      
      do io = 1, nch
         scr(pvn+io-1) = dot_product(scr(pb0 : pb0+nb0-1),
     &                   vnucpsi(:,io))
         scr(pve+io-1) = dot_product(scr(pb0 : pb0+nb0-1),
     &                   vhfpsi(:,io) + vxcpsi(:,io))
         scr(pvs+io-1) = dot_product(scr(pb0 : pb0+nb0-1),
     &                   vspsi(:,io))

         vc_psi(io) = scr(pvn+io-1) + scr(pve+io-1) - scr(pvs+io-1)
      enddo


c      write(luout, '("        sum vn_psi new = ", e14.5)')
c     &     sum(scr(pvn : pvn+nvn-1))
c      write(luout, '("        sum vee_psi new = ", e14.5)')
c     &     sum(scr(pve : pve+nve-1))
c      write(luout, '("        sum vs_psi new = ", e14.5)')
c     &     sum(scr(pvs : pvs+nvs-1))
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_xclvl(dr, lv)

      implicit none
#include "errquit.fh"
      integer dr           ! [output], max. derivative order for the basis functions needed.
      character*(*) lv     ! [output], XC level.
      
      logical xc_gotxc, xc_chkgrad, xc_chktau, xc_hascam, xc_hassrhf
      external xc_gotxc, xc_chkgrad, xc_chktau, xc_hascam, xc_hassrhf
      

      dr = 0
      lv = 'none'
      if (xc_gotxc()) then
         dr = 0
         lv = 'lda'
      endif
      if (xc_chkgrad()) then
         dr = 2
         lv = 'gga'
      endif      
      if (xc_chktau()) then
         dr = 2
         lv = 'mgga'
      endif
      if (xc_hascam()) then
cOLD         call errquit('WFAT: At the moment, WFAT cannot be used with '//
cOLD     &                'CAM functionals.', 1, UERR)
      endif
      if (xc_hassrhf()) then
         call errquit('WFAT: At the moment, WFAT cannot be used with '//
     &                'short-range Hartree-Fock functionals.', 1, UERR)
      endif


      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_collect_vcpsi(geom, bas, nb, ns, nch, nptp, vch,
     &     spinid, chpos, vpt, rdm, zi, smh, ismask, ergmask, amsk, 
     &     bmsk, vc_psi)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer bas
      integer nb
      integer ns
      integer nch
      integer nptp
      double precision vch(nb,nch)
      integer spinid(nch)
      double precision chpos(3,nch)
      double precision vpt(3,nptp)
      double precision rdm(nb,nb,ns)
      double precision zi
      double precision smh
      logical ismask
      double precision ergmask
      double precision amsk         !< [Input]. Only used when mask=.true.
      double precision bmsk         !< [Input]. Only used when mask=.true.
      double precision vc_psi(nptp,nch)

      integer nscr, i, l_sc, k_sc, l_vc, k_vc

      nscr = 21*nb + 4*ns + nb * (3*nb+1+ns) + 4*nch
      if (.not. MA_alloc_get(MT_DBL, nscr, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nch, 'vc_psi', l_vc, k_vc))
     &     call errquit('WFAT: vc_psi allocation failed', 59, MA_ERR)
      
      do i = 1, nptp
         call wfat_dft_vcpsi_1pt(geom, bas, nb, ns, nch, vch, spinid, 
     &           chpos, vpt(:,i), rdm, zi, smh, dbl_mb(k_sc), ismask,
     &           ergmask, amsk, bmsk, dbl_mb(k_vc))
         vc_psi(i,:) = dbl_mb(k_vc:k_vc+nch-1)
      enddo

      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: scratch deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc_psi deallocation failed', 61, MA_ERR)

      end
c=======================================================
      
