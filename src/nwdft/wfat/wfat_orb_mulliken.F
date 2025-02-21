c=======================================================
c=======================================================
      subroutine wfat_mulliken(rtdb, geom, bas, nb, ns, rdm, nv, vin,
     &                         sid, rc)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"      
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
cOLD      double precision vnucvin(nb,nv) ! [output]
cOLD      double precision vhfvin(nb,nv)  ! [output]
cOLD      double precision vxcvin(nb,nv)  ! [output]
cOLD      double precision vsvin(nb,nv)   ! [output]
      
      character*16 atomname
      character*9 rdm_name(2)
      integer l_ov, l_pe, l_ve, l_vxc, l_hy
      integer k_ov, k_pe, k_ve, k_vxc, k_hy
      integer i_ve, i_vxc
      integer l_inpop, l_nupop, l_hfpop, l_xcpop, l_efpop
      integer l_vspop, l_vcpop
      integer k_inpop, k_nupop, k_hfpop, k_xcpop, k_efpop
      integer k_vspop, k_vcpop
      integer i_inpop, i_nupop, i_hfpop, i_xcpop, i_efpop
      integer i_vspop, i_vcpop
      integer g_rd(2)
      integer ninpop, nnupop, nhfpop, nxcpop, nefpop
      integer nvspop, nvc, nvcpop
      integer nov, npe, nve, nvxc, nhy
      integer is, i, ia, ndid, nat
      double precision chfx, numnel, exc(2)

      double precision xc_hfexch
      external xc_hfexch


cforcheck      call wfat_test_fock(rtdb, geom, bas, nb, ns, nel, g_mov, 
cforcheck     &     charge, multip, tol2e, oskel, ixcfit, iolgc)
      

      rdm_name = (/ 'rdm_alpha', 'rdm_beta' /)
      ndid = ga_nodeid()
      

      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 51, GEOM_ERR)
      

c     ==== |vin> Mulliken populations ====
      ninpop = nat * nv
      if (.not. ma_alloc_get(mt_dbl, ninpop, 'vin_pop',
     &     l_inpop, k_inpop))
     &     call errquit('WFAT: Cannot allocate vin_pop.', 1, MA_ERR)
      nov = nb * nb
      if (.not. ma_alloc_get(mt_dbl, nov, 'overlap', l_ov, k_ov))
     &     call errquit('WFAT: Cannot allocate overlap.', 1, MA_ERR)
      call int_1eov_all(bas, nb, dbl_mb(k_ov))
      call wfat_expval_mulliken(nv, bas, dbl_mb(k_ov), vin, geom,
     &     dbl_mb(k_inpop))
cOLD      call wfat_expval_mulliken(nv, bas, vin, geom, dbl_mb(k_inpop))
      

c     ==== vnuc|vin> Mulliken populations ====
      nnupop = nat * nv
      if (.not. ma_alloc_get(mt_dbl, nnupop, 'vnucvin_pop',
     &     l_nupop, k_nupop))
     &     call errquit('WFAT: Cannot allocate vnucvin_pop.', 1, MA_ERR)
      npe = nb * nb
      if (.not. ma_alloc_get(MT_DBL, npe, 'vnuc', l_pe, k_pe))
     &     call errquit('WFAT: Cannot allocate vnuc in '//
     &     'wfat_vnucvec_coef.', 1, MA_ERR)
      call int_1epe_all(bas, nb, dbl_mb(k_pe))
      call wfat_expval_mulliken(nv, bas, dbl_mb(k_pe), vin, geom,
     &     dbl_mb(k_nupop))
cOLD      call wfat_vnucvec_coef(bas, nb, nv, vin, vnucvin)
cOLD      call wfat_expval_mulliken(nv, bas, vnucvin, geom, dbl_mb(k_nupop))


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
      nhfpop = nat * nv
      if (.not. ma_alloc_get(mt_dbl, nhfpop, 'vhfvin_pop',
     &     l_hfpop, k_hfpop))
     &     call errquit('WFAT: Cannot allocate vhfvin_pop.', 1, MA_ERR)
      nve = nb * nb * ns
      if (.not. ma_alloc_get(mt_dbl, nve, 'HF_vee', l_ve, k_ve))
     &     call errquit('WFAT: Cannot allocate HF_vee.', 1, MA_ERR)
      call wfat_get_vhfmat(geom, bas, nb, ns, g_rd, chfx, tol2e, oskel,
     &     dbl_mb(k_ve))
      do i = 1, nv
         is = sid(i)
         i_ve = k_ve + (is-1)*nb**2
         i_hfpop = k_hfpop + (i-1)*nat
         call wfat_expval_mulliken(1, bas, dbl_mb(i_ve), vin(:,i), geom,
     &        dbl_mb(i_hfpop))
      enddo
cOLD      call wfat_vhfvec_coef(geom, bas, nb, ns, g_rd(1:ns), nv, vin,
cOLD     &                      sid, chfx, tol2e, oskel, vhfvin)
cOLD      call wfat_expval_mulliken(nv, bas, vhfvin, geom, dbl_mb(k_hfpop))

      
c     ==== Vxc|vin> Mulliken populations ====
      nxcpop = nat * nv
      if (.not. ma_alloc_get(mt_dbl, nxcpop, 'vxcvin_pop',
     &     l_xcpop, k_xcpop))
     &     call errquit('WFAT: Cannot allocate vxcvin_pop.', 1, MA_ERR)
      nvxc = nb * nb * ns
      if (.not. ma_alloc_get(mt_dbl, nvxc, 'XC_vee', l_vxc, k_vxc))
     &     call errquit('WFAT: Cannot allocate XC_vee.', 1, MA_ERR)
      call wfat_get_vxcmat(rtdb, geom, bas, nb, ns, g_rd, ixcfit, iolgc,
     &     numnel, exc, dbl_mb(k_vxc))
      do i = 1, nv
         is = sid(i)
         i_vxc = k_vxc + (is-1)*nb**2
         i_xcpop = k_xcpop + (i-1)*nat
         call wfat_expval_mulliken(1, bas, dbl_mb(i_vxc), vin(:,i), 
     &        geom, dbl_mb(i_xcpop))
      enddo
cOLD      call wfat_vxcvec_coef(rtdb, geom, bas, nb, ns, g_rd, nv, vin, sid,
cOLD     &                      ixcfit, iolgc, numnel, exc, vxcvin)     
cOLD      call wfat_expval_mulliken(nv, bas, vxcvin, geom, dbl_mb(k_xcpop))
      
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

      
c     ==== Veff|vin> Mulliken populations ====
cOLD      nef = nb * nv
cOLD      if (.not. ma_alloc_get(mt_dbl, nef, 'veffvin', l_ef, k_ef))
cOLD     &     call errquit('WFAT: Cannot allocate veffvin.', 1, MA_ERR)
      nefpop = nat * nv
      if (.not. ma_alloc_get(mt_dbl, nefpop, 'veffvin_pop',
     &     l_efpop, k_efpop))
     &     call errquit('WFAT: Cannot allocate veffvin_pop.', 1, MA_ERR)
cOLD      do i = 1, nv
cOLD         i_ef = k_ef + (i-1)*nb
cOLD         dbl_mb(i_ef:i_ef+nb-1) = vnucvin(:,i) + vhfvin(:,i) +
cOLD     &                            vxcvin(:,i)
cOLD      enddo
cOLD      call wfat_expval_mulliken(nv, bas, dbl_mb(k_ef), geom,
cOLD     &                       dbl_mb(k_efpop))
      dbl_mb(k_efpop : k_efpop+nefpop-1) =
     &     dbl_mb(k_nupop : k_nupop+nnupop-1) +
     &     dbl_mb(k_hfpop : k_hfpop+nhfpop-1) +
     &     dbl_mb(k_xcpop : k_xcpop+nxcpop-1)
      
      
c     ==== V_s|vin> ====
      nvspop = nat * nv
      if (.not. ma_alloc_get(mt_dbl, nvspop, 'vsvin_pop',
     &     l_vspop, k_vspop))
     &     call errquit('WFAT: Cannot allocate vsvin_pop.', 1, MA_ERR)
      nhy = nb * nb
      if (.not. ma_alloc_get(MT_DBL, nhy, 'vhyd', l_hy, k_hy))
     &     call errquit('WFAT: Cannot allocate vhyd in '//
     &     'wfat_vhydvec_coef.', 1, MA_ERR)
      do i = 1, nv
         call int_1ehyd_all(bas, nb, rc(:,i), 1.0d0, dbl_mb(k_hy))
         i_vspop = k_vspop + (i-1)*nat
         call wfat_expval_mulliken(1, bas, dbl_mb(k_hy), vin(:,i), geom,
     &        dbl_mb(i_vspop))
      enddo
cOLD      do i = 1, nv
cOLD         call wfat_vhydvec_coef(bas, nb, rc(:,i), 1.0d0, 1, vin(:,i),
cOLD     &        vsvin(:,i))
cOLD      enddo
cOLD      call wfat_expval_mulliken(nv, bas, vsvin, geom, dbl_mb(k_vspop))
      

c     ==== V_c|vin> ====
cOLD      nvc = nb * nv
cOLD      if (.not. ma_alloc_get(mt_dbl, nef, 'vcvin', l_vc, k_vc))
cOLD     &     call errquit('WFAT: Cannot allocate vcvin.', 1, MA_ERR)
      nvcpop = nat * nv
      if (.not. ma_alloc_get(mt_dbl, nvcpop, 'vcvin_pop',
     &     l_vcpop, k_vcpop))
     &     call errquit('WFAT: Cannot allocate vcvin_pop.', 1, MA_ERR)
cOLD      do i = 1, nv
cOLD         i_vc = k_vc + (i-1)*nb
cOLD         i_ef = k_ef + (i-1)*nb
cOLD         dbl_mb(i_vc:i_vc+nb-1) = dbl_mb(i_ef:i_ef+nb-1) - vsvin(:,i)
cOLD      enddo
cOLD      call wfat_expval_mulliken(nv, bas, dbl_mb(k_vc), geom,
cOLD     &     dbl_mb(k_vcpop))
      dbl_mb(k_vcpop : k_vcpop+nvcpop-1) =
     &     dbl_mb(k_efpop : k_efpop+nefpop-1) -
     &     dbl_mb(k_vspop : k_vspop+nvspop-1)
      
      
c     ==== Print Mulliken populations of ionized channels (orbitals) ====      
      if (ndid==0) then
         write(luout, '(2/)')
         call util_print_centered(luout,
     &        'Orbital-dependent Mulliken population analysis for '//
     &        'the ionized channels',
     &        50, .true.)
         write(luout, '(" ")')

         do i = 1, nv
            write(luout, '(4x, "Channel : ", i4)') i
            call wfat_print_repeat(luout, 107, 6, '-')
            write(luout, '(x)')
            write(luout, '(6x, a4, a5, 7a14)') 
     &           "No.", "Atom",
     &           "Initial",
     &           "Nuclear",
     &           "2-electron",
     &           "DFT XC",
     &           "Effective",
     &           "v_s",
     &           "v_c"
            call wfat_print_repeat(luout, 107, 6, '-')
            write(luout, '(x)')
            
            do ia = 1, nat
               if (.not. geom_cent_tag(geom, ia, atomname))
     &              call errquit('WFAT: Getting atomname failed',
     &              52, GEOM_ERR)
               i_inpop = k_inpop - 1 + (i-1)*nat + ia
               i_nupop = k_nupop - 1 + (i-1)*nat + ia
               i_hfpop = k_hfpop - 1 + (i-1)*nat + ia
               i_xcpop = k_xcpop - 1 + (i-1)*nat + ia
               i_efpop = k_efpop - 1 + (i-1)*nat + ia
               i_vspop = k_vspop - 1 + (i-1)*nat + ia
               i_vcpop = k_vcpop - 1 + (i-1)*nat + ia
               write(luout, '(6x, i4, a5, 7f14.6)')
     &              ia, trim(atomname),
     &              dbl_mb(i_inpop),
     &              dbl_mb(i_nupop),
     &              dbl_mb(i_hfpop),
     &              dbl_mb(i_xcpop),
     &              dbl_mb(i_efpop),
     &              dbl_mb(i_vspop),
     &              dbl_mb(i_vcpop)
            enddo
            write(luout, '(" ")')
         enddo
      endif
      call ga_sync()


      if (.not. ma_free_heap(l_ov))
     &     call errquit('WFAT: Cannot deallocate overlap.', 1, MA_ERR)
      if (.not. ma_free_heap(l_inpop))
     &     call errquit('WFAT: Cannot deallocate vin_pop.', 1, MA_ERR)
      if (.not. ma_free_heap(l_pe))
     &     call errquit('WFAT: Cannot deallocate vnuc.', 1, MA_ERR)
      if (.not. ma_free_heap(l_nupop))
     &     call errquit('WFAT: Cannot deallocate vnucvin_pop.',
     &     1, MA_ERR)
      if (.not. ma_free_heap(l_ve))
     &     call errquit('WFAT: Cannot deallocate HF_vee.', 1, MA_ERR)
      if (.not. ma_free_heap(l_hfpop))
     &     call errquit('WFAT: Cannot deallocate vhfvin_pop.',
     &     1, MA_ERR)
      if (.not. ma_free_heap(l_vxc))
     &     call errquit('WFAT: Cannot deallocate XC_vee.', 1, MA_ERR)
      if (.not. ma_free_heap(l_xcpop))
     &     call errquit('WFAT: Cannot deallocate vxcvin_pop.',
     &     1, MA_ERR)
      if (.not. ma_free_heap(l_efpop))
     &     call errquit('WFAT: Cannot deallocate veffvin_pop.',
     &     1, MA_ERR)
      if (.not. ma_free_heap(l_hy))
     &     call errquit('WFAT: Cannot deallocate vhyd.', 1, MA_ERR)
      if (.not. ma_free_heap(l_vspop))
     &     call errquit('WFAT: Cannot deallocate vsvin_pop.',
     &     1, MA_ERR)
      if (.not. ma_free_heap(l_vcpop))
     &     call errquit('WFAT: Cannot deallocate vcfvin_pop.',
     &     1, MA_ERR)
      

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_expval_mulliken(nvec, bas, mat, vcf, geom,
     &                             pop)
      
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"      
      integer nvec                 ! [input]
      integer bas                  ! [input]
      double precision mat(1:*)    ! [input], dimension nb*nb
      double precision vcf(1:*)    ! [input], dimension nb*nvec
      integer geom                 ! [input]
      double precision pop(1:*)    ! [output], dimension nat*nvec

      integer i, icf1, icf2, ia, ib, ib1, ib2, icff, im1, im2, iia
      integer nat, nb
      double precision hh
      

      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 51, GEOM_ERR)
      if (.not. bas_numbf(bas, nb)) 
     &     call errquit('WFAT: Getting nb failed',
     &     53, BASIS_ERR)

cOLD      nov = nb * nb
cOLD      if (.not. ma_alloc_get(mt_dbl, nov, 'overlap', l_ov, k_ov))
cOLD     &     call errquit('WFAT: Cannot allocate overlap.', 1, MA_ERR)
cOLD      call int_1eov_all(bas, nb, dbl_mb(k_ov))


c     ==== Loop over input vectors ====
      do i = 1, nvec
         icf1 = (i-1)*nb + 1       ! icf1 = first element of the column of the bra vec
         icf2 = (i-1)*nb + nb    ! icf2 = last element of the column of the bra vec

c        ==== Loop over atoms ====
         do ia = 1, nat
            if (.not. bas_ce2bfr(bas, ia, ib1, ib2))
     &           call errquit('WFAT: Getting basis indices failed',
     &           60, BASIS_ERR)

c           ==== Loop over columns of overlap matrix ====
c           ====       corresponding to atom ia      ====
            hh = 0.0d0
            do ib = ib1, ib2
               icff = (i-1)*nb + ib      ! icff = column of ket vec
               im1 = (ib-1)*nb + 1       ! im1 = first element of the column of overlap matrix
               im2 = (ib-1)*nb + nb      ! im2 = last element of the column of overlap matrix
               hh = hh + vcf(icff) *
     &              dot_product(vcf(icf1:icf2), mat(im1:im2))
            enddo

            iia = (i-1)*nat + ia
            pop(iia) = hh
         enddo
      enddo


cOLD      if (.not. ma_free_heap(l_ov))
cOLD     &     call errquit('WFAT: Overlap deallocation failed.', 1, MA_ERR)


      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_matelm_mulliken(bas, mat, vcf1, vcf2, geom,
     &                                pop)
c     DESCRIPTION :
c     Partition the matrix element of an operator O (mat) between |psi> 
c     (vcf1) and |phi> (vcf2) in a similar way as Mulliken population,
c        <psi | O | phi> = \sum_{I=1}^{N_A}  N_A^{psi,phi}
c     where N_A^{psi,phi} is the contribution of atom A to the matrix
c     element and is returned in pop. The basis in which |psi> and |phi>
c     are represented can be either orthonormal or non-orthogonal but
c     normalized. If the overlap between |psi> and |phi> is sought, then
c     supply the overlap matrix for the argument mat.
      
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"      
      integer bas                  ! [input]
      double precision mat(1:*)    ! [input], dimension nb*nb
      double precision vcf1(1:*)   ! [input], dimension nb
      double precision vcf2(1:*)   ! [input], dimension nb
      integer geom                 ! [input]
      double precision pop(1:*)    ! [output], dimension nat

      integer ia, ib, ib1, ib2, icff, im1, im2
      integer nat, nb
      double precision hh
      

      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 51, GEOM_ERR)
      if (.not. bas_numbf(bas, nb)) 
     &     call errquit('WFAT: Getting nb failed',
     &     53, BASIS_ERR)


c     ==== Loop over atoms ====
      do ia = 1, nat
         if (.not. bas_ce2bfr(bas, ia, ib1, ib2))
     &        call errquit('WFAT: Getting basis indices failed',
     &        60, BASIS_ERR)

c        ==== Loop over columns of overlap matrix ====
c        ====       corresponding to atom ia      ====
         hh = 0.0d0
         do ib = ib1, ib2
            icff = ib      ! icff = column of ket vec
            im1 = (ib-1)*nb + 1       ! im1 = first element of the column of overlap matrix
            im2 = (ib-1)*nb + nb      ! im2 = last element of the column of overlap matrix
            hh = hh + vcf2(icff) *
     &           dot_product(vcf1(1:nb), mat(im1:im2))
         enddo

         pop(ia) = hh
      enddo


      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_ga_matelm_mulliken(geom, bas, g_v, g_w, g_a, pop)
c     d = < g_v | g_a | g_w >
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"      
      integer geom
      integer bas
      integer g_v
      integer g_w
      integer g_a
      double precision pop(1:*)

      integer typev, m_v, n_v, typew, m_w, n_w, typea, m_a, n_a
      integer nb, nat, ib1, ib2, ia, g_tmp


      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 51, GEOM_ERR)
      if (.not. bas_numbf(bas, nb)) 
     &     call errquit('WFAT: Getting nb failed',
     &     53, BASIS_ERR)      

      call ga_inquire(g_v, typev, m_v, n_v)
      call ga_inquire(g_w, typew, m_w, n_w)
      call ga_inquire(g_a, typea, m_a, n_a)


      if (typev == typew) then
         if (typew /= typea) then
            call errquit('wfat_ga_mattelm_mulliken: Type of argument '//
     &           '''g_w'' and of argument ''g_a'' must be the same.',
     &           1, UERR)
         endif
      else
         call errquit('wfat_ga_mattelm_mulliken: Type of argument '//
     &        '''g_v'' and of argument ''g_w'' must be the same.',
     &        1, UERR)
      endif
      

      if (.not. (m_v == 1 .or. n_v == 1)) then
         call errquit('wfat_ga_mattelm_mulliken: Argument ''g_v'' '//
     &        'must be a 1D array, i.e. a column or row vector.',
     &        1, UERR)
      endif
      if (m_v*n_v /= nb) then
         call errquit('wfat_ga_mattelm_mulliken: Number of elements '//
     &        'of g_v must be equal to the number of basis in bas.',
     &        1, UERR)
      endif
      if (m_w /= nb) then
         call errquit('wfat_ga_mattelm_mulliken: Number of elements '//
     &        'of g_w must be equal to the number of basis in bas.',
     &        1, UERR)
      endif
      if (n_w /= 1) then
         call errquit('wfat_ga_mattelm_mulliken: Argument gw must '//
     &        'be a column vector.', 1, UERR)
      endif
      if (m_a /= nb) then
         call errquit('wfat_ga_mattelm_mulliken: The size of g_a '//
     &        'must be nb by nb where nb is the no. of bases in bas.',
     &        1, UERR)
      endif
      if (m_a /= n_a) then
         call errquit('wfat_ga_mattelm_mulliken: The argument g_a '//
     &        'must be a nb by nb square array where nb is the no. '//
     &        'of bases in bas.', 1, UERR)
      endif
      


c      if (.not. ga_duplicate(g_v, g_tmp, 'tmp'))   ! g_tmp is a row vector.
c     &     call errquit('wfat_ga_mattelm_mulliken: Cannot create tmp',
c     &     1, GA_ERR)

cc     ==== g_v * g_a ====
c      call ga_dgemm('n', 'n', m_v, n_a, m_a, 1.0d0, g_v, g_a, 0.0d0,
c     &     g_tmp)

      
      if (.not. ga_duplicate(g_w, g_tmp, 'tmp'))   ! g_tmp is a column vector.
     &     call errquit('wfat_ga_mattelm_mulliken: Cannot create tmp',
     &     1, GA_ERR)
      

c     ==== g_a * g_w ====
      call ga_dgemm('n', 'n', m_a, n_w, m_w, 1.0d0, g_a, g_w, 0.0d0,
     &     g_tmp)
      
      
c     ==== g_v * (g_a * g_w) ====
      do ia = 1, nat
         if (.not. bas_ce2bfr(bas, ia, ib1, ib2))
     &        call errquit('WFAT: Getting basis indices failed',
     &        60, BASIS_ERR)
         if (m_v == 1) then
            pop(ia) = ga_ddot_patch(g_v,   'n', 1,   1,   ib1, ib2,
     &                              g_tmp, 'n', ib1, ib2, 1,   1)
         else if (n_v == 1) then
            pop(ia) = ga_ddot_patch(g_v,   'n', ib1, ib2, 1, 1,
     &                              g_tmp, 'n', ib1, ib2, 1, 1)
         endif
      enddo
      

      if (.not. ga_destroy(g_tmp))
     &     call errquit('wfat_ga_mattelm_mulliken: Cannot destroy tmp',
     &     1, GA_ERR)

      end
c=======================================================      

      
cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_expval_mulliken(nvec, bashan, vcf, geomhan,
cOLD     &                             pop)
cOLD      
cOLD      implicit none
cOLD#include "errquit.fh"
cOLD#include "mafdecls.fh"
cOLD#include "bas.fh"
cOLD#include "geom.fh"
cOLD#include "stdio.fh"      
cOLD      integer nvec                 ! [input]
cOLD      integer bashan               ! [input]
cOLD      double precision vcf(1:*)    ! [input], dimension nbas*nvec
cOLD      integer geomhan              ! [input]
cOLD      double precision pop(1:*)    ! [output], dimension natom*nvec
cOLD
cOLD      integer i, icf1, icf2, ia, ib, ib1, ib2, icff, iov1, iov2, iia
cOLD      integer natom, nbas, nov
cOLD      integer l_ov, k_ov
cOLD      double precision hh
cOLD      
cOLD
cOLD      if (.not. geom_ncent(geomhan, natom))
cOLD     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
cOLD      if (.not. bas_numbf(bashan, nbas)) 
cOLD     &     call errquit('int_1eov_all: Getting nbas0 failed',
cOLD     &     53, BASIS_ERR)
cOLD
cOLD      nov = nbas * nbas
cOLD      if (.not. ma_alloc_get(mt_dbl, nov, 'overlap', l_ov, k_ov))
cOLD     &     call errquit('WFAT: Cannot allocate overlap.', 1, MA_ERR)
cOLD      call int_1eov_all(bashan, nbas, dbl_mb(k_ov))
cOLD
cOLD
cOLDc      write(luout, '("Printing overlap matrix")')
cOLDc      call wfat_dprintmat(nbas, nbas, 'e', dbl_mb(k_ov))
cOLD
cOLD
cOLDc     ==== Loop over input vectors ====
cOLD      do i = 1, nvec
cOLD         icf1 = (i-1)*nbas + 1       ! icf1 = first element of the column of the bra vec
cOLD         icf2 = (i-1)*nbas + nbas    ! icf2 = last element of the column of the bra vec
cOLD
cOLDc        ==== Loop over atoms ====
cOLD         do ia = 1, natom
cOLD            if (.not. bas_ce2bfr(bashan, ia, ib1, ib2))
cOLD     &           call errquit('WFAT: Getting basis indices failed',
cOLD     &           60, BASIS_ERR)
cOLD
cOLDc           ==== Loop over columns of overlap matrix ====
cOLDc           ====       corresponding to atom ia      ====
cOLD            hh = 0.0d0
cOLD            do ib = ib1, ib2
cOLD               icff = (i-1)*nbas + ib                    ! icff = column of ket vec
cOLD               iov1 = k_ov - 1 + (ib-1)*nbas + 1         ! iov1 = first element of the column of overlap matrix
cOLD               iov2 = k_ov - 1 + (ib-1)*nbas + nbas      ! iov2 = last element of the column of overlap matrix
cOLD               hh = hh + vcf(icff) *
cOLD     &              dot_product(vcf(icf1:icf2), dbl_mb(iov1:iov2))
cOLD            enddo
cOLD
cOLD            iia = (i-1)*natom + ia
cOLD            pop(iia) = hh
cOLD         enddo
cOLD      enddo
cOLD
cOLD
cOLD      if (.not. ma_free_heap(l_ov))
cOLD     &     call errquit('WFAT: Overlap deallocation failed.', 1, MA_ERR)
cOLD
cOLD
cOLD      end
cOLDc=======================================================      
      
