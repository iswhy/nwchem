c===============================================================
c===============================================================
      subroutine mewfat_dipdiag(rtdb, geom, bas, nchar)

      implicit none
#include "inp.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "mewfat_qm_param.fh"
      integer rtdb
      integer geom
      integer bas
      integer nchar

      character*255 title, basname, scftype
      character(len=nchar) s
      character*2 si
      character*7 nm
      logical l, status
      integer l_oc, l_eg, l_mo1, l_mo2, l_mp, l_ovl, l_ds
      integer k_oc, k_eg, k_mo1, k_mo2, k_mp, k_ovl, k_ds
      integer i1, i2, ndets, nb, ns_, nmo, i, j, ii, is, nmp, max_ndets
      integer ns, nel(2), NMMT, g_ovl(2), i_mp, i_ds, i_ovl, i_mo1
      integer i_mo2, q, ndid
      double precision orig(3), dsl(3), osl, dthr
      parameter(max_ndets=10, NMMT=4)
      integer g_mo(2,max_ndets)      

      logical movecs_read_header, movecs_read
      integer ga_create_atom_blocked
      external movecs_read_header, movecs_read, ga_create_atom_blocked


      ndid = ga_nodeid()
      ns = 2

c      write(luout, '("n dipdiag = ", i4)') nchar

      if (.not. rtdb_cget(rtdb, 'mewfat:dipdiag', 1, s))
     &     call errquit('MEWFAT: Cannot read dipdiag', 32, RTDB_ERR)

c     ==== Get number of bases ====
      i1 = 0
      l = inp_strtok(s, ' ', i1, i2)
      status = movecs_read_header(s(i1:i2), title, basname, scftype,
     &     nb, ns_, nmo, 2)
      if (.not. ma_alloc_get(mt_dbl, nb, 'occ', l_oc, k_oc))
     &     call errquit('MEWFAT: Allocating occ failed', 14, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nb, 'erg', l_eg, k_eg))
     &     call errquit('MEWFAT: Allocating erg failed', 15, MA_ERR)


c     ==== Get the MO of the degenerate Slater determinants ====
      i = 0
      i1 = 0
      if (ndid == 0) write(luout, '("Reading movecs files :")')
      do while (inp_strtok(s, ' ', i1, i2))
         i = i + 1
         if (ndid == 0) write(luout, '(i3, 2x, a)') i, s(i1:i2)

         if (i > max_ndets) call errquit('MEWFAT: No. of degenerate '//
     &        'Slater determinants exceeds its maximum.', 78, UERR)

c        ==== Create global arrays for the MO ====
         write(si, '(i2.2)') i
         do is = 1, ns
            if (is == 1) then
               nm = 'mo-'//si//'-a'
               if (.not. ga_create(mt_dbl, nb, nb, nm, nb, 0,
     &              g_mo(is,i)))
     &              call errquit("MEWFAT: couldnt create "//nm,
     &              0, GA_ERR)
            else if (is == 2) then
               nm = 'mo_'//si//'-b'
               if (.not. ga_duplicate(g_mo(1,i), g_mo(2,i), nm))
     &              call errquit('MEWFAT: Couldnt create '//nm,
     &              0, GA_ERR)
            endif
            call ga_zero(g_mo(is,i))
         enddo

c        ==== Read the MO from movecs files ====
         do is = 1, ns
            if (ns_ == 2 .or. is == 1) then
               if (.not. movecs_read(s(i1:i2), is, dbl_mb(k_oc),
     $              dbl_mb(k_eg), g_mo(is,i)))
     $              call errquit("WFAT: Failed to read movecs data",
     &              0, UERR)
               if (i == 1) nel(is) = nint(sum(dbl_mb(k_oc:k_oc+nb-1)))      ! i = 1 only because all Slater determinants are assumed to have the same no. of electrons.
            else if (ns_ == 1 .and. is == 2) then
               call ga_copy(g_mo(1,i), g_mo(2,i))
               if (i == 1) then
                  if (mod(nel(1),2) == 0) then
                     nel(1) = nel(1) / 2
                  else
                     call errquit('WFAT: nel is odd while the two '//
     &                    'spin channels are identical.', 20, UERR)
                  endif
                  nel(2) = nel(1)
               endif
            endif
         enddo
                  
      enddo

      ndets = i

      
      if (.not. MA_alloc_get(MT_DBL, nb*nb*2, 'movec1', l_mo1, k_mo1))
     &     call errquit('WFAT: movec1 allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nb*nb*2, 'movec2', l_mo2, k_mo2))
     &     call errquit('WFAT: movec2 allocation failed', 84, MA_ERR)
      
      
c     ==== Dipole matrices in AO basis relative to origin ====
      nmp = nb * nb * NMMT
      if (.not. MA_alloc_get(MT_DBL, nmp, 'mpmat', l_mp, k_mp))
     &     call errquit('WFAT: mpmat allocation failed', 41, MA_ERR)
      orig = 0.0d0
      call int_mpole_all(bas, nb, 1, orig, dbl_mb(k_mp))

      
      do is = 1, 2
         g_ovl(is) = ga_create_atom_blocked(geom, bas, 'MO_overlap')
      enddo
      if (.not. ma_alloc_get(mt_dbl, nb*nb*2, 'moov', l_ovl, k_ovl))
     &     call errquit('WFAT: Cannot allocate moov', 1, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, ndets**2*3, 'dipsl', l_ds, k_ds))
     &     call errquit('WFAT: Cannot allocate dipsl', 1, MA_ERR)
      
      do j = 1, ndets         ! columns of dipole matrix
         do i = 1, ndets             ! rows of dipole matrix
c           ==== Overlap matrix between the two MO sets ====
            call mewfat_calc_moovl(geom, bas, oskel, g_mo(:,j), 
     &           g_mo(:,i), g_ovl)

c           ==== Move MOs and overlap matrix to local buffers ====
            do is = 1, 2
               i_ovl = k_ovl + (is-1)*nb*nb
               call ga_get(g_ovl(is), 1, nb, 1, nb, dbl_mb(i_ovl), nb)
               i_mo1 = k_mo1 + (is-1)*nb*nb
               call ga_get(g_mo(is,i), 1, nb, 1, nb, dbl_mb(i_mo1), nb)
               i_mo2 = k_mo2 + (is-1)*nb*nb
               call ga_get(g_mo(is,j), 1, nb, 1, nb, dbl_mb(i_mo2), nb)
            enddo

c           ==== Calculate total dipole matrix element ====
            i_mp = k_mp + (2-1)*nb**2
            call mewfat_dipslater(nb, nel, dbl_mb(k_mo1), dbl_mb(k_mo2),
     &           dbl_mb(k_ovl), dbl_mb(i_mp), dsl)
            call mewfat_ovlslater(nb, nel, dbl_mb(k_mo1), dbl_mb(k_mo2),
     &           dbl_mb(k_ovl), osl)
            if (ndid == 0) then
               write(luout, '(2x, "Overlap between det. ", i2, '//
     &              '" and ", i2, " = ", es14.6)') i, j, osl
               write(luout, '(2x, "Dipole between det. ", i2, '//
     &              '" and ", i2, " (x y z) = ", 3es14.6)') i, j, dsl
            endif

            do q = 1, 3
               i_ds = k_ds + (q-1)*ndets**2 + (j-1)*ndets + (i-1)
               dbl_mb(i_ds) = dsl(q)

c              == Hermitian: (j,i)=(i,j) ==
c               i_ds = k_ds + (q-1)*ndets**2 + (i-1)*ndets + (j-1)
c               dbl_mb(i_ds) = dsl(q)
            enddo
         enddo
      enddo


c     ==== Print the diagonalizing Slater determinants ====
      if (.not. rtdb_get(rtdb, 'wfat:dthr', mt_dbl, 1, dthr))
     &     call errquit('WFAT: Failed reading dthr from rtdb',
     &     1, RTDB_ERR)
      call mewfat_print_dipdiag(rtdb, ndets, dbl_mb(k_ds), dthr)
      

c     ==== Deallocations ====
      do is = 1, 2
         if (.not. ga_destroy(g_ovl(is)))
     &        call errquit('WFAT: Cannot destroy MO_overlap.',
     &        1, GA_ERR)
      enddo
      if (.not. ma_free_heap(l_oc))
     &     call errquit('WFAT: Cannot deallocate occ', 64, MA_ERR)
      if (.not. ma_free_heap(l_eg))
     &     call errquit('WFAT: Cannot deallocate erg', 64, MA_ERR)
      if (.not. ma_free_heap(l_mo1))
     &     call errquit('WFAT: Cannot deallocate movec1', 64, MA_ERR)
      if (.not. ma_free_heap(l_mo2))
     &     call errquit('WFAT: Cannot deallocate movec2', 64, MA_ERR)
      if (.not. ma_free_heap(l_mp))
     &     call errquit('WFAT: Cannot deallocate mpmat', 64, MA_ERR)
      if (.not. ma_free_heap(l_ovl))
     &     call errquit('WFAT: Cannot deallocate moov', 64, MA_ERR)
      if (.not. ma_free_heap(l_ds))
     &     call errquit('WFAT: Cannot deallocate dipsl', 64, MA_ERR)

      do i = 1, ndets
         do is = 1, ns
            if (.not. ga_destroy(g_mo(is,i)))
     &        call errquit('WFAT: Cannot destroy MO.', 1, GA_ERR)
         enddo
      enddo
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine mewfat_print_dipdiag(rtdb, ndets, dsl, dthr)
c     DESCRIPTION:
c     Prints the matrix elements of the unitary transform that diagonalizes
c     the lab frame z-component dipole matrix at each orientation angle into
c     a file.
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "wfat.fh"
      integer rtdb
      integer ndets
      double precision dsl(ndets,ndets,3)
      double precision dthr

      character*255 fname
      logical isdiag
      integer ig, ib, k, i, j, i_it, ngm, nbt, l_it, k_it, l_ev, k_ev
      double precision bt_start, bt_last, dbt, gm_start, gm_last, dgm
      double precision bt, gm, protmat(3,3), cu(3)

      
c     ==== Obtaining info on beta ====
      call wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)


c     ==== Obtaining info on gamma ====
      call wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)

      if (.not. MA_alloc_get(MT_DBL, ndets**2, 'inth', l_it, k_it))
     &     call errquit('WFAT: inth allocation failed', 2, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, ndets, 'eigval', l_ev, k_ev))
     &     call errquit('WFAT: eigval allocation failed', 2, MA_ERR)

c     ==== File name of the unitary matrix file ====
      call util_file_name('wfatd', .false., .false., fname)
      
      open(WF_MEDIPDIAG_UNIT, file=fname, action='write')
      do ig = 1, ngm            ! Loop over gamma
         gm = gm_start + (ig-1)*dgm
         do ib = 1, nbt            ! Loop over beta
            bt = bt_start + (ib-1)*dbt
            
            write(WF_MEDIPDIAG_UNIT, '(2es14.6, 3x)', advance='no')
     &           bt, gm

c           ==== Construct the lab frame z-component dipole matrix ====
            call wfat_passrot(bt, gm, protmat)
            cu = protmat(:,3)
            dbl_mb(k_it:k_it+ndets**2-1) = 0.0d0
            do k = 1, 3
               call wfat_matrix_sumabsorb(ndets, ndets, dbl_mb(k_it),
     &              dsl(:,:,k), 1.0d0, cu(k))
            enddo

c           ==== Diagonalization ====
            call wfat_isdiag(ndets, dbl_mb(k_it), dthr, isdiag)
            if (isdiag) then
               call wfat_get_identity(ndets, dbl_mb(k_it))
            else
               call wfat_dsyev_wrap(ndets, dbl_mb(k_it), dbl_mb(k_ev))
            endif
            do j = 1, ndets      ! Loop over diagonalizing states
               do i = 1, ndets      ! Loop over non-diagonalizing determinants
                  i_it = k_it + (j-1)*ndets + (i-1)
                  write(WF_MEDIPDIAG_UNIT, '(es17.6)', advance='no')
     &                 dbl_mb(i_it)
               enddo
               write(WF_MEDIPDIAG_UNIT, '("     ")', advance='no')
            enddo
            write(WF_MEDIPDIAG_UNIT, '("")')
            
         enddo
         write(WF_MEDIPDIAG_UNIT, '("")')
      enddo


      if (.not. MA_free_heap(l_it))
     &     call errquit('WFAT: inth deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_ev))
     &     call errquit('WFAT: eigval deallocation failed', 3, MA_ERR)

      close(WF_MEDIPDIAG_UNIT)

      end
c===============================================================
      
