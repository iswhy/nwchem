c=================================================
c=================================================
      subroutine mewfat_calc_detcoef(geom, bas, nb, sx, sy, n0, n1,
     &           cx, cy, g_mo0, g_mo1, oskel, ocmo0, ocmo1,
     &           pcf, rcf, qcf, scf)

!     g_pcf = det( S_x[0|(xj)] )          = P(j)
!     g_rcf = det( S_y[0|(xj)] )          = R
!     g_qcf = det( S_x[(xk')|(xj),(xk)] ) = Q(j,k,k')
!     g_scf = det( S_y[(yk')|(xj),(yk)] ) = S(k,k')
!     If sx=1 then sy=2. Likewise, if sx=2 then sy=1.
      
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer geom
      integer bas
      integer nb
      integer sx                                  ! [input], if sx=1, then sy must be 2, and vice versa.
      integer sy                                  ! [input], if sy=1, then sx must be 2, and vice versa.
      integer n0(2)                               ! [input], no. of electrons in the neutral.
      integer n1(2)                               ! [input], no. of electrons in the cation.
      integer cx                                  ! [input], must be equal to max(1,n0(sx)-1).
      integer cy                                  ! [input], must be equal to max(1,n0(sy)).
      integer g_mo0(2)                            ! [input], neutral MO's.
      integer g_mo1(2)                            ! [input], cation MO's.
      logical oskel                               ! [input]   
      double precision ocmo0(nb, n0(1)+n0(2))     ! [output], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))     ! [output], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                ! [output], its elements index neutral MO's.
      double precision rcf                        ! [output], a scalar.
      double precision qcf(n0(sx), n0(sx), cx)    ! [output], 1st and 2nd dimensions index neutral MO's while the 3rd indexes cation MO's.
      double precision scf(cy, cy)                ! [output], rows of scf index neutral MO's while its columns index cation MO's.

      integer ndid, is, g_ovl(2), i_ovl, nn, nsc, i1, i2
      integer l_ovl, l_iv, l_sc, k_ovl, k_iv, k_sc

      integer ga_create_atom_blocked
      external ga_create_atom_blocked


      ndid = ga_nodeid()
      if (cx /= max(1,n0(sx)-1)) then
         if (ndid == 0) write(luout, '(3x, "cx = ", i4)') cx
         if (ndid == 0) write(luout, '(3x, "max(1,n0(sx)-1) = ", i4)')
     &        max(1,n0(sx)-1)
         call errquit('mewfat_calc_detcoef: The argument cx must be '//
     &        'equal to max(1,n0(sx)-1).', 1, UERR)
      endif

      if (cy /= max(1,n0(sy))) then
         if (ndid == 0) write(luout, '(3x, "cy = ", i4)') cy
         if (ndid == 0) write(luout, '(3x, "max(1,n0(sy)) = ", i4)')
     &        max(1,n0(sy))
         call errquit('mewfat_calc_detcoef: The argument cy must be '//
     &        'equal to max(1,n0(sy)).', 1, UERR)
      endif
      

c     ==== Check no. of electrons ====
      call mewfat_check_nelec(sx, sy, n0, n1, 100)
      

c     ==== Calculate overlap between MO's ====
      do is = 1, 2
         g_ovl(is) = ga_create_atom_blocked(geom, bas, 'MO_overlap')
         call ga_zero(g_ovl(is))
      enddo
      call mewfat_calc_moovl(geom, bas, oskel, g_mo0, g_mo1, g_ovl)
      if (.not. ma_alloc_get(mt_dbl, nb*nb*2, 'moov', l_ovl, k_ovl))
     &     call errquit('WFAT: Cannot allocate moov', 1, MA_ERR)
      do is = 1, 2
         i_ovl = k_ovl + (is-1)*nb*nb
         call ga_get(g_ovl(is), 1, nb, 1, nb, dbl_mb(i_ovl), nb)
      enddo


c     ==== Obtain occupied MO's of neutral and cation ====
      do is = 1, 2
         if (n0(is) > 0) then
            i1 = (is-1)*n0(1) + 1
            i2 = i1 + n0(is) - 1
            call ga_get(g_mo0(is), 1, nb, 1, n0(is), ocmo0(:,i1:i2), nb)
         endif

         if (n1(is) > 0) then
            i1 = (is-1)*n1(1) + 1
            i2 = i1 + n1(is) - 1
            call ga_get(g_mo1(is), 1, nb, 1, n1(is), ocmo1(:,i1:i2), nb)
         endif
      enddo


c     ==== Calculate determinant coefficients ====
      nn = maxval(n0)
      if (.not. ma_alloc_get(mt_int, nn, 'ipvt', l_iv, k_iv)) 
     &     call errquit('WFAT: failed to allocate ipvt', 0, MA_ERR)
      nsc = nn*(nn+1)         ! The '+1' is for the 'work' space.
      if (.not. ma_alloc_get(mt_dbl, nsc, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: Cannot allocate scratch.', 1, MA_ERR)

c     ==== P det. coefficients ====
      if (n0(sx) == 1) then
         pcf = 1.0d0
      else
         call mewfat_calc_pcf(sx, n0, nb, dbl_mb(k_ovl), dbl_mb(k_sc),
     &        int_mb(k_iv), pcf)
      endif

c     ==== R det. coefficients ====
      if (n0(sy) == 0) then
         rcf = 1.0d0
      else
         call mewfat_calc_rcf(sy, n0, nb, dbl_mb(k_ovl), dbl_mb(k_sc),
     &        int_mb(k_iv), rcf)
      endif

c     ==== Q det. coefficients ====
      if (n0(sx) == 1) then
         qcf = 0.0d0
      else if (n0(sx) == 2) then
         call mewfat_calc_qcf_n2(qcf)
      else
         call mewfat_calc_qcf(sx, n0, nb, dbl_mb(k_ovl), dbl_mb(k_sc),
     &        int_mb(k_iv), qcf)
      endif

c     ==== S det. coefficients ====
      if (n0(sy) == 0) then
         scf = 0.0d0
      else if (n0(sy) == 1) then
         scf = 1.0d0
      else
         call mewfat_calc_scf(sy, n0, nb, dbl_mb(k_ovl), dbl_mb(k_sc),
     &        int_mb(k_iv), scf)
      endif


c     ==== Deallocations ====
      do is = 1, 2
         if (.not. ga_destroy(g_ovl(is)))
     &        call errquit('WFAT: Cannot destroy MO_overlap.',
     &        1, GA_ERR)
      enddo
      if (.not. ma_free_heap(l_ovl))
     &     call errquit('WFAT: Cannot deallocate moov', 64, MA_ERR)
      if (.not. ma_free_heap(l_iv))
     &     call errquit('WFAT: Cannot deallocate ipvt', 64, MA_ERR)
      if (.not. ma_free_heap(l_sc))
     &     call errquit('WFAT: Cannot deallocate scratch', 64, MA_ERR)
            
      end
c=================================================

      
c=================================================
c=================================================
      subroutine mewfat_calc_moovl(geom, bas, oskel, g_mo1, g_mo2,
     &                             g_ovl)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "bas.fh"      
      integer geom              ! [input]
      integer bas               ! [input]
      logical oskel             ! [input]
      integer g_mo1(2)          ! [input], column orbitals.
      integer g_mo2(2)          ! [input], row orbitals.
      integer g_ovl(2)          ! [output], its rows index g_mo2 vectors while its columns index g_mo1 vectors.

      integer nb, g_s, g_tmp, is

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      

      if (.not. bas_numbf(bas, nb))
     &     call errquit('WFAT: Getting nbas0 failed in ', 53, UERR)
      

      g_s = ga_create_atom_blocked(geom, bas, 'bas_overlap')
      call ga_zero(g_s)
      call int_1e_ga(bas, bas, g_s, 'overlap', oskel)

      if (.not. ga_duplicate(g_mo1(1), g_tmp, 'tmp'))
     &     call errquit('WFAT: Cannot create tmp in mewfat_calc_moovl.',
     &     1, GA_ERR)
      

      do is = 1, 2
         call ga_zero(g_ovl(is))
         call ga_zero(g_tmp)
         call ga_dgemm('n', 'n', nb, nb, nb, 1.0d0, g_s, g_mo1(is), 
     &        0.0d0, g_tmp)
         call ga_dgemm('t', 'n', nb, nb, nb, 1.0d0, g_mo2(is), g_tmp, 
     &        0.0d0, g_ovl(is))
      enddo

      
      if (.not. ga_destroy(g_tmp))
     &     call errquit('WFAT: Cannot destroy tmp in '//
     &     'mewfat_calc_moovl.', 1, GA_ERR)
      if (.not. ga_destroy(g_s))
     &     call errquit('WFAT: Cannot destroy bas_overlap in '//
     &     'mewfat_calc_moovl.', 1, GA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_calc_pcf(sx, nel, nb, ovl, scr, ipvt, pcf)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"      
      integer sx                            ! [input], sx = 1 or 2, spin of the ionized electron.
      integer nel(2)                        ! [input], no. of electrons in the neutral.
      integer nb                            ! [input], no. of basis.
      double precision ovl(nb, nb, 2)       ! [input], The rows of ovl must correspond to cation orbitals. the columns must correspond to neutral orbitals.
      double precision scr(nel(sx)-1, 1:*)  ! [in/output], Minimum no. of columns should be nel(sx).
      integer ipvt(1:*)                     ! [in/output], Minimum no. of elements should be nel(sx)-1.
      double precision pcf(nel(sx))         ! [output], its elements index neutral MO's.

      integer i, nel_, iwork, info
      double precision det(2)

      
      if (sx /= 1 .and. sx /= 2)
     &     call errquit('WFAT: Unknown value for the argument ''sx'' '//
     &     'of mewfat_calc_pcf.', 1, UERR)
      if (nel(sx) == 1)
     &     call errquit('mewfat_calc_pcf: Argument nel(sx) cannot be '//
     &     'equal to unity. When nel(sx)=1, pcf only has one element '//
     &     'and is equal to unity.', 1, UERR)
      

      iwork = nel(sx)
      nel_ = nel(sx) - 1            ! nel_ x nel_ is the size of the modified overlap matrix

      
c     ==== Build overlap matrices where ====
c     ====    col i has been removed    ====
      do i = 1, nel(sx)              ! Loop over elements of pcf, which indexes neutral MO's.
c        ==== Left part of modified overlap matrix ====
         if (i > 1) then
            scr(1:nel_, 1:i-1) = ovl(1:nel_, 1:i-1, sx)
         endif
         
c        ==== Right part of modified overlap matrix ====
         if (i < nel(sx)) then
            scr(1:nel_, i:nel_) = ovl(1:nel_, i+1:nel(sx), sx)
         endif

c        ==== Calculate the determinant ====
         call dgefa(scr(:, 1:nel_), nel_, nel_, ipvt, info)
         if (info /= 0) call errquit('WFAT: Error in dgefa', 0, UERR)
         call dgedi(scr(:, 1:nel_), nel_, nel_, ipvt, det, 
     &              scr(:, iwork), 10)
         pcf(i) =  det(1) * 10.0d0**det(2)
      enddo

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_calc_rcf(sy, nel, nb, ovl, scr, ipvt, rcf)

      implicit none
#include "errquit.fh"
      integer sy                         ! [input], sy = 1 or 2, spin of the ionized electron.
      integer nel(2)                     ! [input], no. of electrons in the neutral.
      integer nb                         ! [input], no. of basis.
      double precision ovl(nb, nb, 2)    ! [input], The rows of ovl must correspond to cation orbitals. the columns must correspond to neutral orbitals.
      double precision scr(nel(sy), 1:*) ! [in/output], Minimum no. of columns should be nel(sy)+1.
      integer ipvt(1:*)                  ! [in/output], Minimum no. of elements should be nel(sy).
      double precision rcf               ! [output], The desired coefficient.

      integer iwork, nel_, info
      double precision det(2)


      if (nel(sy) == 0)
     &     call errquit('mewfat_calc_rcf: Argument nel(sy) cannot be '//
     &     'equal to zero. When nel(sy)=0, rcf is equal to unity.',
     &     1, UERR)
      

      iwork = nel(sy) + 1
      nel_ = nel(sy)
      
      scr(1:nel_, 1:nel_) = ovl(1:nel_, 1:nel_, sy)

      
      call dgefa(scr(:, 1:nel_), nel_, nel_, ipvt, info)
      if (info /= 0) call errquit('WFAT: Error in dgefa', 0, UERR)
      call dgedi(scr(:, 1:nel_), nel_, nel_, ipvt, det,
     &           scr(:, iwork), 10)
      rcf =  det(1) * 10.0d0**det(2)
      
      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_calc_qcf(sx, nel, nb, ovl, scr, ipvt, qcf)

      implicit none
#include "errquit.fh"
      integer sx                                          ! [input], sx = 1 or 2, spin of the ionized electron.
      integer nel(2)                                      ! [input], no. of electrons in the neutral.
      integer nb                                          ! [input], no. of basis.
      double precision ovl(nb, nb, 2)                     ! [input], The rows of ovl must correspond to cation orbitals. the columns must correspond to neutral orbitals.
      double precision scr(nel(sx)-2, 1:*)                ! [in/output], Minimum no. of columns should be nel(sx)-1.
      integer ipvt(1:*)                                   ! [in/output], Minimum no. of elements should be nel(sx)-2.
      double precision qcf(nel(sx), nel(sx), nel(sx)-1)   ! [output], 1st and 2nd dimensions index neutral MO's while the 3rd indexes cation MO's.

      integer iwork, nel_, kk, j, k, jj, os, info
      double precision det(2)


      if (nel(sx) <= 2)
     &     call errquit('mewfat_calc_qcf: Argument nel(sx) cannot be '//
     &     'smaller than 3. When nel(sx)=1, qcf only has one element '//
     &     'and is equal to zero. When nel(sx)=2, call '//
     &     '''mewfat_calc_qcf_n2'' instead.', 1, UERR)
      

      iwork = nel(sx) - 1
      nel_ = nel(sx) - 2
      

      qcf = 0.0d0
      do kk = 1, nel(sx)-1     ! Loop over the 3rd dimension of qcf, which indexes the cation MO's.
         do j = 1, nel(sx)     ! Loop over the 2nd dimension of qcf, which indexes the neutral MO's.
            do k = 1, j-1      ! k < j. Loop over the 1st dimension of qcf, which indexes the neutral MO's.
               
c              ==== Build overlap matrix where row kk, ====
c              ==== col j, and col k have been removed ====
               do jj = 1, nel_
                  if (jj < k) os = 0
                  if (jj >= k .and. jj < j-1) os = 1
                  if (jj >= j-1) os = 2

c                 ==== Upper part of col. jj ====
                  if (kk > 1) then
                     scr(1:kk-1, jj) = ovl(1:kk-1, jj+os, sx)
                  endif
                  
c                 ==== Lower part of col. jj ====
                  if (kk < nel(sx)-1) then
                     scr(kk:nel_, jj) = ovl(kk+1:nel(sx)-1, jj+os, sx)
                  endif                  
               enddo

c              ==== Calculate the determinant ====
               call dgefa(scr(:, 1:nel_), nel_, nel_, ipvt, info)
               if (info /= 0) call errquit('WFAT: Error in dgefa',
     &                                     0, UERR)
               call dgedi(scr(:, 1:nel_), nel_, nel_, ipvt, det,
     &                    scr(:, iwork), 10)               
               qcf(k,j,kk) = det(1) * 10.0d0**det(2)
               
               qcf(j,k,kk) = qcf(k,j,kk)
            enddo
         enddo
      enddo

      end
c=================================================

      
c=================================================
c=================================================
      subroutine mewfat_calc_qcf_n2(qcf)

      implicit none
      double precision qcf(2,2)

      qcf = 0.0d0
      qcf(1,2) = 1.0d0
      qcf(2,1) = 1.0d0

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_calc_scf(sy, nel, nb, ovl, scr, ipvt, scf)

      implicit none
#include "errquit.fh"
      integer sy                               ! [input], complementary spin channel
      integer nel(2)                           ! [input], no. of electrons in the neutral.
      integer nb                               ! [input], no. of bases.
      double precision ovl(nb, nb, 2)          ! [input], rows of ovl index cation MO's while its columns index neutral MO's.
      double precision scr(nel(sy)-1, 1:*)     ! [in/output], Minimum no. of columns should be nel(sy).
      integer ipvt(1:*)                        ! [in/output], Minimum no. of elements should be nel(sy)-1.
      double precision scf(nel(sy), nel(sy))   ! [output], rows of scf index neutral MO's while its columns index cation MO's.

      integer iwork, nel_, i, j, jj, os, info
      double precision det(2)


      if (nel(sy) <= 1)
     &     call errquit('mewfat_calc_scf: Argument nel(sy) cannot be '//
     &     'smaller than 2. When nel(sy)=0, scf only has one element '//
     &     'and is equal to zero. When nel(sy)=1, it also has one  '//
     &     'element and is equal to unity.', 1, UERR)
      
      
      iwork = nel(sy)
      nel_ = nel(sy) - 1          ! nel_ x nel_ is the size of the modified overlap matrix
            
      do i = 1, nel(sy)           ! Loop over columns of scf, which indexes cation MO's. 
         do j = 1, nel(sy)        ! Loop over rows of scf, which indexes neutral MO's.

c           ==== Build overlap matrix where row i ====
c           ====    and col j have been removed   ====
            do jj = 1, nel_     ! Loop over cols of modified overlap matrix
               if (jj < j) os = 0
               if (jj >= j) os = 1

c              ==== Upper part of col. jj ====
               if (i > 1) then
                  scr(1:i-1, jj) = ovl(1:i-1, jj+os, sy)
               endif
               
c              ==== Lower part of col. jj ====
               if (i < nel(sy)) then
                  scr(i:nel_, jj) = ovl(i+1:nel(sy), jj+os, sy)
               endif
            enddo

c           ==== Calculate the determinant ====
            call dgefa(scr(:, 1:nel_), nel_, nel_, ipvt, info)
            if (info /= 0) call errquit('WFAT: Error in dgefa', 0, UERR)
            call dgedi(scr(:, 1:nel_), nel_, nel_, ipvt, det,
     &                 scr(:, iwork), 10)
            scf(j,i) =  det(1) * 10.0d0**det(2)
         enddo
      enddo
            
      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_copy_2dblock(ma, na, ma_, na_, a, ord,
     &                             mb, nb, mb_, nb_, b)

      implicit none
#include "errquit.fh"      
      integer ma
      integer na
      integer ma_(2)
      integer na_(2)
      double precision a(ma, na)
      integer ord(2)
      integer mb
      integer nb
      integer mb_(2)
      integer nb_(2)
      double precision b(mb, nb)

      integer mma, nna, mmb, nnb
      

      mma = ma_(2) - ma_(1) + 1    ! row
      nna = na_(2) - na_(1) + 1    ! col
      if (mma > ma) then
         call errquit('WFAT: Row size of the block of source matrix '//
     &        'exceeds the row size of the latter', 1, UERR)
      else if (nna > na) then
         call errquit('WFAT: Column size of the block of source '//
     &        'matrix exceeds the column size of the latter', 1, UERR)
      endif


      mmb = mb_(2) - mb_(1) + 1    ! row
      nnb = nb_(2) - nb_(1) + 1    ! col
      if (mmb > mb) then
         call errquit('WFAT: Row size of the block of target matrix '//
     &        'exceeds the row size of the latter', 1, UERR)
      else if (nnb > nb) then
         call errquit('WFAT: Column size of the block of target '//
     &        'matrix exceeds the column size of the latter', 1, UERR)
      endif


      if (mma*nna /= mmb*nnb)
     &     call errquit('WFAT: No. of elements of the block of '//
     &     'source matrix must be the same as that of the target '//
     &     'matrix', 1, UERR)
      
      
      b(mb_(1):mb_(2), nb_(1):nb_(2)) =
     &     reshape(a(ma_(1):ma_(2), na_(1):na_(2)), (/ mmb, nnb /),
     &     ORDER=ord)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_check_spin(sx, sy, i)

c     DESCRIPTION:
C        Only pass this check when sx=1 and sy=2, or when sx=2 and sy=1.
      implicit none
#include "errquit.fh"
      integer sx, sy, i


      if (sx == sy) 
     &     call errquit('WFAT: Argument sx and sy '//
     &     'cannot be equal', i, UERR)
      if (sx /= 1 .and. sx /= 2)
     &     call errquit('WFAT: Argument sx '//
     &     'must be equal to either 1 or 2.', i, UERR)
      if (sy /= 1 .and. sy /= 2)
     &     call errquit('WFAT: Argument sy '//
     &     'must be equal to either 1 or 2.', i, UERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_check_nelec(sx, sy, n0, n1, i)

      implicit none
#include "errquit.fh"
      integer sx, sy, n0(2), n1(2), i


      call mewfat_check_spin(sx, sy, 100)

      if (n1(sx) /= n0(sx)-1)
     &     call errquit('WFAT: Inconsistency in the no. of '//
     &     'electrons in the ionized channel of neutral '//
     &     'and cation', i, UERR)
      if (n1(sy) /= n0(sy))
     &     call errquit('WFAT: Inconsistency in the no. of '//
     &     'electrons in the complementary channel of neutral '//
     &     'and cation', i, UERR)
      
      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_print_detcoef(sx, n0, cx, cy, pcf, rcf, qcf,
     &     scf)

      implicit none
#include "stdio.fh"
      integer sx
      integer n0(2)
      integer cx                                  ! [input], must be equal to max(1,n0(sx)-1).
      integer cy                                  ! [input], must be equal to max(1,n0(sy)).
      double precision pcf(n0(sx))                ! [input], its elements index neutral MO's.
      double precision rcf                        ! [input], a scalar.
      double precision qcf(n0(sx), n0(sx), cx)    ! [output], 1st and 2nd dimensions index neutral MO's while the 3rd indexes cation MO's.
      double precision scf(cy, cy)                ! [input], rows of scf index neutral MO's while its columns index nation MO's.

      integer i, j, g, j_in, j_fin, ng
      double precision qmean, qstd


      write(luout, '(2/)')
      call util_print_centered(luout, 'Determinant coefficients', 40,
     &     .true.)
      write(luout, '(" ")')

      
      write(luout, '("   *** P determinant coefficients ***")')
      write(luout, '("     ")', advance='no')
      do i = 1, n0(sx)
         write(luout, '(i3, ":", f12.8, "  ")', advance='no') i, pcf(i)
         if (mod(i,5) == 0) then
            write(luout, '("")')
            write(luout, '("     ")', advance='no')
         endif
         if (i == n0(sx)) write(luout, '("")')
      enddo
      write(luout, '("")')


      write(luout, '("   *** R determinant coefficients ***")')
      write(luout, '("     ")', advance='no')
      write(luout, '(i3, ": ", f13.8, " ")') 1, rcf
      write(luout, '("")')


      write(luout, '("   *** S determinant coefficients ***")')
      ng = int(cy/5) + 1
      j_in = 1
      j_fin = min(j_in + 5 - 1, cy)
      do g = 1, ng     !Loop over groups

c        ==== Columns title ====
         write(luout, '(5x, "        ")', advance='no')
         do j = j_in, j_fin
            write(luout, '(5x, " col", i3, ":")', advance='no') j
         enddo
         write(luout, '(" ")')

c        ==== Determinant values ====
         do i = 1, cy     ! Loop over rows
            write(luout, '(5x, " row", i3, ":")', advance='no') i
            do j = j_in, j_fin     ! Loop over columns
               write(luout, '(f13.8)', advance='no') scf(i,j)
            enddo
            write(luout, '("")')
         enddo
         write(luout, '("")')

         j_in = j_in + 5
         j_fin = min(j_fin + 5, cy)
      enddo


      write(luout, '("   *** Q determinant coefficients ***")')
      write(luout, '("     The Q determinant coefficients is a 3D '//
     &     'tensor with many elements. Instead of printing all of '//
     &     'them, a statistics below is provided.")')
      write(luout, '(5x, "Max. value = ", f13.8)') maxval(qcf)
      write(luout, '(5x, "Min. value = ", f13.8)') minval(qcf)
      qmean = sum(qcf) / (n0(sx)*n0(sx)*cx)
      write(luout, '(5x, "Mean value = ", f13.8)') qmean
      qstd = sqrt( sum((qcf-qmean)**2) / (n0(sx)*n0(sx)*cx) )
      write(luout, '(5x, "Std. deviation = ", f13.8)') qstd
      write(luout, '("")')

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_dipslater(nb, nel, mo1, mo2, ovl, dip, dsl)
c     DESCRIPTION:
c     Calculates <Sl(MO1)|D|Sl(MO2)>
c     where Sl(MO) is a Slater determinant consisting of molecular orbitals MO and D is
c     one of x, y, and z position operators.
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nb
      integer nel(2)
      double precision mo1(nb,nb,2)   ! [input], AO coefficients of MO1.
      double precision mo2(nb,nb,2)   ! [input], AO coefficients of MO2.
      double precision ovl(nb,nb,2)   ! [input], overlap matrix between MO1 (row) and MO2 (col).
      double precision dip(nb,nb,3)   ! [input], dipole matrix in AO basis.
      double precision dsl(3)         ! [output], dipole matrix elements between the two determinants represented by MO1 and MO2.

      integer l_tm, l_d, l_iv, l_mi, l_ovs
      integer k_tm, k_d, k_iv, k_mi, k_ovs
      integer is, js, nel_, ntmp, i_d, q, niv, i1, i2, info, g1, g2
      double precision det0(2), det, h(3)

      integer delkron
      external delkron
      

      dsl = 0.0d0
      do is = 1, 2

c        ==== Temporary array ====
         nel_ = nel(is) - 1
         ntmp = max(nb*nel(is), nel_, nel(js))
         if (.not. ma_alloc_get(mt_dbl, ntmp, 'tmp', l_tm, k_tm))
     &        call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
         
c        ==== Dipole matrix between MO1 (row) and MO2 (col) ====
         if (.not. ma_alloc_get(mt_dbl, 3*nel(is)**2, 'dipmo',
     &        l_d, k_d))
     &        call errquit('WFAT: Cannot allocate dipmo', 1, MA_ERR)
         do q = 1, 3
            i_d = k_d + (q-1)*nel(is)**2
            call wfat_dsymm_wrap(nb, nel(is), 1.0d0, dip(:,:,q), 
     &           mo2(:,1:nel(is),is), 0.0d0, dbl_mb(k_tm))
            call wfat_dgemm_wrap('t', 'n', nb, nel(is), nb, nel(is), 
     &           1.0d0, mo1(:,1:nel(is),is), dbl_mb(k_tm), 0.0d0,
     &           dbl_mb(i_d))
         enddo
                           
         if (is == 1) js = 2
         if (is == 2) js = 1
         niv = max(nel_, nel(js))
         if (.not. ma_alloc_get(mt_int, niv, 'ipvt', l_iv, k_iv))
     &        call errquit('WFAT: Cannot allocate ipvt', 18, MA_ERR)
         if (.not. ma_alloc_get(mt_dbl, nel_**2, 'minor', l_mi, k_mi))
     &        call errquit('WFAT: Cannot allocate minor', 18, MA_ERR)
         h = 0.0d0
         do i2 = 1, nel(is)
            do i1 = 1, nel(is)
c              ==== Minor matrix for the current (i1,i2) ====
               call wfat_minor(i1, i2, nel(is), nel(is), 
     &              ovl(1:nel(is),1:nel(is),is), dbl_mb(k_mi))

c              ==== Calculate the determinant of ====
c              ==== the current spin block minor ====
               call dgefa(dbl_mb(k_mi), nel_, nel_, int_mb(k_iv), info)
               if (info /= 0)
     &              call errquit('WFAT: Error in dgefa', 45, UERR)
               call dgedi(dbl_mb(k_mi), nel_, nel_, int_mb(k_iv), det0, 
     &                    dbl_mb(k_tm), 10)
               det = det0(1) * 10.0d0**det0(2)

c              ==== Accumulate ====
               g1 = (-1)**(delkron(2,is)*nel(1)+i1-1)
               g2 = (-1)**(delkron(2,is)*nel(1)+i2-1)
               do q = 1, 3
                  i_d = k_d + (q-1)*nel(is)**2 + (i2-1)*nel(is) + (i1-1)
                  h(q) = h(q) + dbl_mb(i_d) * det * g1 * g2
               enddo
               
            enddo
         enddo

         if (.not. MA_free_heap(l_d))
     &        call errquit('WFAT: dipmo deallocation failed', 7, MA_ERR)
         if (.not. MA_free_heap(l_mi))
     &        call errquit('WFAT: minor deallocation failed', 7, MA_ERR)

         
c        ==== Calculate the determinant of the other spin block ====
         if (.not. ma_alloc_get(mt_dbl, nel(js)**2, 'ovl(js)',
     &        l_ovs, k_ovs))
     &        call errquit('WFAT: Cannot allocate ovl(js)', 58, MA_ERR)
         call dcopy(nel(js)**2, ovl(1:nel(js),1:nel(js),js), 1,
     &        dbl_mb(k_ovs), 1)
         call dgefa(dbl_mb(k_ovs), nel(js), nel(js), int_mb(k_iv), info)
         if (info /= 0)
     &        call errquit('WFAT: Error in dgefa', 46, UERR)
         call dgedi(dbl_mb(k_ovs), nel(js), nel(js), int_mb(k_iv), det0,
     &              dbl_mb(k_tm), 10)
         det = det0(1) * 10.0d0**det0(2)

         
c        ==== Deallocation ====
         if (.not. MA_free_heap(l_iv))
     &        call errquit('WFAT: ipvt deallocation failed', 7, MA_ERR)
         if (.not. MA_free_heap(l_tm))
     &        call errquit('WFAT: tmp deallocation failed', 7, MA_ERR)
         if (.not. MA_free_heap(l_ovs))
     &        call errquit('WFAT: ovl(js) deallocation failed',
     &        7, MA_ERR)

         
c        ==== Accumulate ====
         dsl = dsl + det * h
      enddo
      
      end
c=================================================
      
      
c=================================================
c=================================================
      subroutine mewfat_ovlslater(nb, nel, mo1, mo2, ovl, osl)
c     DESCRIPTION:
c     Calculates <Sl(MO1)|Sl(MO2)>
c     where Sl(MO) is a Slater determinant consisting of molecular orbitals MO.
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nb
      integer nel(2)
      double precision mo1(nb,nb,2)   ! [input], AO coefficients of MO1.
      double precision mo2(nb,nb,2)   ! [input], AO coefficients of MO2.
      double precision ovl(nb,nb,2)   ! [input], overlap matrix between MO1 (row) and MO2 (col).
      double precision osl            ! [output], overlap matrix elements between the two determinants represented by MO1 and MO2.

      integer is, l_ovs, k_ovs, l_iv, k_iv, info, niv, ntmp, l_tm, k_tm
      double precision det0(2), det(2)

      
      niv = maxval(nel)
      if (.not. ma_alloc_get(mt_int, niv, 'ipvt', l_iv, k_iv))
     &     call errquit('WFAT: Cannot allocate ipvt', 23, MA_ERR)
      ntmp = niv
      if (.not. ma_alloc_get(mt_dbl, ntmp, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 24, MA_ERR)
         
      do is = 1, 2
         if (.not. ma_alloc_get(mt_dbl, nel(is)**2, 'ovl(is)',
     &        l_ovs, k_ovs))
     &        call errquit('WFAT: Cannot allocate ovl(is)', 58, MA_ERR)
         call dcopy(nel(is)**2, ovl(1:nel(is),1:nel(is),is), 1,
     &        dbl_mb(k_ovs), 1)
         call dgefa(dbl_mb(k_ovs), nel(is), nel(is), int_mb(k_iv), info)
         if (info /= 0)
     &        call errquit('WFAT: Error in dgefa', 46, UERR)
         call dgedi(dbl_mb(k_ovs), nel(is), nel(is), int_mb(k_iv), det0,
     &              dbl_mb(k_tm), 10)
         det(is) = det0(1) * 10.0d0**det0(2)
         if (.not. MA_free_heap(l_ovs))
     &        call errquit('WFAT: ovl(js) deallocation failed',
     &        7, MA_ERR)
      enddo

      osl = det(1) * det(2)

      if (.not. MA_free_heap(l_iv))
     &     call errquit('WFAT: ipvt deallocation failed', 11, MA_ERR)
      if (.not. MA_free_heap(l_tm))
     &     call errquit('WFAT: tmp deallocation failed', 11, MA_ERR)
      
      end
      
