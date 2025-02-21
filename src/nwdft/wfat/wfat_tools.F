c=======================================================
c=======================================================
      subroutine wfat_matrix_sumabsorb(m, n, a, b, f, g)

      implicit none
      integer m               ! [input]
      integer n               ! [input]
      double precision a(m,n) ! [in/output]
      double precision b(m,n) ! [input]
      double precision f      ! [input]
      double precision g      ! [input]

      a = f*a + g*b
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dtranspose(m, n, a, b)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"      
      integer m, n
      double precision a(m,n)
      double precision b(n,m)

      integer l_tm, k_tm
      integer i, i_tm


      if (.not. ma_alloc_get(mt_dbl, n*m, 'tmp', l_tm, k_tm))      ! 1)
     &     call errquit('WFAT: Cannot allocate tmp in wfat_dtranspose',
     &     1, MA_ERR)
      do i = 1, m
         i_tm = k_tm + (i-1)*n
         dbl_mb(i_tm:i_tm+n-1) = a(i,1:n)
      enddo

      do i = 1, m
         i_tm = k_tm + (i-1)*n
         b(1:n,i) = dbl_mb(i_tm:i_tm+n-1)
      enddo
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp in '//
     &     'wfat_dtranspose', 1, MA_ERR)

      end
      
c     NOTES:
c       1) Temporary array is needed to anticipate a situation where a and b 
c          point to the same memory address, i.e. when a and b in the calling
c          procedure are the same array.
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dsyev_wrap(n, a, w)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      integer n                 ! [input] order of matrix a 
      double precision a(n,n)   ! [in/output] On input, the matrix to be diagonalized. On output, the eigenvectors.
      double precision w(n)     ! [output] eigenvalues

      integer lwork, info
      integer l_wo, k_wo
      

c     ==== Workspace query ====
      if (.not. ma_alloc_get(MT_DBL, 1, 'work', l_wo, k_wo))
     &     call errquit('WFAT: Cannot allocate work for workspace '//
     &     'query step of wfat_dsyev_wrap.', 1, UERR)

      call dsyev('v', 'u', n, a, n, w, dbl_mb(k_wo), -1, info)
      lwork = int( dbl_mb(k_wo) )

      if (.not. ma_free_heap(l_wo))
     &     call errquit('WFAT: Cannot deallocate work in workspace '//
     &     'query step of wfat_dsyev_wrap.', 1, UERR)


c     ==== Do the diagonalization ====
      if (.not. ma_alloc_get(MT_DBL, lwork, 'work', l_wo, k_wo))
     &     call errquit('WFAT: Cannot allocate work for the '//
     &     'solving step of wfat_dsyev_wrap.', 1, UERR)
      
      call dsyev('v', 'u', n, a, n, w, dbl_mb(k_wo), lwork, info)

      if (info /= 0)
     &     call errquit('WFAT: Error when solving for the solution '//
     &     'using wfat_dsyev_wrap.', 1, UERR)
      
      if (.not. ma_free_heap(l_wo))
     &     call errquit('WFAT: Cannot deallocate work in the '//
     &     'solving step of wfat_dsyev_wrap.', 1, UERR)

      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dsymv_wrap(n, alpha, a, x, beta, y)

      implicit none
      integer n
      double precision alpha
      double precision a(n,n)
      double precision x(n)
      double precision beta
      double precision y(n)


      call dsymv('u', n, alpha, a, n, x, 1, beta, y, 1)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dgemm_wrap(ta, tb, ma, na, mb, nb, alpha, a, b,
     &     beta, c)

      implicit none
#include "errquit.fh"
      character*1 ta
      character*1 tb
      integer ma
      integer na
      integer mb
      integer nb
      double precision alpha
      double precision a(ma,na)
      double precision b(mb,nb)
      double precision beta
      double precision c(1:*)
      
      integer m, n, k, lda, ldb, ldc
      

      if (ta /= 'n' .and. ta /= 't')
     &     call errquit('wfat_dgemm_wrap: Unknown value of argument ta',
     &     1, UERR)
      if (tb /= 'n' .and. tb /= 't')
     &     call errquit('wfat_dgemm_wrap: Unknown value of argument tb',
     &     1, UERR)

c     m = no. of rows of op(A)
c     n = no. of columns of op(B)
c     k = no. of columns of op(A) = no. of rows of op(B)
      
      if (ta == 'n') then
         m = ma
         k = na
      else
         m = na
         k = ma
      endif

      if (tb == 'n') then
         n = nb
         if (k /= mb)
     &        call errquit('wfat_dgemm_wrap: Inconsistent array sizes',
     &        1, UERR)
      else
         n = mb
         if (k /= nb)
     &        call errquit('wfat_dgemm_wrap: Inconsistent array sizes',
     &        1, UERR)
      endif

      
      if (ta == 'n') then
         lda = max(1,m)
      else
         lda = max(1,k)
      endif

      if (tb == 'n') then
         ldb = max(1,k)
      else
         ldb = max(1,n)
      endif
      
      ldc = max(1,m)

      
      call dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta,
     &     c, ldc)

c     The codes above are written such that the following conditions are always satisfied:
c       - When ta='n', LDA=ma and ka=k=na. When ta='t', LDA=ma and ka=m=na.
c       - When tb='n', LDB=mb and kb=n=nb. When tb='t', LDB=mb and kb=k=nb.
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dsymm_wrap(m, n, alpha, a, b, beta, c)
c     c = alpha*a*b + beta*c
      implicit none
      integer m                   ! [input], size of a and no. of rows of b
      integer n                   ! [input], no. of columns of b
      double precision alpha      ! [input], 
      double precision a(m,m)     ! [input], symmetric matrix
      double precision b(m,n)     ! [input], any matrix
      double precision beta       ! [input], 
      double precision c(m,n)     ! [output],   


      call dsymm('l', 'u', m, n, alpha, a, m, b, m, beta, c, m)


      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dsysv_wrap(n, nrhs, a, b)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"      
      integer n
      integer nrhs
      double precision a(n,n)
      double precision b(n,nrhs)

      integer l_wo, k_wo
      integer lwork, ipiv(n), info
      

c     ==== Workspace query ====
      if (.not. ma_alloc_get(MT_DBL, 1, 'work', l_wo, k_wo))
     &     call errquit('WFAT: Cannot allocate work for workspace '//
     &     'query step of wfat_dsysv_wrap.', 1, UERR)

      call dsysv('u', n, nrhs, a, n, ipiv, b, n, dbl_mb(k_wo), -1, 
     &     info)

c      if (.not. info /= 0)
c     &     call errquit('WFAT: Error during workspace query in '//
c     &     'wfat_dsysv_wrap.', 1, UERR)
      lwork = int( dbl_mb(k_wo) )

      if (.not. ma_free_heap(l_wo))
     &     call errquit('WFAT: Cannot deallocate work in workspace '//
     &     'query step of wfat_dsysv_wrap.', 1, UERR)

      
c     ==== Solving the matrix equation ====
      if (.not. ma_alloc_get(MT_DBL, lwork, 'work', l_wo, k_wo))
     &     call errquit('WFAT: Cannot allocate work for the '//
     &     'solving step of wfat_dsysv_wrap.', 1, UERR)
      
      call dsysv('u', n, nrhs, a, n, ipiv, b, n, dbl_mb(k_wo), lwork,
     &     info)

      if (info /= 0)
     &     call errquit('WFAT: Error when solving for the solution '//
     &     'using wfat_dsysv_wrap.', 1, UERR)
      
      if (.not. ma_free_heap(l_wo))
     &     call errquit('WFAT: Cannot deallocate work in the '//
     &     'solving step of wfat_dsysv_wrap.', 1, UERR)

      end
c=======================================================            

      
c=======================================================
c=======================================================
      subroutine wfat_trace(n, a, tr)

      implicit none
      integer n
      double precision a(n,n)
      double precision tr

      integer i
      
      tr = 0.0d0
      do i = 1, n
         tr = tr + a(i,i)
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_project_mat(m, n, v, p)
c     DESCRIPTION :
c     Return the matrix of projection operator defined by
c        p = \sum_{i=1}^n  v(:,i) * v(:,i)^T
c     where the * sign indicates matrix product.
      
      implicit none
      integer m
      integer n
      double precision v(m,n)
      double precision p(m,m)


      p = matmul(v, transpose(v))

      end
c=======================================================      
      

c=======================================================
c=======================================================
      subroutine wfat_ne_expval_1eopr(nb, ns, nel, vcf, opr,
     &           val1e, valne)

      implicit none
      integer nb
      integer ns
      integer nel(ns)
      double precision vcf(nb, nb, ns)
      double precision opr(nb, nb)
      double precision val1e(nb, ns)
      double precision valne

      integer is, io
      

      valne = 0.0d0
      do is = 1, ns
         do io = 1, nb
            call wfat_1e_expval_1eopr(nb, vcf(:,io,is), opr,
     &                            val1e(io,is))
            if (io <= nel(is)) valne = valne + val1e(io,is)
         enddo
      enddo
      
      if (ns == 1) valne = 2 * valne

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_1e_expval_1eopr(nb, vcf, opr, val1e)

      implicit none
      integer nb
      double precision vcf(nb)
      double precision opr(nb, nb)     ! [input], opr must be a symmetric matrix.
      double precision val1e

      integer i, j, f
      

      val1e = 0.0d0
      do j = 1, nb
         do i = 1, j
            if (i /= j) then
               f = 2
            else
               f = 1
            endif

            val1e = val1e + f * vcf(j) * vcf(i) * opr(i,j)
         enddo
      enddo
      

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_mat_element(n, u, v, o, val)
c     val = <u|o|v>
c     o must be a symmetric matrix.
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer n
      double precision u(n)
      double precision v(n)
      double precision o(n,n)
      double precision val

      integer l_sc, k_sc


      if (.not. ma_alloc_get(mt_dbl, n, 'scr', l_sc, k_sc))
     &     call errquit('WFAT: Cannot allocate scr in wfat_mat_element',
     &     1, MA_ERR)
      
      call wfat_dsymv_wrap(n, 1.0d0, o, v, 0.0d0, dbl_mb(k_sc))
      val = dot_product(u, dbl_mb(k_sc:k_sc+n-1))

      if (.not. ma_free_heap(l_sc))
     &     call errquit('WFAT: Cannot deallocate scr in '//
     &     'wfat_mat_element', 1, MA_ERR )

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_ga_get_lo(g_a, vilo, vjlo)

      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      integer g_a
      integer vilo(*)
      integer vjlo(*)

      integer nbl, nblock(2), ip, ilo, ihi, jlo, jhi, g_tmp, i, j

      ip = ga_nodeid()
      call ga_nblock(g_a, nblock)
      call ga_distribution(g_a, ip, ilo, ihi, jlo, jhi)

      nbl = maxval(nblock)
      if (.not. ga_create(mt_int, 1, nbl, 'tmp', 1, nbl, g_tmp))
     &     call errquit('WFAT: Cannot create tmp global array.', 
     &     1, GA_ERR)

      call ga_zero(g_tmp)
      vilo(1:nblock(1)) = 0
      do i = 1, nblock(1)
         if (ip == i - 1) call ga_put(g_tmp, 1, 1, i, i, ilo, 1)
      enddo
      call ga_sync()
      call ga_get(g_tmp, 1, 1, 1, nblock(1), vilo, 1)

      call ga_zero(g_tmp)
      vjlo(1:nblock(2)) = 0
      j = 1
      do i = 1, nblock(1)*nblock(2), nblock(1)
         if (ip == i - 1) call ga_put(g_tmp, 1, 1, j, j, jlo, 1)
         j = j + 1
      enddo
      call ga_sync()
      call ga_get(g_tmp, 1, 1, 1, nblock(2), vjlo, 1)

      if (.not. ga_destroy(g_tmp))
     &     call errquit('WFAT: Cannot destroy tmp.', 21, GA_ERR)
      call ga_sync()

      end
c=======================================================


c===============================================================
c===============================================================
      subroutine wfat_ga_mat_element(g_v, g_w, g_a, d)
c  d = < g_v | g_a | g_w >
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"      
      integer g_v
      integer g_w
      integer g_a
      double precision d

      logical l
      integer typev, typew, typea
      integer m_v, n_v, m_w, n_w, m_a, n_a
      integer g_tmp


      call ga_inquire(g_v, typev, m_v, n_v)
      call ga_inquire(g_w, typew, m_w, n_w)
      call ga_inquire(g_a, typea, m_a, n_a)


      if (typev == typew) then
         if (typew /= typea) then
            call errquit('WFAT: Type of argument ''g_w'' and of '//
     &           'argument ''g_a'' of ''wfat_ga_vrepresent'' must be '//
     &           'the same.', 1, UERR)
         endif
      else
         call errquit('WFAT: Type of argument ''g_v'' and of '//
     &        'argument ''g_w'' of ''wfat_ga_vrepresent'' must be '//
     &        'the same.', 1, UERR)
      endif
      

      if (n_w /= 1) then
         call errquit('WFAT: Argument ''g_w'' of '//
     &        '''wfat_ga_mat_element'' must be a column vector.',
     &        1, UERR)
      endif
      if (.not. (m_v == 1 .or. n_v == 1)) then
         call errquit('WFAT: Argument ''g_v'' of '//
     &        '''wfat_ga_mat_element'' must be a 1D array, '//
     &        'i.e. a column or row vector.', 1, UERR)
      endif
      
      if (n_a /= m_w) then
         call errquit('WFAT: In ''wfat_ga_mat_element'', the number '//
     &    'of columns of ''g_a'' must be equal to the number of '//
     &    'rows of ''g_w''.', 1, UERR)
      endif

      if (m_v*n_v /= m_a) then
         call errquit('WFAT: In ''wfat_ga_mat_element'', the number '//
     &        'of elements of ''g_v'' must be equal to the number of '//
     &        'rows of ''g_a''.', 1, UERR)
      endif
      

         

      if (m_a == n_a) then
         l = ga_duplicate(g_w, g_tmp, 'tmp')
      else
         l = ga_create(mt_dbl, m_a, 1, 'tmp', 0, 0, g_tmp)
      endif
      if (.not. l)
     &     call errquit('wfat_ga_mat_element: Cannot create tmp.',
     &     1, GA_ERR)
      call ga_zero(g_tmp)

      
c     ==== g_a * g_w ====
      call ga_dgemm('n', 'n', m_a, n_w, m_w, 1.0d0, g_a, g_w, 0.0d0,
     &     g_tmp)

c     ==== g_v * (g_a * g_w) ====
      d = ga_ddot(g_v, g_tmp)

      
      if (.not. ga_destroy(g_tmp))
     &     call errquit('wfat_ga_vrepresent: Cannot destroy tmp.',
     &     1, GA_ERR)
      
      end
c===============================================================
      

c===============================================================
c===============================================================
      subroutine wfat_ga_vrepresent(g_v, g_a, g_b)
c   g_b = g_v^T . g_a . g_v
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"      
      integer g_v     ! [input]
      integer g_a     ! [input]
      integer g_b     ! [output]

      integer typev, m_v, n_v, typea, m_a, n_a, typeb, m_b, n_b
      integer g_tmp
      
      
      call ga_inquire(g_v, typev, m_v, n_v)
      call ga_inquire(g_a, typea, m_a, n_a)
      call ga_inquire(g_b, typeb, m_b, n_b)

      if (typea == typeb) then
         if (typeb /= typev) then
            call errquit('WFAT: Type of argument g_b and of '//
     &           'argument g_v of wfat_ga_vrepresent must be the same.',
     &           1, UERR)
         endif
      else
         call errquit('WFAT: Type of argument g_a and of argument '//
     &        'g_b of wfat_ga_vrepresent must be the same.', 1, UERR)
      endif
      
      if (n_a /= m_v) then
         call errquit('WFAT: The no. of columns of the matrix '//
     &        'represented by argument g_a of wfat_ga_vrepresent '//
     &        'must be equal the no. of rows of the matrix '//
     &        'represented by g_v', 1, GA_ERR)
      endif
      if (m_a /= n_a) then
         call errquit('WFAT: The matrix represented by argument g_a '//
     &        'of wfat_ga_vrepresent must be square', 1, GA_ERR)
      endif
      if (m_b /= n_b) then
         call errquit('WFAT: The matrix represented by argument g_b '//
     &        'of wfat_ga_vrepresent must be square', 1, GA_ERR)
      endif
      if (n_b /= n_v) then
         call errquit('WFAT: The no. of columns of the matrix '//
     &        'represented by argument g_b of wfat_ga_vrepresent '//
     &        'must be equal the no. of columns of the matrix '//
     &        'represented by g_v', 1, GA_ERR)
      endif
      

      if (.not. ga_create(mt_dbl, m_a, n_v, 'tmp', 0, 0, g_tmp))
     &     call errquit('wfat_ga_vrepresent: Cannot create tmp.',
     &     1, GA_ERR)
      call ga_zero(g_tmp)
      
      call ga_dgemm('n', 'n', m_a, n_v, n_a, 1.0d0, g_a, g_v, 0.0d0,
     &     g_tmp)
      call ga_dgemm('t', 'n', m_b, n_b, m_a, 1.0d0, g_v, g_tmp, 0.0d0,
     &     g_b)

      if (.not. ga_destroy(g_tmp))
     &     call errquit('wfat_ga_vrepresent: Cannot destroy tmp.',
     &     1, GA_ERR)

      end
c===============================================================
      

c===============================================================
c===============================================================
      subroutine wfat_represent_sym(m, n, v, a, b)
c     b = v^T . a . v
c     where a is a symmetric matrix
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer m
      integer n
      double precision v(m,n)
      double precision a(m,m)       ! a is symmetric
      double precision b(n,n)

      integer l_tm, k_tm
      

      if (.not. ma_alloc_get(MT_DBL, m*n, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
      
      call wfat_dsymm_wrap(m, n, 1.0d0, a, v, 0.0d0, dbl_mb(k_tm))
      call wfat_dgemm_wrap('t', 'n', m, n, m, n, 1.0d0, v, dbl_mb(k_tm),
     &     0.0d0, b)

      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_represent_diag(m, n, v, a, b)
c     b = v^T . a . v
c     where a is a diagonal matrix.
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer m
      integer n
      double precision v(m,n)
      double precision a(m)       ! a is symmetric
      double precision b(n,n)

      integer l_tm, k_tm, i, j, i_tm
      

      if (.not. ma_alloc_get(MT_DBL, m*n, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)

      do j = 1, n
         do i = 1, m
            i_tm = k_tm + (j-1)*m + (i-1)
            dbl_mb(i_tm) = a(i) * v(i,j)
         enddo
      enddo

      call wfat_dgemm_wrap('t', 'n', m, n, m, n, 1.0d0, v, dbl_mb(k_tm),
     &     0.0d0, b)

      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_represent(m, n, v, a, t, b)

      implicit none
      integer m
      integer n
      double precision v(m,n)
      double precision a(m,m)
      double precision t(m,n)
      double precision b(n,n)

      t = matmul(a, v)
      b = matmul(transpose(v), t)

      end
c===============================================================
      
            
c=======================================================
c=======================================================
      subroutine wfat_get_i2cart(ic, ic1, ic2)
c     Based on subroutine xc_drhonuc in $NWCHEM_TOP/src/nwdft/xc/xc_rhogen.F
      
      implicit none
#include "errquit.fh"
      integer ic
      integer ic1, ic2


      if (ic == 1) then
         ic1 = 1
         ic2 = 1
      else if (ic == 2) then
         ic1 = 1
         ic2 = 2
      else if (ic == 3) then 
         ic1 = 1
         ic2 = 3
      else if (ic == 4) then 
         ic1 = 2
         ic2 = 2
      else if (ic == 5) then 
         ic1 = 2
         ic2 = 3
      else if (ic == 6) then 
         ic1 = 3
         ic2 = 3
      else
         call errquit('WFAT: The chosen value for argument ic '//
     &        'of wfat_get_i2cart is unavailable.', 1, UERR)
      endif




cOLD      if (ic == 1) then
cOLD         ic1 = 1
cOLD         ic2 = 1
cOLD      else if (ic == 2) then
cOLD         ic1 = 2
cOLD         ic2 = 2
cOLD      else if (ic == 3) then 
cOLD         ic1 = 3
cOLD         ic2 = 3
cOLD      else if (ic == 4) then 
cOLD         ic1 = 1
cOLD         ic2 = 2
cOLD      else if (ic == 5) then 
cOLD         ic1 = 2
cOLD         ic2 = 3
cOLD      else if (ic == 6) then 
cOLD         ic1 = 1
cOLD         ic2 = 3
cOLD      else
cOLD         call errquit('WFAT: The chosen value for argument ic '//
cOLD     &        'of wfat_get_i2cart is unavailable.', 1, UERR)
cOLD      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_cart2i(ic1, ic2, ic)

      implicit none
#include "errquit.fh"
      integer ic1, ic2
      integer ic      



      if (ic1 == 1 .and. ic2 == 1) then
         ic = 1          
      else if (ic1 == 2 .and. ic2 == 2) then
         ic = 4          
      else if (ic1 == 3 .and. ic2 == 3) then
         ic = 6           
      else if ((ic1 == 1 .and. ic2 == 2) .or.
     &         (ic1 == 2 .and. ic2 == 1)) then
         ic = 2           
      else if ((ic1 == 2 .and. ic2 == 3) .or.
     &         (ic1 == 3 .and. ic2 == 2)) then
         ic = 5           
      else if ((ic1 == 1 .and. ic2 == 3) .or.
     &         (ic1 == 3 .and. ic2 == 1)) then
         ic = 3 
      else
         call errquit('WFAT: The chosen values for arguments ic1 '//
     &        'and ic2 of wfat_get_cart2i are unavailable.', 1, UERR)
      endif




      
cOLD      if (ic1 == 1 .and. ic2 == 1) then
cOLD         ic = 1          
cOLD      else if (ic1 == 2 .and. ic2 == 2) then
cOLD         ic = 2          
cOLD      else if (ic1 == 3 .and. ic2 == 3) then
cOLD         ic = 3           
cOLD      else if ((ic1 == 1 .and. ic2 == 2) .or.
cOLD     &         (ic1 == 2 .and. ic2 == 1)) then
cOLD         ic = 4           
cOLD      else if ((ic1 == 2 .and. ic2 == 3) .or.
cOLD     &         (ic1 == 3 .and. ic2 == 2)) then
cOLD         ic = 5           
cOLD      else if ((ic1 == 1 .and. ic2 == 3) .or.
cOLD     &         (ic1 == 3 .and. ic2 == 1)) then
cOLD         ic = 6 
cOLD      else
cOLD         call errquit('WFAT: The chosen values for arguments ic1 '//
cOLD     &        'and ic2 of wfat_get_cart2i are unavailable.', 1, UERR)
cOLD      endif

      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_get_digits(val, dig)

      implicit none
      integer val
      integer dig

      integer aval, i, rem, nmax
      parameter (nmax=10)
      

      aval = abs(val)

      do i = 1, nmax
         rem = mod(aval, 10**i)
         if (rem == aval) then
            dig = i
            exit
         endif
      enddo

      if (val < 0) dig = dig + 1

      end
c=======================================================      


c=======================================================
c=======================================================
      subroutine wfat_dprintmat(m, n, e, a)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
      integer m
      integer n
      character*1 e
      double precision a(m,n)

      character*6 nn
      character*30 efmt
      character*7 efmt0
      integer i


      if (n <= 0 .or. n >= 1000000) then
         call errquit('WFAT: The value of argument n of '//
     &        'wfat_dprintmat must be positive and less than '//
     &        '1000000.', 1, UERR)
      endif

      write(nn, '(i6)') n
      if (e == 'f') then
         efmt0 = "(i10)"
         efmt = "(i5, "//trim(adjustl(nn))//"f10.5)"
      else if (e == 'e') then
         efmt0 = "(i14)"
         efmt = "(i5, "//trim(adjustl(nn))//"e14.5)"
      else
         call errquit('WFAT: The choice of argument e of '//
     &                'wfat_dprintmat is unavailable.', 1, UERR)
      endif

      write(luout, '(5x)', advance='no')
      call flush(luout)
      do i = 1, n
         write(luout, efmt0, advance='no') i
         call flush(luout)
      enddo
      write(luout, '(1x)')
      
      do i = 1, m
         write(luout, efmt) i, a(i, 1:n)
      enddo

      write(luout, '(" ")') 
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_repeat(iunit, n, os, c)

      implicit none
#include "errquit.fh"      
      integer iunit, n, os
      character*1 c

      integer i


      if (os > 100)
     &     call errquit('WFAT: Argument os of wfat_print_repeat '//
     &     'cannot be larger than 100.', 1, UERR)

      call util_flush(iunit)
      do i = 1, os
         write(iunit, '(" ")', advance='no')
         call util_flush(iunit)
      enddo

      do i = 1, n
         write(iunit, '(a1)', advance='no') c
         call util_flush(iunit)
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_tf_to_yn(l, e, a)

      implicit none
      logical l       ! [input], input logical to be converted.
      logical e       ! [input], .true. for 'yes/no', .false. for 'y/n'.
      character*(*) a ! [output], result of conversion.


      if (l) then
         if (e)       a = 'yes'
         if (.not. e) a = 'y'
      else
         if (e)       a = 'no'
         if (.not. e) a = 'n'
      endif

      end
c=======================================================            


c=======================================================
c=======================================================
      subroutine wfat_def_intglname(rtdb, lvl, fname)

      implicit none
#include "errquit.fh"
      integer rtdb
      character*2 lvl
      character*(*) fname


      if (lvl == 'oe') then
         call util_file_prefix('owigl', fname)
      else if (lvl == 'me') then
         call util_file_prefix('mwigl', fname)
      else
         call errquit('WFAT: Unknown value of argument ''lvl'' '//
     &        'of ''wfat_def_intglname''', 1, UERR)
      endif
      
      call wfat_prepend_rtdbpath(rtdb, fname)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_def_movecname(lvl, rtdb, fname)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      character*(*) lvl
      integer rtdb
      character*(*) fname


      if (lvl == 'oe') then
         call util_file_prefix('movecs', fname)
      else if (lvl == 'me0') then
         call util_file_prefix('ion0.movecs', fname)
      else if (lvl == 'me1') then
         call util_file_prefix('ion1.movecs', fname)
      else
         call errquit('WFAT: Argument ''lvl'' of wfat_def_movecname '//
     &        'is unknown.', 1, UERR)
      endif
      call wfat_prepend_rtdbpath(rtdb, fname)

      end
c=======================================================      

      
c=======================================================
c=======================================================
      subroutine wfat_prepend_rtdbpath(rtdb, fname)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"      
      integer rtdb              ! [input]
      character*(*) fname       ! [input/output]

      character*256 rtdb_name
      integer rtlen, rtlen0, nmlen0, i, pos



      if(.not. rtdb_getfname(rtdb, rtdb_name)) 
     &     call errquit('WFAT: Cannot get RTDB file name', 1, RTDB_ERR)
      
      rtlen = inp_strlen(rtdb_name)
      pos = 0
      do i = 1, rtlen
         if (rtdb_name(i:i) == '/' ) pos = i
      enddo

      rtlen0 = inp_strlen(trim(adjustl(rtdb_name(1:pos))))
      nmlen0 = inp_strlen(trim(adjustl(fname)))
      if (rtlen0 + nmlen0 > len(fname))
     &     call errquit('WFAT: The argument ''fname'' of '//
     &     'wfat_prepend_rtdbpath subroutine is too small to hold '//
     &     'the output string where the path of RTDB will have '//
     &     'been prepended.', 1, UERR)
      fname = trim(adjustl(rtdb_name(1:pos)))//adjustl(fname)
      
cforcheck      write(luout, '("Vector file : ", a256)') adjustl(fname)

      end
c=======================================================


c=======================================================
c=======================================================
      integer function delkron(i,j)
c     DESCRIPTION:
c     Kronecker delta function for integer arguments.
      
      implicit none
      integer i, j

      if (i == j) then
         delkron = 1
      else
         delkron = 0
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_identity(n, e)

      implicit none
      integer n
      double precision e(n,n)

      integer i

      e = 0.0d0
      do i = 1, n
         e(i,i) = 1.0d0
      enddo

      end
c=======================================================

      
c=======================================================
c=======================================================
      double precision function double_fact(n)

      implicit none
      integer n

      integer k
      double precision f
      

      if (n >= 0) then
         if (mod(n,2) == 0) then
            f = 1.0d0
            do k = 1, n/2
               f = f * (2*k)
            enddo
         else if (mod(n,2) == 1) then
            f = 1.0d0
            do k = 1, (n+1)/2
               f = f * (2*k-1)
            enddo
         endif
      else
         if (mod(abs(n),2) == 1) then  ! negative odd integer n.
            f = 1.0d0
            do k = 1, (n+1)/2
               f = f * (2*k-1)
            enddo
            f = (-1)**((abs(n)-1)/2) * abs(n) / f
         endif
      endif
      double_fact = f

      end
c=======================================================

      
c=======================================================
c=======================================================
      double precision function rise_fact(x, n)
      implicit none
      double precision x
      integer n
      
      integer i
      double precision f


      if (n == 0) then
         rise_fact = 1.0d0
      else
         f = 1.0d0
         do i = 1, n
            f = f * (x+i-1)
         enddo
         rise_fact = f
      endif
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      double precision function binom_cf(p, n)
c     DESCRIPTION:
C     This subroutine calculates the binomial coefficient defined in eq. (6)
c     in Ref. [1].
      implicit none
#include "errquit.fh"
#include "stdio.fh"
      double precision p
      integer n

      integer i
      double precision f


c      write(luout, '("bn1")') 
      if (n == 0) then
         binom_cf = 1.0d0
      else if (n > 0) then
         f = 1.0d0
         do i = 1, n
c            write(luout, '("bn2", i6)') i
            f = f * (p-i+1)/dble(i)
         enddo
         binom_cf = f
      else
         call errquit('ERROR: The argument n of bnom_cf cannot '//
     &        'be negative' , 90, UERR)
      endif

c     REFERENCES:
c     [1] T. M. Rassias and H. M. Srivastava, The Orthogonality Property of the 
c         Classical Laguerre Polynomials, APPLIED MATHEMATICS AND COMPUTATION 50:167-173
c         (1992).
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_legendre(ut,N,x,L)
      implicit none
#include "errquit.fh"      
      character*1 ut
      integer N
      double precision x
      double precision L
      
      integer i
      double precision norm, L0, L1

      
      if (ut == 'y') then
         norm = sqrt((2*N+1.0d0)/2)
      else if (ut == 'n') then
         norm = 1.0d0
      else
         call errquit('WFAT: The argument ut for wfat_legendre '//
     &        'must be either ''y'' or ''n''.', 18, UERR)
      end if
      
      
      L0 = 1.0d0
      L1 = x
      
      if (N == 0) then
         L = L0
         return
      else if (N == 1) then
         L = L1
         return
      end if
      
c     Calculate the N-th order Legendre polynomial using the recursive formula
c     N.L_N(x) = (2N-1).x.L_(N-1)(x) - (N-1).L_(N-2)(x)
      do i = 2,N
         L = ((2*i-1)*x*L1 - (i-1)*L0)/i
         L0 = L1
         L1 = L
      end do
      
      L = norm*L

      end
c=======================================================
      

c============================================================
c============================================================
      subroutine wfat_submatrix(mm, nn, a, m, n, irow, icol, b)

      implicit none
      integer mm
      integer nn
      double precision a(mm,nn)
      integer m
      integer n
      integer irow(m)           ! [input], 1 <= irow(i) <= mm
      integer icol(n)           ! [input], 1 <= icol(i) <= nn
      double precision b(m,n)

      integer j, i

      do j = 1, n
         do i = 1, m
            b(i,j) = a( irow(i), icol(j) )
         enddo
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_kron_prod(ma, na, mb, nb, a, b, c)

      implicit none
      integer ma
      integer na
      integer mb
      integer nb
      double precision a(ma,na)
      double precision b(mb,nb)
      double precision c(ma*mb, na*nb)

      integer m, n, i, j, ic1, ic2, jc1, jc2

      m = ma * mb
      n = na * nb
      do j = 1, na
         jc1 = (j-1)*nb + 1
         jc2 = jc1 + nb - 1
         do i = 1, ma
            ic1 = (i-1)*mb + 1
            ic2 = ic1 + mb - 1
            
            c(ic1:ic2, jc1:jc2) = a(i,j) * b
         enddo
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_ga_kron_prod(ma, na, mb, nb, a, b, g_c)

      implicit none
#include "global.fh"
      integer ma
      integer na
      integer mb
      integer nb
      double precision a(ma,na)
      double precision b(mb,nb)
      integer g_c

      integer ip, mp, ilo, ihi, jlo, jhi, i, j, ia, ja, ib, jb
      
      call ga_sync()
      
      ip = ga_nodeid()
      call ga_distribution(g_c, ip, ilo, ihi, jlo, jhi)
      mp = ihi - ilo + 1

      do j = jlo, jhi
         do i = ilo, ihi
c           Given (i,j), find ia, ja, ib, and jb so that c(i,j) = a(ia,ja)*b(ib,jb).
            call wfat_invert_kron_id(i, mb, ia, ib)
            call wfat_invert_kron_id(j, nb, ja, jb)
            call ga_put(g_c, i, i, j, j, a(ia,ja) * b(ib,jb), 1)
         enddo
      enddo
      
      call ga_sync()
      
      end
c============================================================

      
c============================================================
c============================================================
      subroutine wfat_invert_kron_id(i, n, ia, ib)
c     If C = A \otimes B (C is equal to the kronecker product between
c     A and B), then C(i,j) = A(ia,ja)*B(ib,jb). This subroutine is
c     used to get ia and ib given i and the number of rows of B. It
c     can also be used to get ja and jb given j and the number of
c     columns of B.
      implicit none
      integer i      ! [input], index of C in the dimension of interest.
      integer n      ! [input], size of B in the dimension of interest.
      integer ia     ! [output], index of A in the dimension of interest.
      integer ib     ! [output], index of B in the dimension of interest.

      if (mod(i,n) == 0) then
         ia = i/n
      else
         ia = i/n + 1
      endif

      ib = mod(i-1,n) + 1

      end
c============================================================
      

c============================================================
c============================================================
      subroutine wfat_isdiag(n, a, thr, isdiag)

      implicit none
      integer n
      double precision a(n,n)
      double precision thr
      logical isdiag

      integer i, j

      isdiag = .true.
      do j = 1, n
         do i = 1, n
            if (i /= j .and. abs(a(i,j)) > thr) then
               isdiag = .false.
               return
            endif
         enddo
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_info_banner(u, p0)

      implicit none
#include "global.fh"
      integer u     ! [input], Print unit.
      logical p0    ! [input], If .True., only print on process 0, otherwise, print for all processes.

      if (p0) then
         if (ga_nodeid() == 0) then
            write(u, '("")')
            write(u, '(3x, "*** INFO ***")')
         endif
      else
         write(u, '("")')
         write(u, '(3x, "*** INFO ***")')
      endif

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_attention_banner(u, p0)

      implicit none
#include "global.fh"
      integer u     ! [input], Print unit.
      logical p0    ! [input], If .True., only print on process 0, otherwise, print for all processes.

      if (p0) then
         if (ga_nodeid() == 0) then
            write(u, '("")')
            write(u, '(3x, ">>> ATTENTION <<<")')
         endif
      else
         write(u, '("")')
         write(u, '(3x, ">>> ATTENTION <<<")')
      endif

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_num_derivative2(n, dx, f, d2f)

      implicit none
      integer n
      double precision dx
      double precision f(n)
      double precision d2f(n)

      integer i


      do i = 1, n
         if (i == 1) then
c           ==== 3-point forward ====
            d2f(i) = (f(i+2) - 2*f(i+1) + f(i)) / dx**2
         else if (i == 2 .or. i == n-1) then
c           ==== 3-point central ====
            d2f(i) = (f(i+1) - 2*f(i) + f(i-1)) / dx**2
         else if (i == n) then
c           ==== 3-point backward ====
            d2f(i) = (f(i) - 2*f(i-1) + f(i-2)) / dx**2
         else
c           ==== 5-point central ====
            d2f(i) = (-f(i+2) + 16*f(i+1) - 30*f(i) + 16*f(i-1) 
     &           - f(i-2)) / (12*dx**2)
         endif
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_read_channel_id(nch, oid, sid, nqn, nid, mid, s, 
     &     ir, ic, info)
c     DESCSRIPTION:
c     Reads off the element address of a channel (consisting of
c     a spin-orbital and parabolic quantum numbers) from a given
c     list of spin-orbitals and quantum numbers.
c     
c     Acceptable format for s is "A<delim>B<delim>C<delim>D" where A,
c     B, C, and D are spin ID, orbital ID, n parabolic quantum number,
c     and m parabolic quantum number, respectively. <delim> is the
c     delimiter, which must be either a blank space character or a comma.
      implicit none
#include "inp.fh"
#include "stdio.fh"
      integer nch         ! [input],
      integer oid(nch)    ! [input],
      integer sid(nch)    ! [input],
      integer nqn         ! [input],
      integer nid(nqn)    ! [input],
      integer mid(nqn)    ! [input],
      character(len=*) s  ! [input],
      integer ir          ! [output],
      integer ic          ! [output],
      integer info        ! [output],

      logical tt
      integer i1, i2, sid0, oid0, nid0, mid0, ios

      
      info = 0
      i1 = 0

      if (inp_strtok(s, ' ,', i1, i2)) then
         if (s(i1:i2) == 'a') then
            sid0 = 1
         else if (s(i1:i2) == 'b') then
            sid0 = 2
         else
            info = 1
         endif
      else
         info = 2
      endif

      if (inp_strtok(s, ' ,', i1, i2)) then
         read(s(i1:i2), '(i256)', iostat=ios) oid0
         if (ios > 0) info = 3
      else
         info = 4
      endif

      call wfat_get_orb_id(oid0, sid0, nch, oid, sid, ic)
      if (ic < 0) info = 5


      if (inp_strtok(s, ' ,', i1, i2)) then
         read(s(i1:i2), '(i256)', iostat=ios) nid0
         if (ios > 0 .or. nid0 < 0) info = 6
      else
         info = 7
      endif

      if (inp_strtok(s, ' ,', i1, i2)) then
         read(s(i1:i2), '(i256)', iostat=ios) mid0
         if (ios > 0) info = 8
      else
         info = 9
      endif

      call wfat_get_paraqn_id(nid0, mid0, nqn, nid, mid, ir)
      if (ir < 0) info = 10

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_read_orb_id(nch, oid, sid, s, ic, info)
c     DESCSRIPTION:
c     Reads off the element address of a spin-orbital from a given
c     list of spin-orbitals.
c     
c     Acceptable format for s is "A<delim>B<delim>C<delim>D" where A,
c     B, C, and D are spin ID, orbital ID, n parabolic quantum number,
c     and m parabolic quantum number, respectively. <delim> is the
c     delimiter, which must be either a blank space character or a comma.
      implicit none
#include "inp.fh"
#include "stdio.fh"
      integer nch         ! [input],
      integer oid(nch)    ! [input],
      integer sid(nch)    ! [input],
      character(len=*) s  ! [input],
      integer ic          ! [output],
      integer info        ! [output],

      logical tt
      integer i1, i2, sid0, oid0, ios

      
      info = 0
      i1 = 0

      if (inp_strtok(s, ' ,', i1, i2)) then
         if (s(i1:i2) == 'a') then
            sid0 = 1
         else if (s(i1:i2) == 'b') then
            sid0 = 2
         else
            info = 1
         endif
      else
         info = 2
      endif

      if (inp_strtok(s, ' ,', i1, i2)) then
         read(s(i1:i2), '(i256)', iostat=ios) oid0
         if (ios > 0) info = 3
      else
         info = 4
      endif

      call wfat_get_orb_id(oid0, sid0, nch, oid, sid, ic)
      if (ic < 0) info = 5

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_get_orb_id(oid0, sid0, nch, oid, sid, id)
c     DESCRIPTION:
c     Returns the element address of an orbital with orbital ID oid0 and
c     spin ID sid0.
      implicit none
      integer oid0      ! [input],
      integer sid0      ! [input],
      integer nch       ! [input],
      integer oid(nch)  ! [input],
      integer sid(nch)  ! [input],
      integer id        ! [output],

      integer i
      
      id = -1
      do i = 1, nch
         if (oid(i) == oid0 .and. sid(i) == sid0) then
            id = i
            exit
         endif
      enddo

      end
c============================================================
      
            
c============================================================
c============================================================
      subroutine wfat_get_paraqn_id(nid0, mid0, nqn, nid, mid, id)
c     DESCRIPTION:
c     Returns the element address of a parabolic channel with n quantum
c     number nid0 and m quantum number mid0.
      implicit none
      integer nid0       ! [input],
      integer mid0       ! [input],
      integer nqn        ! [input],
      integer nid(nqn)   ! [input],
      integer mid(nqn)   ! [input],
      integer id         ! [output],

      integer i
      
      id = -1
      do i = 1, nqn
         if (nid(i) == nid0 .and. mid(i) == mid0) then
            id = i
            exit
         endif
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_minor(k1, k2, m, n, a, b)
c     DESCRIPTION:
C     Returns the minor of a matrix.
      implicit none
#include "errquit.fh"
      integer k1                     ! [input], index of the removed row.
      integer k2                     ! [input], index of the removed column.
      integer m                      ! [input], number of rows of input matrix.
      integer n                      ! [input], number of columns of input matrix
      double precision a(m,n)        ! [input], input matrix.
      double precision b(m-1,n-1)    ! [output], minor of the input matrix.

      if (m == 1 .or. n == 1)
     &     call errquit('WFAT: Argument m and n of wfat_minor must be'//
     &     ' greater than 1.', 95, UERR)
      if (k1 < 1 .or. k1 > m)
     &     call errquit('WFAT: Argument k1 is out of the allowed '//
     &     'range.', 96, UERR)
      if (k2 < 1 .or. k2 > n)
     &     call errquit('WFAT: Argument k2 is out of the allowed '//
     &     'range.', 96, UERR)

c     ==== Upper left ====
      if (k1 > 1 .and. k2 > 1) b(1:k1-1, 1:k2-1) = a(1:k1-1, 1:k2-1)

c     ==== Upper right ====
      if (k1 > 1 .and. k2 < n) b(1:k1-1, k2:n-1) = a(1:k1-1, k2+1:n)

c     ==== Lower left ====
      if (k1 < m .and. k2 > 1) b(k1:m-1, 1:k2-1) = a(k1+1:m, 1:k2-1)

c     ==== Lower right ====
      if (k1 < m .and. k2 < n) b(k1:m-1, k2:n-1) = a(k1+1:m, k2+1:n)
      
      end
c============================================================
      
