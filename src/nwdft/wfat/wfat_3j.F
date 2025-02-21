!==================================================
!==================================================
      subroutine wfat_3j_jscan(tj1, tm1, tj2, tm2, cf)
c     DESCRIPTION:
c     Calculates the Wigner 3-j coefficients
c       / j1  j2      j3   \
c       \ m1  m2  -(m1+m2) /
c     for $j3 = [j_min, j1+j2]$, where
c     $j_min = max(|j1-j2|,|m1+m2|)$.
      implicit none
      integer tj1
      integer tm1
      integer tj2
      integer tm2
      double precision cf(1:*)

      integer i
      double precision j1, j2, j3, m1, m2, j_min, a1, a2, b

      
      j1 = tj1 / 2.0d0
      j2 = tj2 / 2.0d0
      m1 = tm1 / 2.0d0
      m2 = tm2 / 2.0d0

      
c     ==== Determine the type of the first 3-j symbol ====
      j_min = max(abs(j1-j2), abs(m1+m2))
      if ( idnint(2*j_min) == abs(tj1-tj2) ) then
c        / j1  j2   |j1-j2| \
c        \ m1  m2  -(m1+m2) /
         call wfat_3j_special_type1(tj1, tm1, tj2, tm2, cf(1))
         if (abs(tj1-tj2) == tj1+tj2) return
c        / j1  j2  |j1-j2|+1 \
c        \ m1  m2   -(m1+m2) /
         call wfat_3j_special_type3(tj1, tm1, tj2, tm2, cf(2))
      else if ( idnint(2*j_min) == abs(tm1+tm2) ) then
c        / j1  j2   |m1+m2| \
c        \ m1  m2  -(m1+m2) /
         call wfat_3j_special_type2(tj1, tm1, tj2, tm2, -(tm1+tm2),
     &        cf(1))
         if (abs(tm1+tm2) == tj1+tj2) return
c        / j1  j2   |m1+m2|+1 \
c        \ m1  m2    -(m1+m2) /
         call wfat_3j_special_type4(tj1, tm1, tj2, tm2, -(tm1+tm2),
     &        cf(2))
      endif
      

c     ==== Begin the recurrence ====
      do i = 3, idnint((j1+j2) - j_min + 1)
         j3 = j_min + i - 1
         call wfat_3j_a(j1, j2, j3, -m1-m2, a1)
         call wfat_3j_a(j1, j2, j3-1, -m1-m2, a2)
         call wfat_3j_b(j1, m1, j2, m2, j3-1, -m1-m2, b)
         cf(i) = (-b*cf(i-1) - j3*a2*cf(i-2)) / (j3-1) / a1
      enddo
       
      end
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_jscan000(j1, j2, cf)
c     DESCRIPTION:
c     Calculates the Wigner 3-j coefficients
c       / j1  j2  j3 \
c       \  0   0   0 /
c     for $j3 = [|j1-j2|, j1+j2]$.
      implicit none
      integer j1
      integer j2
      double precision cf(1:*)

      integer i, n, j3
      double precision  c
      

      n = 2*min(j1,j2) + 1
      cf(1:n) = 0.0d0
      
c     ==== Determine the type of the first 3-j symbol ====
c     / j1  j2   |j1-j2| \
c     \  0   0      0    /
      call wfat_3j_special_type5(j1, j2, cf(1))
     

c     ==== Begin the recurrence ====
      do i = 3, n, 2
         j3 = abs(j1-j2) + i - 1
         call wfat_3j_c(j1, j2, j3, c)
         cf(i) = -c * cf(i-2)
      enddo

      end
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_special_type1(tj1, tm1, tj2, tm2, cf)
c     Calculates the 3-j coefficients:
c     / j1  j2   |j1-j2| \
c     \ m1  m2  -(m1+m2) /
      implicit none 
#include "errquit.fh"
      integer tj1    ! [input], 2*j1
      integer tm1    ! [input], 2*m1
      integer tj2    ! [input], 2*j2
      integer tm2    ! [input], 2*m2
      double precision cf     ! [output], the 3-j coefficients

      integer i, n, j
      double precision f, j1, m1, j2, m2, j3, v(9), j1_, m1_, j2_, m2_


      if (mod(tj1,2) /= mod(abs(tm1),2))
     &     call errquit('WFAT: Both the arguments tj1 and tm1 of '//
     &     'wfat_3j_special_type1 must either be odd or even.', 6, UERR)
      if (mod(tj2,2) /= mod(abs(tm2),2))
     &     call errquit('WFAT: Both the arguments tj2 and tm2 of '//
     &     'wfat_3j_special_type1 must either be odd or even.', 6, UERR)

      
      if (abs(tj1-tj2) < abs(tm1+tm2)) then
         cf = 0.0d0
         return
      endif
      

c     Calculate
c        <(j1,-m1)(j2,-m2)|(|j1-j2|,-m1-m2)>
c     if j1>=j2, and
c        <(j2,-m2)(j1,-m1)|(|j1-j2|,-m1-m2)>
c     if j1<j2.
      if (tj1 >= tj2) then
         j1 = tj1 / 2.0d0
         m1 = -tm1 / 2.0d0
         j2 = tj2 / 2.0d0
         m2 = -tm2 / 2.0d0
      else
         j1 = tj2 / 2.0d0
         m1 = -tm2 / 2.0d0
         j2 = tj1 / 2.0d0
         m2 = -tm1 / 2.0d0
      endif
      j3 = abs(j1-j2)
      
      v(1) = j1 + m1
      v(2) = j1 - m1
      v(3) = 2*j2
      v(4) = 2*j1 - 2*j2 + 1
      v(5) = 2*j1 + 1
      v(6) = j2 + m2
      v(7) = j2 - m2
      v(8) = j1 - j2 + m1 + m2
      v(9) = j1 - j2 - m1 - m2
      do i = 1, 9
         if (abs(v(i)) < 1.0d-12) v(i) = 1.0d0
      enddo

      
c     ==== Sec. 8.5, eq. (13) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8)
     &        / v(9)          
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**idnint(j2+m2) * sqrt(f)


c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      j1_ = tj1 / 2.0d0
      m1_ = tm1 / 2.0d0
      j2_ = tj2 / 2.0d0
      m2_ = tm2 / 2.0d0
      if (tj1 >= tj2) then
         cf = (-1)**idnint(j3-m1_-m2_+tj1) / sqrt(2*j3+1) * cf
      else
c        == Take care of the interchange between j1 and j2 ==
c        ==       using Sec. 8.4, eq. (5) in Ref. [1]      ==
         cf = (-1)**idnint(j1_+3*j2_+2*j3-m1_-m2_) / sqrt(2*j3+1) * cf
      endif
      
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================

      
!==================================================
!==================================================
      subroutine wfat_3j_special_type2(tj1, tm1, tj2, tm2, tm3, cf)
c     Calculates the 3-j coefficients:
c     / j1  j2  |m3| \  or  / j1  j2  |m3| \
c     \ m1  m2   m3  /      \ m1  m2   -m3 /
      implicit none 
#include "errquit.fh"
      integer tj1    ! [input], 2*j1
      integer tm1    ! [input], 2*m1
      integer tj2    ! [input], 2*j2
      integer tm2    ! [input], 2*m2
      integer tm3    ! [input], 2*m3
      double precision cf     ! [output], the 3-j coefficients

      integer i, n, j
      double precision f, j1, m1, j2, m2, j3, v(9)


      if (mod(tj1,2) /= mod(abs(tm1),2))
     &     call errquit('WFAT: Both the arguments tj1 and tm1 of '//
     &     'wfat_3j_special_type2 must either be odd or even.', 6, UERR)
      if (mod(tj2,2) /= mod(abs(tm2),2))
     &     call errquit('WFAT: Both the arguments tj2 and tm2 of '//
     &     'wfat_3j_special_type2 must either be odd or even.', 6, UERR)
      

      if ( tm1+tm2+tm3 /= 0 .or. (abs(tm3) < abs(tj1-tj2) .or.
     &     abs(tm3) > tj1+tj2) ) then
         cf = 0.0d0
         return
      endif


c     Calculate
c        <(j1,-m1)(j2,-m2)|(|m3|,m3)>
c     if m3>=0, and
c        <(j1,m1)(j2,m2)|(|m3|,-m3)>
c     if m3<0.
      j1 = tj1 / 2.0d0
      if (tm3 >= 0) then
         m1 = -tm1 / 2.0d0
      else
         m1 = tm1 / 2.0d0
      endif
      j2 = tj2 / 2.0d0
      if (tm3 >= 0) then
         m2 = -tm2 / 2.0d0
      else
         m2 = tm2 / 2.0d0
      endif
      j3 = abs(tm3) / 2.0d0

      
      v(1) = 2*j3 + 1
      v(2) = j1 + j2 - j3
      v(3) = j1 + m1
      v(4) = j2 + m2
      v(5) = j1 + j2 + j3 + 1
      v(6) = j1 - j2 + j3
      v(7) = -j1 + j2 + j3
      v(8) = j1 - m1
      v(9) = j2 - m2
      do i = 1, 9
         if (abs(v(i)) < 1.0d-12) v(i) = 1.0d0
      enddo

      
c     ==== Sec. 8.5, eq. (35) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8)
     &        / v(9)          
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**idnint(j1-m1) * sqrt(f)

      
c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      if (tm3 >= 0) then
         cf = (-1)**idnint(2*j3+tj1) / sqrt(2*j3+1) * cf
      else
c        == Take care of the sign swap of m1 and m2 ==
c        ==   using Sec. 8.4, eq. (11) in Ref. [1]  ==
         cf = (-1)**idnint(3*j1+j2+tm3/2.0d0) / sqrt(2*j3+1) * cf
      endif
      
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_special_type3(tj1, tm1, tj2, tm2, cf)
c     / j1  j2  |j1-j2|+1 \
c     \ m1  m2   -(m1+m2) /
      implicit none
#include "errquit.fh"
      integer tj1    ! [input], 2*j1
      integer tm1    ! [input], 2*m1
      integer tj2    ! [input], 2*j2
      integer tm2    ! [input], 2*m2
      double precision cf     ! [output], the 3-j coefficients

      integer i, n, j
      double precision f, j1, m1, j2, m2, j3, v0, v(9), j1_, m1_, j2_
      double precision m2_
      

      if (mod(tj1,2) /= mod(abs(tm1),2))
     &     call errquit('WFAT: Both the arguments tj1 and tm1 of '//
     &     'wfat_3j_special_type3 must either be odd or even.', 6, UERR)
      if (mod(tj2,2) /= mod(abs(tm2),2))
     &     call errquit('WFAT: Both the arguments tj2 and tm2 of '//
     &     'wfat_3j_special_type3 must either be odd or even.', 6, UERR)

      
      if ((tm1 == 0 .and. tm2 == 0) .or. 
     &     (abs(tj1-tj2)+2 < abs(tm1+tm2))) then
         cf = 0.0d0
         return
      endif
      

c     Calculate
c        <(j1,-m1)(j2,-m2)|(|j1-j2|+1,-m1-m2)>
c     if j1>=j2, and
c        <(j2,-m2)(j1,-m1)|(|j1-j2|+1,-m1-m2)>
c     if j1<j2.
      if (tj1 >= tj2) then
         j1 = tj1 / 2.0d0
         m1 = -tm1 / 2.0d0
         j2 = tj2 / 2.0d0
         m2 = -tm2 / 2.0d0
      else
         j1 = tj2 / 2.0d0
         m1 = -tm2 / 2.0d0
         j2 = tj1 / 2.0d0
         m2 = -tm1 / 2.0d0
      endif
      j3 = abs(j1-j2)+1

      v0 = 2*j1 - 2*j2 + 3
      v(1) = 2*j2 - 1             
      v(2) = 2*j1 - 2*j2 + 1      
      v(3) = j1 + m1              
      v(4) = j1 - m1              
      v(5) = 2*j1 + 2             
      v(6) = j2 + m2              
      v(7) = j2 - m2              
      v(8) = j1 - j2 + m1 + m2 + 1
      v(9) = j1 - j2 - m1 - m2 + 1
      do i = 1, 9
         if (abs(v(i)) < 1.0d-12) v(i) = 1.0d0
      enddo

      
c     ==== Sec. 8.5, eq. (17) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8) 
     &        / v(9)
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**idnint(j2+m2+1) * 2*(j1*m2+j2*m1+m2) * sqrt(v0*f)


c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      j1_ = tj1 / 2.0d0
      m1_ = tm1 / 2.0d0
      j2_ = tj2 / 2.0d0
      m2_ = tm2 / 2.0d0
      if (tj1 >= tj2) then
         cf = (-1)**idnint(j3-m1_-m2_+tj1) / sqrt(2*j3+1) * cf
      else
c        == Take care of the interchange between j1 and j2 ==
c        ==       using Sec. 8.4, eq. (5) in Ref. [1]      ==
         cf = (-1)**idnint(j1_+3*j2_+2*j3-m1_-m2_) / sqrt(2*j3+1) * cf
      endif
      
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================
      
      
!==================================================
!==================================================
      subroutine wfat_3j_special_type4(tj1, tm1, tj2, tm2, tm3, cf)
c     Calculates the 3-j coefficients:
c     / j1  j2  |m3|+1 \  or  / j1  j2  |m3|+1 \
c     \ m1  m2    m3   /      \ m1  m2    -m3  /
      implicit none
#include "errquit.fh"
      integer tj1    ! [input], 2*j1
      integer tm1    ! [input], 2*m1
      integer tj2    ! [input], 2*j2
      integer tm2    ! [input], 2*m2
      integer tm3    ! [input], 2*m3
      double precision cf     ! [output], the 3-j coefficients

      integer i, n, j
      double precision f, j1, m1, j2, m2, j3, c, v0, v(9)


      if (mod(tj1,2) /= mod(abs(tm1),2))
     &     call errquit('WFAT: Both the arguments tj1 and tm1 of '//
     &     'wfat_3j_special_type4 must either be odd or even.', 6, UERR)
      if (mod(tj2,2) /= mod(abs(tm2),2))
     &     call errquit('WFAT: Both the arguments tj2 and tm2 of '//
     &     'wfat_3j_special_type4 must either be odd or even.', 6, UERR)
      

      if ( tm1+tm2+tm3 /= 0 .or. (abs(tm3)+2 < abs(tj1-tj2) .or.
     &     abs(tm3)+2 > tj1+tj2) ) then
         cf = 0.0d0
         return
      endif


c     Calculate
c        <(j1,-m1)(j2,-m2)|(|m3|+1,m3)>
c     if m3>=0, and
c        <(j1,m1)(j2,m2)|(|m3|+1,-m3)>
c     if m3<0.
      j1 = tj1 / 2.0d0
      if (tm3 >= 0) then
         m1 = -tm1 / 2.0d0
      else
         m1 = tm1 / 2.0d0
      endif
      j2 = tj2 / 2.0d0
      if (tm3 >= 0) then
         m2 = -tm2 / 2.0d0
      else
         m2 = tm2 / 2.0d0
      endif
      j3 = abs(tm3) / 2.0d0 + 1
      c = j3 - 1

      v0 = 2*c + 3
      v(1) = 2*c + 1
      v(2) = j1 + j2 - c - 1
      v(3) = j1 + m1
      v(4) = j2 + m2
      v(5) = j1 + j2 + c + 2
      v(6) = j1 - j2 + c + 1
      v(7) = -j1 + j2 + c + 1
      v(8) = j1 - m1
      v(9) = j2 - m2
      do i = 1, 9
         if (abs(v(i)) < 1.0d-12) v(i) = 1.0d0
      enddo

      
c     ==== Sec. 8.5, eq. (44) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8)
     &        / v(9)          
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**idnint(j1-m1) * ((j2-m2)*(j2+m2+1) - (j1-m1)*(j1+m1+1))
     &     * sqrt(v0*f)


c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      if (tm3 >= 0) then
         cf = (-1)**idnint(tm3+1.0d0+tj1) / sqrt(2*j3+1) * cf
      else
c        == Take care of the sign swap of m1 and m2 ==
c        ==   using Sec. 8.4, eq. (11) in Ref. [1]  ==
         cf = (-1)**idnint(3*j1+j2+tm3/2.0d0) / sqrt(2*j3+1) * cf
      endif
      
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_special_type5(j1, j2, cf)
c     Calculates the 3-j coefficients:
c     / j1  j2   |j1-j2| \
c     \  0   0      0    /
      implicit none
      integer j1
      integer j2
      double precision cf

      integer i, n, j, j1_, j2_, j3, phase_def
      double precision f, v(9)


c     Calculate
c        <(j1,0)(j2,0)|(|j1-j2|,0)>
c     if j1>=j2, and
c        <(j2,0)(j1,0)|(|j1-j2|,0)>
c     if j1<j2.
      j3 = abs(j1-j2)
      phase_def = (-1)**(j3+2*j1)    ! The phase from the definition of 3j in terms of CG.
      if (j1 >= j2) then
         j1_ = j1
         j2_ = j2
      else
         j1_ = j2
         j2_ = j1
      endif


      v(1) = dble(j1_)
      v(2) = dble(j1_)
      v(3) = dble(2*j2_)
      v(4) = dble(2*j1_ - 2*j2_ + 1)
      v(5) = dble(j2_)
      v(6) = dble(j2_)
      v(7) = dble(j1_ - j2_)
      v(8) = dble(j1_ - j2_)
      v(9) = dble(2*j1_ + 1)
      do i = 1, 9
         if (abs(v(i)) < 1.0d-12) v(i) = 1.0d0
      enddo


c     ==== Sec. 8.5, eq. (34) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8)
     &        / v(9)          
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**j2_ * sqrt(f)
      

c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      if (j1 >= j2) then
         cf = phase_def / sqrt(2.0d0*j3+1) * cf
      else
c        == Take care of the interchange between j1 and j2 ==
c        ==       using Sec. 8.4, eq. (5) in Ref. [1]      ==
         cf = (-1)**(j1+j2+j3) * phase_def / sqrt(2.0d0*j3+1) * cf
      endif
         
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_a(j1, j2, j3, m3, a)

      implicit none
      double precision j1
      double precision j2
      double precision j3
      double precision m3
      double precision a

      a = sqrt(j3**2-(j1-j2)**2) * sqrt((j1+j2+1)**2-j3**2) *
     &     sqrt(j3**2-m3**2)

      end
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_b(j1, m1, j2, m2, j3, m3, b)

      implicit none
      double precision j1
      double precision m1
      double precision j2
      double precision m2
      double precision j3
      double precision m3
      double precision b

      b = -(2*j3+1) * (j1*(j1+1)*m3 - j2*(j2+1)*m3 - j3*(j3+1)*(m2-m1))
      
      end
!==================================================
      
      
!==================================================
!==================================================
      subroutine wfat_3j_c(j1, j2, j3, c)

      implicit none
      integer j1
      integer j2
      integer j3
      double precision c
      
      double precision c1, c2, c3, c4

      c1 = (j3-1)**2 - (j1-j2)**2
      c2 = (j1+j2+1)**2 - (j3-1)**2
      c3 = j3**2 - (j1-j2)**2
      c4 = (j1+j2+1)**2 - j3**2
      c = sqrt(c1*c2/c3/c4)

      end
!==================================================
      
