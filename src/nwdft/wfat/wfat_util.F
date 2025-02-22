c=======================================================
c=======================================================
      subroutine wfat_ival_atleast(val, minim, catgr, valname,
     &           w, r)

      implicit none
#include "errquit.fh"
      integer val
      integer minim
      character*(*) catgr
      character*(*) valname
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.

      character*20 minim_c


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_ival_atleast'' is unavailable.',
     &        0, UERR)
      endif
     

      write(minim_c, '(i10)') minim
      if (val < minim) then
         r = .false.
         if (w == 'i') then
            call errquit(trim(adjustl(catgr))//': '//
     &           'The variable named '//trim(adjustl(valname))//
     &           ' must be greater than or equal to '//
     &           trim(adjustl(minim_c))//'.', 0, UERR)
         else if (w == 'o') then
            return
         endif
      else
         r = .true.
         return
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_ival_atmost(val, maxim, catgr, valname,
     &           w, r)

      implicit none
#include "errquit.fh"
      integer val
      integer maxim
      character*(*) catgr
      character*(*) valname
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.

      character*20 maxim_c


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_ival_atmost'' is unavailable.',
     &        0, UERR)
      endif
      

      write(maxim_c, '(i10)') maxim
      if (val > maxim) then
         r = .false.
         if (w == 'i') then
            call errquit(trim(adjustl(catgr))//': '//
     &           'The variable named '//trim(adjustl(valname))//
     &           ' must be less than or equal to '//
     &           trim(adjustl(maxim_c))//'.', 0, UERR)
         else if (w == 'o') then
            return
         endif
      else
         r = .true.
         return
      endif

      end
c=======================================================      


c=======================================================
c=======================================================
      subroutine wfat_ival_inrange(val, minim, maxim, catgr, valname,
     &           w, r)

      implicit none
#include "errquit.fh"
      integer val
      integer minim
      integer maxim
      character*(*) catgr
      character*(*) valname
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.

      character*20 minim_c, maxim_c


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_ival_inrange'' is unavailable.',
     &        0, UERR)
      endif
      

      write(minim_c, '(i10)') minim
      write(maxim_c, '(i10)') maxim
      
      if (val < minim .or. val > maxim) then
         r = .false.
         if (w == 'i') then
            call errquit(trim(adjustl(catgr))//': '//
     &           'The variable named '//trim(adjustl(valname))//
     &           ' must be in between '//
     &           trim(adjustl(minim_c))//' and '//
     &           trim(adjustl(maxim_c))//'.', 0, UERR)
         else if (w == 'o') then
            return
         endif
      else
         r = .true.
         return
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_ival_outrange(val, minim, maxim, catgr, valname,
     &           w, r)

      implicit none
#include "errquit.fh"
      integer val
      integer minim
      integer maxim
      character*(*) catgr
      character*(*) valname
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.

      character*20 minim_c, maxim_c


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_ival_outrange'' is unavailable.',
     &        0, UERR)
      endif
      

      write(minim_c, '(i10)') minim
      write(maxim_c, '(i10)') maxim
      
      if (val > minim .and. val < maxim) then
         r = .false.
         if (w == 'i') then
            call errquit(trim(adjustl(catgr))//': '//
     &           'The variable named '//trim(adjustl(valname))//
     &           ' must lie outside of the range bounded by '//
     &           trim(adjustl(minim_c))//' and '//
     &           trim(adjustl(maxim_c))//'.', 0, UERR)
         else if (w == 'o') then
            return
         endif
      else
         r = .true.
         return
      endif

      end
c=======================================================      


c=======================================================
c=======================================================
      subroutine wfat_ivec_atleast(nvec, vec, minim, catgr, vecname,
     &           w, r, k)

      implicit none
#include "errquit.fh"
      integer nvec             ! [input] The number of elements of vec.
      integer vec(nvec)        ! [input] The vector to be tested.
      integer minim            ! [input] The smallest value allowed.
      character*(*) catgr      ! [input] Error category for error printing purpose, only relevant when w='i'.
      character*(*) vecname    ! [input] The name of the vector variable, only relevant when w='i'.
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.
      integer k                ! [output] In the case of w='o', returns -1 if requirement is satisfied,
                               !          otherwise k contains the index of the first element not satisfying 
                               !          the requirement. When w='i', returns -1 if requirement is satisfied,
                               !          irrelevant otherwise since the execution is terminated inside this 
                               !          subroutine.

      character*7 i_c
      character*3 suf
      character*20 minim_c
      logical r0
      integer i


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_ivec_atleast'' is unavailable.',
     &        0, UERR)
      endif
     

      write(minim_c, '(i10)') minim
      
      do i = 1, nvec
         call wfat_ival_atleast(vec(i), minim, catgr, vecname, 'o', r0)

         if (.not. r0) then
            r = .false.
            k = i
            if (w == 'i') then
               if (mod(i-1,10) == 0) then
                  suf = '-st'
               else if (mod(i-2,10) == 0) then
                  suf = '-nd'
               else if (mod(i-3,10) == 0) then
                  suf = '-rd'
               else
                  suf = '-th'
               endif
               write(i_c, '(i7)') i
               call errquit(trim(adjustl(catgr))//': '//
     &              'The elements of a vector variable named '//
     &              trim(adjustl(vecname))//' must be greater than '//
     &              'or equal to '//trim(adjustl(minim_c))//'. '//
     &              'The '//trim(adjustl(i_c))//suf//' element '//
     &              'does not satisfy this.', 0, UERR)

            else if (w == 'o') then
               return
            endif
         endif

      enddo

      k = -1
      r = .true.

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_ivec_atmost(nvec, vec, maxim, catgr, vecname,
     &           w, r, k)

      implicit none
#include "errquit.fh"
      integer nvec             ! [input] The number of elements of vec.
      integer vec(nvec)        ! [input] The vector to be tested.
      integer maxim            ! [input] The smallest value allowed.
      character*(*) catgr      ! [input] Error category for error printing purpose, only relevant when w='i'.
      character*(*) vecname    ! [input] The name of the vector variable, only relevant when w='i'.
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.
      integer k                ! [output] In the case of w='o', returns -1 if requirement is satisfied,
                               !          otherwise k contains the index of the first element not satisfying 
                               !          the requirement. When w='i', returns -1 if requirement is satisfied,
                               !          irrelevant otherwise since the execution is terminated inside this 
                               !          subroutine.

      character*7 i_c
      character*3 suf
      character*20 maxim_c
      logical r0
      integer i


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_ivec_atmost'' is unavailable.',
     &        0, UERR)
      endif
     

      write(maxim_c, '(i10)') maxim
      
      do i = 1, nvec
         call wfat_ival_atmost(vec(i), maxim, catgr, vecname, 'o', r0)

         if (.not. r0) then
            r = .false.
            k = i
            if (w == 'i') then
               if (mod(i-1,10) == 0) then
                  suf = '-st'
               else if (mod(i-2,10) == 0) then
                  suf = '-nd'
               else if (mod(i-3,10) == 0) then
                  suf = '-rd'
               else
                  suf = '-th'
               endif
               write(i_c, '(i7)') i
               call errquit(trim(adjustl(catgr))//': '//
     &              'The elements of a vector variable named '//
     &              trim(adjustl(vecname))//' must be less than '//
     &              'or equal to '//trim(adjustl(maxim_c))//'. '//
     &              'The '//trim(adjustl(i_c))//suf//' element '//
     &              'does not satisfy this.', 0, UERR)

            else if (w == 'o') then
               return
            endif
         endif

      enddo

      k = -1
      r = .true.

      end
c=======================================================      


c=======================================================
c=======================================================
      subroutine wfat_ivec_inrange(nvec, vec, minim, maxim, catgr,
     &           vecname, w, r, k)

      implicit none
#include "errquit.fh"
      integer nvec             ! [input] The number of elements of vec.
      integer vec(nvec)        ! [input] The vector to be tested.
      integer minim            ! [input] The smallest value allowed.
      integer maxim            ! [input] The largest value allowed.
      character*(*) catgr      ! [input] Error category for error printing purpose, only relevant when w='i'.
      character*(*) vecname    ! [input] The name of the vector variable, only relevant when w='i'.
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.
      integer k                ! [output] In the case of w='o', returns -1 if requirement is satisfied,
                               !          otherwise k contains the index of the first element not satisfying 
                               !          the requirement. When w='i', returns -1 if requirement is satisfied,
                               !          irrelevant otherwise since the execution is terminated inside this 
                               !          subroutine.

      character*7 i_c
      character*3 suf
      character*20 minim_c, maxim_c
      logical r0
      integer i


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_ivec_inrange'' is unavailable.',
     &        0, UERR)
      endif
     

      write(minim_c, '(i10)') minim
      write(maxim_c, '(i10)') maxim
      
      do i = 1, nvec
         call wfat_ival_inrange(vec(i), minim, maxim, catgr, vecname,
     &                          'o', r0)

         if (.not. r0) then
            r = .false.
            k = i
            if (w == 'i') then
               if (mod(i-1,10) == 0) then
                  suf = '-st'
               else if (mod(i-2,10) == 0) then
                  suf = '-nd'
               else if (mod(i-3,10) == 0) then
                  suf = '-rd'
               else
                  suf = '-th'
               endif
               write(i_c, '(i7)') i
               call errquit(trim(adjustl(catgr))//': '//
     &              'The elements of a vector variable named '//
     &              trim(adjustl(vecname))//
     &              ' must be greater than or equal to '//
     &              trim(adjustl(minim_c))//
     &              ' AND less than or equal to '//
     &              trim(adjustl(maxim_c))//
     &              '. The '//trim(adjustl(i_c))//suf//' element '//
     &              'does not satisfy this.', 0, UERR)

            else if (w == 'o') then
               return
            endif
         endif

      enddo

      k = -1
      r = .true.

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_ivec_outrange(nvec, vec, minim, maxim, catgr,
     &           vecname, w, r, k)

      implicit none
#include "errquit.fh"
      integer nvec             ! [input] The number of elements of vec.
      integer vec(nvec)        ! [input] The vector to be tested.
      integer minim            ! [input] The smallest value allowed.
      integer maxim            ! [input] The largest value allowed.
      character*(*) catgr      ! [input] Error category for error printing purpose, only relevant when w='i'.
      character*(*) vecname    ! [input] The name of the vector variable, only relevant when w='i'.
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.
      integer k                ! [output] In the case of w='o', returns -1 if requirement is satisfied,
                               !          otherwise k contains the index of the first element not satisfying 
                               !          the requirement. When w='i', returns -1 if requirement is satisfied,
                               !          irrelevant otherwise since the execution is terminated inside this 
                               !          subroutine.

      character*7 i_c
      character*3 suf
      character*20 minim_c, maxim_c
      logical r0
      integer i


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_ivec_outrange'' is unavailable.',
     &        0, UERR)
      endif
     

      write(minim_c, '(i10)') minim
      write(maxim_c, '(i10)') maxim
      
      do i = 1, nvec
         call wfat_ival_outrange(vec(i), minim, maxim, catgr, vecname,
     &                          'o', r0)

         if (.not. r0) then
            r = .false.
            k = i
            if (w == 'i') then
               if (mod(i-1,10) == 0) then
                  suf = '-st'
               else if (mod(i-2,10) == 0) then
                  suf = '-nd'
               else if (mod(i-3,10) == 0) then
                  suf = '-rd'
               else
                  suf = '-th'
               endif
               write(i_c, '(i7)') i
               call errquit(trim(adjustl(catgr))//': '//
     &              'The elements of a vector variable named '//
     &              trim(adjustl(vecname))//
     &              ' must be less than or equal to '//
     &              trim(adjustl(minim_c))//
     &              ' OR greater than or equal to '//
     &              trim(adjustl(maxim_c))//
     &              '. The '//trim(adjustl(i_c))//suf//' element '//
     &              'does not satisfy this.', 0, UERR)

            else if (w == 'o') then
               return
            endif
         endif

      enddo

      k = -1
      r = .true.

      end
c=======================================================      

      
c=======================================================
c=======================================================
      subroutine wfat_dval_atleast(val, minim, e, catgr, valname,
     &           style, w, r)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      double precision val
      double precision minim
      character*1 e
      character*(*) catgr
      character*(*) valname
      character*1 style
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.

      character*20 ss
      character*20 minim_c
      logical e_q, e_g


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_dval_atleast'' is unavailable.',
     &        0, UERR)
      endif
      

      if (style == 'e') then
         write(minim_c, '(es15.6)') minim
      else if (style == 'f') then
         write(minim_c, '(f11.6)') minim
      else
         call errquit('The value of the argument ''style'' of the '//
     &        'subroutine ''wfat_dval_atleast'' is unavailable.',
     &        0, UERR)
      endif

      e_q = .true.
      e_g = .true.
      
c     ==== Greater than or equal to ====
      if (e == 'q') then
         ss = ' or equal to '
         if (val < minim) e_q = .false.
c     ==== Greater than ====         
      else if (e == 'g') then
         ss = ' '
         if (val <= minim) e_g = .false.
      else
         call errquit('The value of the argument ''e'' of the '//
     &        'subroutine ''wfat_dval_atleast'' is unavailable.',
     &        0, UERR)
      endif


      if ((.not. e_q) .or. (.not. e_g)) then
         r = .false.
         if (w == 'i') then
            write(luout, '("value of ", a, " = ", es14.6)') valname, val
            call errquit(trim(adjustl(catgr))//': '//
     &                   'The variable named '//trim(adjustl(valname))//
     &                   ' must be greater than '//trim(adjustl(ss))//
     &                   ' '//trim(adjustl(minim_c))//'.', 0, UERR)
         else if (w == 'o') then
            return
         endif
      else
         r = .true.
         return
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dval_atmost(val, maxim, e, catgr, valname,
     &           style, w, r)

      implicit none
#include "errquit.fh"
      double precision val
      double precision maxim
      character*1 e
      character*(*) catgr
      character*(*) valname
      character*1 style
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.

      character*20 ss
      character*20 maxim_c
      logical e_q, e_l


      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_dval_atmost'' is unavailable.',
     &        0, UERR)
      endif
      

      if (style == 'e') then
         write(maxim_c, '(es15.6)') maxim
      else if (style == 'f') then
         write(maxim_c, '(f11.6)') maxim
      else
         call errquit('The value of the argument ''style'' of the '//
     &        'subroutine ''wfat_dval_atmost'' is unavailable.',
     &        0, UERR)
      endif

      e_q = .true.
      e_l = .true.

c     ==== Less than or equal to ====
      if (e == 'q') then
         ss = ' or equal to '
         if (val > maxim) e_q = .false.
c     ==== Less than ====         
      else if (e == 'l') then
         ss = ' '
         if (val >= maxim) e_l = .false.
      else
         call errquit('The value of the argument ''e'' of the '//
     &        'subroutine ''wfat_dval_atmost'' is unavailable.',
     &        0, UERR)
      endif


      if ((.not. e_q) .or. (.not. e_l)) then
         r = .false.
         if (w == 'i') then
            call errquit(trim(adjustl(catgr))//': '//
     &                   'The variable named '//trim(adjustl(valname))//
     &                   ' must be less than '//trim(adjustl(ss))//' '//
     &                   trim(adjustl(maxim_c))//'.', 0, UERR)
         else if (w == 'o') then
            return
         endif
      else
         r = .true.
         return
      endif
            
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dval_inrange(val, minim, e1, maxim, e2,
     &     catgr, valname, style, w, r)

      implicit none
#include "errquit.fh"
      double precision val     ! [input]
      double precision minim   ! [input]
      character*1 e1           ! [input] Must be 'q' for '>=' or 'g' for '>'   
      double precision maxim   ! [input]
      character*1 e2           ! [input] Must be 'q' for '<=' or 'l' for '<'   
      character*(*) catgr      ! [input]
      character*(*) valname    ! [input]
      character*1 style        ! [input] Must be 'e' for exponential format or 'f' for decimal format.
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.

      character*20 ss1, ss2
      character*20 minim_c, maxim_c
      logical e1_q, e1_g, e2_q, e2_l

      
      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_dval_inrange'' is unavailable.',
     &        0, UERR)
      endif
      

      if (style == 'e') then
         write(minim_c, '(es15.6)') minim
         write(maxim_c, '(es15.6)') maxim
      else if (style == 'f') then
         write(minim_c, '(f11.6)') minim
         write(maxim_c, '(f11.6)') maxim
      else
         call errquit('The value of the argument ''style'' of the '//
     &        'subroutine ''wfat_dval_inrange'' is unavailable.',
     &        0, UERR)
      endif

      e1_q = .true.
      e1_g = .true.
      e2_q = .true.
      e2_l = .true.

c     ==== Greater than or equal to minim ====
      if (e1 == 'q') then
         ss1 = ' or equal to '
         if (val < minim) e1_q = .false.
c     ==== Greater than minim ====         
      else if (e1 == 'g') then
         ss1 = ' '
         if (val <= minim) e1_g = .false.
      else
         call errquit('The value of the argument ''e1'' of the '//
     &        'subroutine ''wfat_dval_inrange'' is unavailable.',
     &        0, UERR)
      endif

c     ==== Less than or equal to maxim ====
      if (e2 == 'q') then
         ss2 = ' or equal to '
         if (val > maxim) e2_q = .false.
c     ==== Less than maxim ====
      else if (e2 == 'l') then
         ss2 = ' '
         if (val >= maxim) e2_l = .false.
      else
         call errquit('The value of the argument ''e2'' of the '//
     &        'subroutine ''wfat_dval_inrange'' is unavailable.',
     &        0, UERR)
      endif


      if ( .not. (e1_q .and. e2_q)
     &     .or.
     &     .not. (e1_q .and. e2_l)
     &     .or.
     &     .not. (e1_g .and. e2_q)
     &     .or.
     &     .not. (e1_g .and. e2_l)) then
         r = .false.
         if (w == 'i') then
            call errquit(trim(adjustl(catgr))//': '//
     &                   'The variable named '//trim(adjustl(valname))//
     &                   ' must be greater than '//trim(adjustl(ss1))//
     &                   ' '//trim(adjustl(minim_c))//
     &                   ' AND less than '//trim(adjustl(ss2))//
     &                   trim(adjustl(maxim_c))//'.', 0, UERR)
         else if (w == 'o') then
            return
         endif
      else
         r = .true.
         return
      endif
      
      
      end
c=======================================================            


c=======================================================
c=======================================================
      subroutine wfat_dval_outrange(val, minim, e1, maxim, e2,
     &     catgr, valname, style, w, r)

      implicit none
#include "errquit.fh"
      double precision val     ! [input]
      double precision minim   ! [input]
      character*1 e1           ! [input] Must be 'q' for '>=' or 'g' for '>'   
      double precision maxim   ! [input]
      character*1 e2           ! [input] Must be 'q' for '<=' or 'l' for '<'   
      character*(*) catgr      ! [input]
      character*(*) valname    ! [input]
      character*1 style        ! [input] Must be 'e' for exponential format or 'f' for decimal format.
      character*1 w            ! [input] Must be 'i' if program termination (in case requirement is not
                               !         satisfied) is desired inside this routine or 'o' if an error message
                               !         (contained in argument 'r') is desired to be catched by the
                               !         calling program unit.
      logical r                ! [output] In the case of w='o', returns .true. if requirement is satisfied,
                               !          .false. otherwise. When w='i', returns .true. if requirement is
                               !          satisfied, irrelevant otherwise since the execution is terminated
                               !          inside this subroutine.

      character*20 ss1, ss2
      character*20 minim_c, maxim_c
      logical e1_q, e1_g, e2_q, e2_l

      
      if (w /= 'i' .and. w /= 'o') then
         call errquit('The value of the argument ''w'' of the '//
     &        'subroutine ''wfat_dval_outrange'' is unavailable.',
     &        0, UERR)
      endif
      

      if (style == 'e') then
         write(minim_c, '(es15.6)') minim
         write(maxim_c, '(es15.6)') maxim
      else if (style == 'f') then
         write(minim_c, '(f11.6)') minim
         write(maxim_c, '(f11.6)') maxim
      else
         call errquit('The value of the argument ''style'' of the '//
     &        'subroutine ''wfat_dval_outrange'' is unavailable.',
     &        0, UERR)
      endif

      e1_q = .true.
      e1_g = .true.
      e2_q = .true.
      e2_l = .true.

c     ==== Less than or equal to minim ====
      if (e1 == 'q') then
         ss1 = ' or equal to '
         if (val > minim) e1_q = .false.
c     ==== Less than minim ====         
      else if (e1 == 'g') then
         ss1 = ' '
         if (val >= minim) e1_g = .false.
      else
         call errquit('The value of the argument ''e1'' of the '//
     &        'subroutine ''wfat_dval_outrange'' is unavailable.',
     &        0, UERR)
      endif

c     ==== Greater than or equal to maxim ====
      if (e2 == 'q') then
         ss2 = ' or equal to '
         if (val < maxim) e2_q = .false.
c     ==== Greater than maxim ====
      else if (e2 == 'l') then
         ss2 = ' '
         if (val <= maxim) e2_l = .false.
      else
         call errquit('The value of the argument ''e2'' of the '//
     &        'subroutine ''wfat_dval_outrange'' is unavailable.',
     &        0, UERR)
      endif


      if ( .not. (e1_q .or. e2_q)
     &     .or.
     &     .not. (e1_q .or. e2_l)
     &     .or.
     &     .not. (e1_g .or. e2_q)
     &     .or.
     &     .not. (e1_g .or. e2_l)) then
         r = .false.
         if (w == 'i') then
            call errquit(trim(adjustl(catgr))//': '//
     &                   'The variable named '//trim(adjustl(valname))//
     &                   ' must be less than '//trim(adjustl(ss1))//
     &                   ' '//trim(adjustl(minim_c))//
     &                   ' OR greater than '//trim(adjustl(ss2))//
     &                   trim(adjustl(maxim_c))//'.', 0, UERR)
         else if (w == 'o') then
            return
         endif
      else
         r = .true.
         return
      endif
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_check_invec(k, n, vec)

      implicit none
      integer k
      integer n
      integer vec(n)

      integer i
      

      wfat_check_invec = .false.
      do i = 1, n
         if (vec(i) == k) then
            wfat_check_invec = .true.
            exit
         endif
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_ordinal_suffix(i, s)

      implicit none
      integer i
      character*2 s

      integer j
      

      j = mod(i,10)

      if (j == 1) then
         s = 'st'
      else if (j == 2) then
         s = 'nd'
      else if (j == 3) then
         s = 'rd'
      else
         s = 'th'
      endif


      if ( mod(i-11, 100)==0 .or.
     &     mod(i-12, 100)==0 .or.
     &     mod(i-13, 100)==0) s = 'th'

      end
c=======================================================


c=======================================================
c=======================================================
      double precision function wfat_wvnm_afreqau(vin)

      implicit none
      double precision vin
      
      double precision pi_const, lspeed, nm2bohr
      parameter (pi_const=3.14159265359d0, lspeed=137.035999084d0)
      parameter (nm2bohr=18.897261255d0)      


      wfat_wvnm_afreqau = 2*pi_const*lspeed / (vin*nm2bohr)

      end
c=======================================================

      
c=======================================================
c=======================================================
      double precision function wfat_tfs_tau(t, vin)

      implicit none
#include "errquit.fh"      
      integer t
      double precision vin

      double precision fs2au
      parameter (fs2au=41.3413733352d0)


      if (t == 1) then
         wfat_tfs_tau = fs2au * vin
      else if (t == 2) then
         wfat_tfs_tau = vin / fs2au
      else
         call errquit('WFAT: Unknown value of argument ''t'' '//
     &        'of wfat_tfs_tau.', 1, UERR)
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      double precision function wfat_convert_freq(unit_in, val_in,
     &     unit_out)

      implicit none
#include "errquit.fh"      
      character*(*) unit_in
      double precision val_in
      character*(*) unit_out

      double precision dd, val_au, factor1, factor2
      double precision au2au, ev2au, thz2au, nm2au, bohr2au, au2ev
      double precision au2thz, au2nm, au2bohr

      double precision pi_const, lspeed, nm2bohr, fs2au
      parameter
     &     (pi_const=3.14159265359d0, lspeed=137.035999084d0,
     &     nm2bohr=18.897261255d0, fs2au=41.3413733352d0)
      
      parameter
     &     (au2au=1.0d0,
     &     ev2au=0.036749308136649d0,
     &     thz2au=2*pi_const*1.0d12 / (fs2au*1.0d15),       ! = 2*pi*2.4188843265857e-5
     &     nm2au=2*pi_const*lspeed / nm2bohr,
     &     bohr2au=2*pi_const*lspeed,
     &     au2ev=1.0d0/ev2au,
     &     au2thz=1.0d0/thz2au,
     &     au2nm=1.0d0/nm2au,
     &     au2bohr=1.0d0/bohr2au)


      if (unit_in/='au' .and. unit_in/='ev' .and. unit_in/='thz' .and.
     &     unit_in/='nm' .and. unit_in/='bohr')
     &     call errquit('wfat_convert_freq: Value of ''unit_in'' is '//
     &     'unavailable', 1, UERR)
      if (unit_out/='au' .and. unit_out/='ev' .and. unit_out/='thz' 
     &     .and. unit_out/='nm' .and. unit_out/='bohr')
     &     call errquit('wfat_convert_freq: Value of ''unit_out'' is '//
     &     'unavailable', 1, UERR)
      
      
      if (unit_in == 'nm' .or. unit_in == 'bohr') then
         dd = 1.0d0/val_in
      else
         dd = val_in
      endif
      

      if (unit_in == 'au') factor1 = au2au
      if (unit_in == 'ev') factor1 = ev2au
      if (unit_in == 'thz') factor1 = thz2au
      if (unit_in == 'nm') factor1 = nm2au
      if (unit_in == 'bohr') factor1 = bohr2au

      val_au = factor1 * dd          ! val_au is in atomic unit of angular frequency.

      if (unit_out == 'au') factor2 = au2au
      if (unit_out == 'ev') factor2 = au2ev
      if (unit_out == 'thz') factor2 = au2thz
      if (unit_out == 'nm') factor2 = au2nm
      if (unit_out == 'bohr') factor2 = au2bohr

      wfat_convert_freq = factor2 * val_au


      if (unit_out == 'nm' .or. unit_out == 'bohr') then
         wfat_convert_freq = 1.0d0/wfat_convert_freq
      else
         wfat_convert_freq = wfat_convert_freq
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      double precision function wfat_convert_time(unit_in, val_in,
     &     unit_out)

      implicit none
#include "errquit.fh"      
      character*(*) unit_in
      double precision val_in
      character*(*) unit_out

      double precision au2au, fs2au, au2fs, val_au, factor1, factor2
      parameter
     &     (au2au=1.0d0,
     &     fs2au=41.3413733352d0,
     &     au2fs=1.0d0/fs2au)


      if (unit_in/='au' .and. unit_in/='fs')
     &     call errquit('wfat_convert_time: Value of ''unit_in'' is '//
     &     'unavailable', 1, UERR)
      if (unit_out/='au' .and. unit_out/='fs')
     &     call errquit('wfat_convert_time: Value of ''unit_out'' is '//
     &     'unavailable', 1, UERR)
      

      if (unit_in == 'au') factor1 = au2au
      if (unit_in == 'fs') factor1 = fs2au

      val_au = factor1 * val_in         ! val_au is in atomic unit of time.

      if (unit_out == 'au') factor2 = au2au
      if (unit_out == 'fs') factor2 = au2fs

      wfat_convert_time = factor2 * val_au

      end
c=======================================================


c=======================================================
c=======================================================
      double precision function wfat_convert_length(unit_in, val_in,
     &     unit_out)

      implicit none
#include "errquit.fh"      
      character*(*) unit_in
      double precision val_in
      character*(*) unit_out

      double precision bo2bo, nm2bo, an2bo, bo2nm, bo2an
      double precision val_bo, factor1, factor2
      parameter
     &     (bo2bo=1.0d0,
     &     bo2an=0.529177210903d0,
     &     bo2nm=bo2an/10.0d0,
     &     nm2bo=1.0d0/bo2nm,
     &     an2bo=1.0d0/bo2an)


      if (unit_in/='bohr' .and. unit_in/='angstrom' .and. unit_in/='nm')
     &     call errquit('wfat_convert_length: Value of ''unit_in'' '//
     &     'is unavailable', 1, UERR)
      if (unit_out/='bohr' .and. unit_out/='angstrom' .and.
     &     unit_in/='nm')
     &     call errquit('wfat_convert_length: Value of ''unit_out'' '//
     &     'is unavailable', 1, UERR)
      

      if (unit_in == 'bohr') factor1 = bo2bo
      if (unit_in == 'angstrom') factor1 = an2bo
      if (unit_in == 'nm') factor1 = nm2bo

      val_bo = factor1 * val_in         ! val_bo is in bohr.

      if (unit_out == 'bohr') factor2 = bo2bo
      if (unit_out == 'angstrom') factor2 = bo2an
      if (unit_out == 'nm') factor2 = bo2nm

      wfat_convert_length = factor2 * val_bo

      end
c=======================================================


c=======================================================
c=======================================================
      double precision function wfat_convert_efield(unit_in, val_in,
     &     unit_out)

      implicit none
#include "errquit.fh"
      character*(*) unit_in
      double precision val_in
      character*(*) unit_out

      double precision au2twcm, ff
      parameter (au2twcm=3.5094448d4)

      if (unit_in == 'au' .and. unit_out == 'tw/cm2') then
         ff = au2twcm * val_in**2
      else if (unit_in == 'tw/cm2' .and. unit_out == 'au') then
         ff = sqrt(val_in/au2twcm)
      else if (trim(adjustl(unit_in)) == trim(adjustl(unit_out))) then
         ff = val_in
      else
         call errquit('wfat_convert_efield: Unit option is unavailable',
     &        1, UERR)
      endif

      wfat_convert_efield = ff
      
      end
c=======================================================


c=======================================================
c=======================================================
      double precision function wfat_convert_angle(unit_in, val_in,
     &     unit_out)

      implicit none
#include "errquit.fh"
      character*(*) unit_in
      double precision val_in
      character*(*) unit_out

      double precision factor1, factor2
      double precision rad2rad, deg2rad, rad2deg, val_rad
      parameter
     &     (rad2rad = 1.0d0,
     &     deg2rad=0.01745329251d0,
     &     rad2deg=1.0d0/deg2rad)

      
      if (unit_in/='rad' .and. unit_in/='deg')
     &     call errquit('wfat_convert_angle: Value of ''unit_in'' is '//
     &     'unavailable', 1, UERR)
      if (unit_out/='rad' .and. unit_out/='deg')
     &     call errquit('wfat_convert_angle: Value of ''unit_out'' is'//
     &     ' unavailable', 1, UERR)

      if (unit_in == 'rad') factor1 = rad2rad
      if (unit_in == 'deg') factor1 = deg2rad

      val_rad = factor1 * val_in         ! val_rad is in radians.

      if (unit_out == 'rad') factor2 = rad2rad
      if (unit_out == 'deg') factor2 = rad2deg

      wfat_convert_angle = factor2 * val_rad

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_dsortx(dx, inc, n, indx)
c     Experimental code to sort a double precision array and store its
c     original index. This algorithm is based on the subroutine isortx.
      Implicit NONE
#include "stdio.fh"
      Integer Inc, N
      double precision dx(*)
      Integer indx(*)
c -- sort array ix in ascending/decending order and set up index
c    array indx, such that indx(i)=original position of i
c      do 100 i=1,n
c      indx(i)=i
c 100   continue
C
      Integer I, L, Irra, J, IncX, Ir
      double precision dxrra
C
      if (n==1) return
      
      if (inc.eq.0)then
        write(6,105)inc,n
 105       format(' error isortx: inc, n ',2i10)
        stop 20
      else if(inc.lt.0)then
        incx=-inc
        do 110 i=0,n-1
        dx(i*incx+1)=-dx(i*incx+1)
 110       continue
      else
        incx=inc
      endif
c
      l=n/2+1
      ir=n
 10     continue
        if(l.gt.1) then
          l=l-1
          dxrra=dx((l-1)*incx+1)
          irra=indx(l)
        else
          dxrra=dx((ir-1)*incx+1)
          irra=indx(ir)
          dx((ir-1)*incx+1)=dx(1)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1) then
            dx(1)=dxrra
            indx(1)=irra
            goto 999
          end if
        end if
        i=l
        j=2*l
 20         if(j.le.ir) then
          if(j.lt.ir) then
            if(dx((j-1)*incx+1).lt.dx(j*incx+1)) j=j+1
          end if
          if(dxrra.lt.dx((j-1)*incx+1)) then
            dx((i-1)*incx+1)=dx((j-1)*incx+1)
            indx(i)=indx(j)
            i=j
            j=2*j
          else
            j=ir+1
          end if
        goto 20
        end if
        dx((i-1)*incx+1)=dxrra
        indx(i)=irra
      goto 10
c
 999   continue
      if(inc.lt.0)then
        do 210 i=0,n-1
        dx(i*incx+1)=-dx(i*incx+1)
 210       continue
      endif

cc      write(luout, '("max, min id ins = ", 2i7)')
cc     &     maxval(indx(1:n)), minval(indx(1:n))
      
      return
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_unique_count_i(nl, list, n)
c     DESCRIPTION:
c     Counts the number of different integers in a list.
      implicit none
      integer nl
      integer list(nl)
      integer n

      integer minv, maxv
      
      minv = minval(list) - 1
      maxv = maxval(list)
      n = 0
      do while (minv < maxv)
         n = n + 1
         minv = minval(list, mask=(list>minv))
      enddo
       
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dsortx_simple()

      implicit none


c
c      do i = 1, n
c         do j = 1, n
c            if (j /= k .and. d0 < d(j)) d0 = d(j)
c         enddo
c         d(i) = d0
c
c
c      enddo
c      


      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_mabyte_copy(m, s, i_s)
c     DESCRIPTION:
c     Copies the first m characters in string s into an ma_byte array.
c     The caller is responsible in ensuring that the m elements of the
c     ma_byte array that will be filled have been allocated.
      implicit none
#include "mafdecls.fh"
      integer m            ! [input] length of the string
      character(len=m) s   ! [input] the string
      integer i_s          ! [output] ma index if byte array in byte_mb

      integer i
      
      do i = 1, m
         byte_mb(i_s+i-1) = s(i:i)
      enddo
      
      end
c=======================================================
      
