c=======================================================
c=======================================================
      subroutine wfat_detect_degen(nb, ns, erg, thr, dgid, dgn, norb_dg)
c     DESCRIPTION:
c     dgid = The ID of orbitals taking degeneracy into account. This means that orbitals
c            belonging to the same degenerate subspace will have the same ID.
c     dgn = dgn(i,j) contains the dimension of the degenerate subspace where the i-th
c           orbital in the j-th spin channel belongs to. This means, for example, if in
c           the j-th spin channel, the 3rd, 5th, and 12th orbitals are degenerate, then
c           dgn(3,j) = dgn(5,j) = dgn(12,j).
c     norb_dg = The number of degenerate subspaces within each spin channel.
      
      implicit none
      integer nb
      integer ns
      double precision erg(nb, ns)
      double precision thr
      integer dgid(nb, ns)           ! [output]
      integer dgn(nb, ns)            ! [output]
      integer norb_dg(ns)            ! [output]

      logical isdtct(nb)
      integer is, id, i, j
      double precision derg
      

      do is = 1, ns
         id = 0
         isdtct = .false.
         do i = 1, nb

            if (.not. isdtct(i)) then
               id = id + 1
               dgid(i,is) = id
            
               do j = i+1, nb
                  derg = erg(i,is) - erg(j,is)
    
                  if ( abs(derg) <= thr ) then
                     dgid(j,is) = id
                     isdtct(j) = .true.
                  endif
               enddo
            endif
            
            isdtct(i) = .true.
         enddo
      enddo

      
      do is = 1, ns
         norb_dg(is) = maxval(dgid(:,is))
      enddo


      do is = 1, ns
         do i = 1, nb

            dgn(i,is) = 0
            do j = 1, nb
               if (dgid(j,is) == dgid(i,is)) dgn(i,is) = dgn(i,is) + 1
            enddo

         enddo
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_degen(nb, ns, erg, occ, thr, dgid, dgn)

      implicit none
#include "stdio.fh"
      integer nb
      integer ns
      double precision erg(nb,ns)
      double precision occ(nb,ns)
      double precision thr
      integer dgid(nb,ns)
      integer dgn(nb,ns)

      integer is, i, nel(ns)


      write(luout, '(5x, "Degenerate threshold = ", es15.6)') thr
      write(luout, '(5x, "Orbital degeneracies '//
     &     '(spin-a: 2nd-5th columns, spin-b: 6th-9th columns):")')
      

c     ==== Print column titles ====
      call wfat_print_repeat(luout, 5+ns*43, 5, '-')
      write(luout, '("")')
      write(luout, '(5x, a5)', advance='no') 'No.'
      do is = 1, ns
         nel(is) = int( sum(occ(:,is)) )
         write(luout, '(a23, a10, 2a5)', advance='no')
     &        'Energy (H)', 'Occ.', 'ID', 'Size'
      enddo
      write(luout, '(" ")')
      call wfat_print_repeat(luout, 5+ns*43, 5, '-')
      write(luout, '("")')

            
c     ==== Print ID's ====
      do i = 1, maxval(nel)
         write(luout, '(5x, i5)', advance='no') i
         do is = 1, ns
            if (i > nel(is)) then
               write(luout, '(a43)') ' '
            else
               write(luout, '(es23.12, f10.5, 2i5)', advance='no')
     &              erg(i,is), occ(i,is), dgid(i,is), dgn(i,is)
            endif
         enddo
         write(luout, '(" ")')
      enddo

      end
c=======================================================

      
c============================================================
c============================================================
      subroutine wfat_detect_degen0(nch, orbid, spinid, dgid0, dgn0)
c     dgid0 = The degeneracy ID of the ionized channels whose MO ID and spin channel are
c             specified in orbid and spinid, respectively. Note that ionized channels
c             having the same degenerate ID but belonging to different spin channels are
c             possible. In this case it is not possible to tell that these two channels
c             belong to different spin channels without consulting the corresponding
c             element of spinid. Therefore, when the spin channel of degenerate ionized
c             channels is needed, use spinid alongside dgid0 produced by this subroutine.
      implicit none
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
      integer nch          ! [input],
      integer orbid(nch)   ! [input],
      integer spinid(nch)  ! [input],
      integer dgid0(nch)   ! [output],
      integer dgn0(nch)    ! [output], the degree of degeneracy of the ionized channels.

      integer i, io, is

c     dgid0 is a subset of int_mb(k_modgid:k_modgid+nbas0*nspin0-1)
c     dgn0 is a subset of int_mb(k_dgn:k_dgn+nbas0*nspin0-1)
      do i = 1, nch
         io = orbid(i)
         is = spinid(i)
         dgid0(i) = int_mb(k_modgid + (is-1)*nbas0 + (io-1))
         dgn0(i) = int_mb(k_modgn + (is-1)*nbas0 + (io-1))
      enddo

      end
c============================================================

      
c============================================================
c============================================================
      subroutine wfat_get_degen_mfold(nch, orbid, spinid, dgid0, dgn0, 
     &     ich, mfid, orbid0)
c     Looks for the ionized channels that are within the same degenerate subspace
c     as channel ich. Therefore, given ich, :
c     1) mfid contains the ID's of the ionized channels that are
c        in the same degenerate manifold and has the same spin as ich (it
c        includes ich too) , and
c     2) orbid0 contains the ID's of the MO's that are in the same
c        degenerate manifold and has the same spin as the MO ich corresponds 
c        to (it includes the MO corresponding to ich too).
c     Hence, for a given ich, mfid and orbid0 constitute two equivalent
c     representations of the same set of degenerate MO's.

      implicit none
#include "errquit.fh"
      integer nch          ! [input],
      integer orbid(nch)   ! [input],
      integer spinid(nch)  ! [input],
      integer dgid0(nch)   ! [input], this array is produced by wfat_detect_degen0.
      integer dgn0(nch)    ! [input], this array is produced by wfat_detect_degen0.
      integer ich          ! [input], an integer whose value satisfies 1 <= ich <= nch.
      integer mfid(1:*)    ! [output], 1 <= mfid(i) <= nch, with i runs from 1 to dgn0(ich).
      integer orbid0(1:*)  ! [output], 1 <= orbid0(i) <= nbas0, with i runs from 1 to dgn0(ich).
      
      integer k, jch, j, j_, jo
      
      
      if (ich > nch) call errquit('WFAT: The argument ich of '//
     &     'wfat_get_degen_mfold must not be larger than nch', 34, UERR)
      
      k = 0
      do jch = 1, nch
         if (dgid0(jch) == dgid0(ich) .and. spinid(jch) == spinid(ich))
     &        then
            k = k + 1
            mfid(k) = jch
         endif
      enddo

      if (k /= dgn0(ich)) call errquit('WFAT: An inconsistency is '//
     &     'detected when identifying degenerate manifold. This is '//
     &     'most likely caused by not all degenerate MO''s within '//
     &     'a given degenerate manifold being included in the '//
     &     'ionized channels.', 23, UERR)

      do j = 1, dgn0(ich)
         j_ = mfid(j)
         jo = orbid(j_)
         orbid0(j) = jo
      enddo
      
      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_dgn_nsub_ch(ns, nch, orbid, spinid, nsub, maxn)
c     DESCRIPTION:
c     Computes the number of degenerate subspaces among the ionized channels whose IDs
c     are contained in orbid and whose spin channels are given in spinid. In doing so,
c     channels of different spins are considered as belonging to a different degenerate
c     subspace. It also returns the dimension of the largest subspace.

      implicit none
#include "stdio.fh"
      integer ns           ! [input]
      integer nch          ! [input]
      integer orbid(nch)   ! [input]
      integer spinid(nch)  ! [input]
      integer nsub         ! [output]
      integer maxn         ! [output]

      logical wfat_check_invec
      integer nsub_s(ns), id(nch), dgid0(nch), dgn0(nch), is, ich
      external wfat_check_invec
      

      call wfat_detect_degen0(nch, orbid, spinid, dgid0, dgn0)
      maxn = maxval(dgn0)
      
      do is = 1, ns
         id = -9999
         do ich = 1, nch
c           ==== Collect degenerate IDs belonging to spin channel is ====
            if (spinid(ich) == is) id(ich) = dgid0(ich)
         enddo
         call wfat_unique_count_i(nch, id, nsub_s(is))     ! 1)
         if (wfat_check_invec(-9999, nch, id))
     &        nsub_s(is) = nsub_s(is) - 1          ! 2)
      enddo
      nsub = sum(nsub_s)

      end
c     NOTES:
c     1) nsub_s(is) is the number of degenerate subspaces among the channels in orbid having
c        spin is.
c     2) The substraction by 1 is done because -9999 must be excluded in counting nsub_s(is).
c============================================================

      
c============================================================
c============================================================
      subroutine wfat_dgn_sub(nch, orbid, spinid, nsub, maxn, dgs, dgss,
     &     ndgs)

      implicit none
#include "stdio.fh"
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      integer nsub            ! [input]
      integer maxn            ! [input]
      integer dgs(maxn,nsub)  ! [output]
      integer dgss(nsub)      ! [output]
      integer ndgs(nsub)      ! [output]

      logical done(nch)
      integer i, ich, nn, j, dgid0(nch), dgn0(nch), mfid(nch)
      integer orbid0(nch)

      call wfat_detect_degen0(nch, orbid, spinid, dgid0, dgn0)

      dgs = -9999
      done = .false.
      do i = 1, nsub     ! Loop over degenerate subspaces among the ionized channels.

c        ==== Collect channel IDs that have not been ====
c        ====    taken into a degenerate subspace    ====
         do ich = 1, nch
            if (.not. done(ich)) then
               nn = dgn0(ich)
               call wfat_get_degen_mfold(nch, orbid, spinid, dgid0, 
     &              dgn0, ich, mfid, orbid0)
               dgs(1:nn,i) = mfid(1:nn)
               do j = 1, nn
                  done(mfid(j)) = .true.
               enddo
               exit      ! 1)
            endif
         enddo

c        ==== Spin channel of this degenerate subspace ====
         dgss(i) = spinid(dgs(1,i))

c        ==== Size of this degenerate subspace ====
         ndgs(i) = dgn0(dgs(1,i))
      enddo

      end
c     NOTES:
c     1) Exit otherwise the next set of unused channels ID will overwrite what's been
c        stored in dgs(:,i).
c============================================================
