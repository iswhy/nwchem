c=================================================
c=================================================
      subroutine wfat_init_progress(nptp, g_nptx, npt_tot)

      implicit none
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer nptp              ! [input]
      integer g_nptx            ! [output]
      double precision npt_tot  ! [output]

      integer g_npt_tot

      if (.not. ga_create(mt_dbl, 1, 1, 'nptx', 1, 1, g_nptx))
     &     call errquit('WFAT: Cannot create nptx global array.',
     &     1, GA_ERR)
      call ga_zero(g_nptx)
      if (.not. ga_create(mt_dbl, 1, 1, 'npt_tot', 1, 1, g_npt_tot))
     &     call errquit('WFAT: Cannot create npt_tot global array.',
     &     1, GA_ERR)
      call ga_zero(g_npt_tot)
      call ga_acc(g_npt_tot, 1, 1, 1, 1, dble(nptp), 1, 1.0d0)
      call ga_sync()
      call ga_get(g_npt_tot, 1, 1, 1, 1, npt_tot, 1)
      call ga_sync()
      if (.not. ga_destroy(g_npt_tot))
     &     call errquit('WFAT: Cannot destroy npt_tot.', 1, GA_ERR)

      end
c=================================================

      
c=================================================
c=================================================
      subroutine wfat_intgl_progress(k, pr_count, dt, npt_tot,
     &     g_nptx, nptx)

      implicit none
#include "global.fh"
#include "stdio.fh"
      integer k                 ! [input]
      integer pr_count          ! [input]
      double precision dt       ! [input]
      double precision npt_tot  ! [input], computed by wfat_init_progress
      integer g_nptx            ! [in/output]
      double precision nptx     ! [output]

      character*12 nptx_str, npt_tot_str, dt_str

      call ga_zero(g_nptx)
      call ga_acc(g_nptx, 1, 1, 1, 1, dble(pr_count), 1, 1.0d0)
      call ga_sync()
      if (ga_nodeid() == 0) then
         call ga_get(g_nptx, 1, 1, 1, 1, nptx, 1)
         write(nptx_str, '(i12)') int(nptx)
         write(npt_tot_str, '(i12)') int(npt_tot)
         write(dt_str, '(f12.2)') dt
         write(luout, '(7x, "Progress : ", i3, "% ", a, 3x, '//
     &        '"time needed = ", a, " s")')
     &        k*10, '('//trim(adjustl(nptx_str))//'/'
     &        //trim(adjustl(npt_tot_str))//')', trim(adjustl(dt_str))
      endif
      call ga_sync()

      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_terminate_progress(g_nptx)

      implicit none
#include "global.fh"
#include "errquit.fh"
      integer g_nptx    ! [input]

      if (.not. ga_destroy(g_nptx))
     &     call errquit('WFAT: Cannot destroy nptx.', 1, GA_ERR)

      end
c=================================================
